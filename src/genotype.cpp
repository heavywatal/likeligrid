/*! @file genotype.cpp
    @brief Implementation of GenotypeModel class
*/
#include "genotype.hpp"
#include "util.hpp"

#include <wtl/debug.hpp>
#include <wtl/resource.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zlib.hpp>
#include <wtl/algorithm.hpp>
#include <wtl/math.hpp>
#include <wtl/exception.hpp>

#include <clippson/json.hpp>

namespace likeligrid {

GenotypeModel::GenotypeModel(const std::string& infile, const size_t max_sites)
: filename_(infile) {
    HERE;
    wtl::zlib::ifstream ist(filename_);
    init(ist, max_sites);
}

void GenotypeModel::init(std::istream& ist, const size_t max_sites) {HERE;
    nlohmann::json jso;
    ist >> jso;
    names_ = jso["pathway"].get<std::vector<std::string>>();
    num_pathways_ = names_.size();
    annot_.reserve(num_pathways_);
    for (const std::string& s: jso.value("annotation", std::vector<std::string>{})) {
        annot_.emplace_back(s);
    }
    std::cerr << "annot_: " << annot_ << std::endl;

    const size_t nsam = jso["sample"].size();
    std::vector<bits_t> all_genotypes;
    all_genotypes.reserve(nsam);
    for (const std::string& s: jso.value("sample", std::vector<std::string>{})) {
        all_genotypes.emplace_back(s);
    }

    genot_.reserve(nsam);
    num_genes_ = jso["sample"].at(0u).get<std::string>().size();
    nsam_with_s_.assign(num_genes_ + 1u, 0u);  // at most
    std::valarray<double> s_gene(num_genes_);
    for (const auto& bits: all_genotypes) {
        const size_t s = bits.count();
        ++nsam_with_s_[s];
        if (s > max_sites) continue;
        genot_.push_back(bits);
        for (size_t j=0u; j<num_genes_; ++j) {
            if (bits[j]) ++s_gene[j];
        }
    }
    wtl::rstrip(&nsam_with_s_);
    std::cerr << "Original N_s: " << nsam_with_s_ << std::endl;
    if (max_sites + 1u < nsam_with_s_.size()) {
        nsam_with_s_.resize(max_sites + 1u);
        std::cerr << "Using N_s: " << nsam_with_s_ << std::endl;
    } else {
        std::cerr << "Note: -s is too large" << std::endl;
    }
    const std::valarray<double> w_gene = s_gene / s_gene.sum();
    ln_w_gene_ = std::log(w_gene);
    std::cerr << "s_gene: " << s_gene << std::endl;
    std::cerr << "w_gene: " << w_gene << std::endl;
    std::cerr << "ln_w_gene_: " << ln_w_gene_ << std::endl;

    max_sites_ = nsam_with_s_.size() - 1u;
    effects_.reserve(num_genes_);
    for (size_t j=0u; j<num_genes_; ++j) {
        effects_.emplace_back(translate(j));
    }
    // std::cerr << "effects_: " << effects_ << std::endl;
}

bool GenotypeModel::set_epistasis(const std::pair<size_t, size_t>& pair, const bool pleiotropy) {HERE;
    if (pair.first == pair.second) return false;
    epistasis_pair_ = pair;
    pleiotropy_idx_ = epistasis_idx_ = num_pathways_;
    std::ostringstream oss;
    oss << names_.at(pair.first) << ":" << names_.at(pair.second);
    names_.push_back(oss.str());
    std::cerr << "epistasis: " << names_.back() << std::endl;
    if (pleiotropy) {
        names_.push_back("pleiotropy");
        ++pleiotropy_idx_;
        std::cerr << "pleiotropy: true" << std::endl;
    }
    return epistasis_ = true;
}

double GenotypeModel::calc_loglik(const std::valarray<double>& theta) {
    ln_theta_ = std::log(theta);
    // std::cerr << "denoms_: " << denoms_ << std::endl;
    double loglik = 0.0;
    for (const auto& genotype: genot_) {
        loglik += lnp_sample(genotype);
    }
    ln_denoms_.resize(max_sites_ + 1u);
    ln_denoms_ = -std::numeric_limits<double>::infinity();
    mutate();
    // std::cerr << "lnD: " << ln_denoms_ << std::endl;
    // -inf, 0, D2, D3, ...
    for (size_t s=2u; s<=max_sites_; ++s) {
        loglik -= nsam_with_s_[s] * ln_denoms_[s];
    }
    return loglik;
}

inline std::valarray<size_t> to_indices(const bits_t& bits) {
    std::valarray<size_t> indices(bits.count());
    for (size_t i=0u, j=0u; i<indices.size(); ++j) {
        if (bits[j]) {
            indices[i] = j;
            ++i;
        }
    }
    return indices;
}

inline double slice_sum(const std::valarray<double>& ln_coefs, const bits_t& bits) {
    double lnp = 0.0;
    for (size_t i=0; i<ln_coefs.size(); ++i) {
        if (bits[i]) lnp += ln_coefs[i];
    }
    return lnp;
}

inline double add_lnp(const double ln_bigger, const double ln_smaller) {
    return ln_bigger + std::log1p(std::exp(ln_smaller - ln_bigger));
}

inline double sub_lnp(const double ln_bigger, const double ln_smaller) {
    return ln_bigger + std::log1p(-std::exp(ln_smaller - ln_bigger));
}

double GenotypeModel::lnp_sample(const bits_t& genotype) const {
    double lnp = -std::numeric_limits<double>::infinity();
    const double lnp_basic = slice_sum(ln_w_gene_, genotype);
    auto mut_route = to_indices(genotype);
    do {
        lnp = add_lnp(lnp_basic + sum_ln_theta(mut_route), lnp);
    } while (std::next_permutation(std::begin(mut_route), std::end(mut_route)));
    return lnp;
}

void GenotypeModel::mutate(const bits_t& genotype, const bits_t& pathtype, const double anc_lnp, const double open_lnp) {
    const auto s = genotype.count() + 1u;
    for (size_t j=0u; j<num_genes_; ++j) {
        if (genotype[j]) continue;
        if (ln_w_gene_[j] == -std::numeric_limits<double>::infinity()) continue;
        const bits_t& mut_path = effects_[j];
        double lnp = anc_lnp;
        lnp += ln_w_gene_[j];
        lnp -= open_lnp;
        lnp += ln_theta_if_subset(pathtype, mut_path);
        if (epistasis_) {lnp += ln_theta_if_paired(pathtype, mut_path);}
        ln_denoms_[s] = add_lnp(lnp, ln_denoms_[s]);
        if (s < max_sites_) {
            if (wtl::SIGINT_RAISED()) {throw wtl::KeyboardInterrupt();}
            mutate(bits_t(genotype).set(j), pathtype | mut_path, lnp, sub_lnp(open_lnp, ln_w_gene_[j]));
        }
    }
}

void GenotypeModel::benchmark(const int n) {
    const std::valarray<double> param(0.9, num_pathways_);
    double leaves = wtl::pow(static_cast<double>(num_genes_), static_cast<unsigned int>(max_sites_));
    std::cerr << "# parameters: " << num_pathways_ << std::endl;
    std::cerr << "width: " << num_genes_ << std::endl;
    std::cerr << "depth: " << max_sites_ << std::endl;
    std::cerr << "w ^ d: " << leaves * 1e-6 << " M" <<std::endl;
    std::cerr << wtl::delta_rusage([&param,this]() {calc_loglik(param);}, n);
}

} // namespace likeligrid
