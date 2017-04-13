// -*- mode: c++; coding: utf-8 -*-
/*! @file exact.cpp
    @brief Inplementation of ExactModel class
*/
#include "exact.hpp"
#include "util.hpp"

#include <functional>

#include <json.hpp>

#include <wtl/debug.hpp>
#include <wtl/exception.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/algorithm.hpp>
#include <wtl/numeric.hpp>
#include <wtl/math.hpp>
#include <wtl/os.hpp>

namespace likeligrid {

const std::vector<double> ExactModel::STEPS_ = {0.4, 0.2, 0.1, 0.05, 0.02, 0.01};
const std::vector<size_t> ExactModel::BREAKS_ = {5, 5, 5, 5, 6, 5};
bool ExactModel::SIGINT_RAISED_ = false;

ExactModel::ExactModel(const std::string& infile, const size_t max_sites):
    ExactModel(wtl::izfstream(infile), max_sites) {HERE;}

ExactModel::ExactModel(std::istream&& ist, const size_t max_sites) {HERE;
    nlohmann::json jso;
    ist >> jso;
    names_ = jso["pathway"].get<std::vector<std::string>>();
    const size_t npath = names_.size();
    annot_.reserve(npath);
    for (const std::string& s: jso["annotation"]) {
        annot_.emplace_back(s);
    }
    std::cerr << "annot_: " << annot_ << std::endl;

    const size_t nsam = jso["sample"].size();
    std::vector<bits_t> all_genotypes;
    all_genotypes.reserve(nsam);
    for (const std::string& s: jso["sample"]) {
        all_genotypes.emplace_back(s);
    }

    genot_.reserve(nsam);
    const size_t ngene = all_genotypes[0].size();
    nsam_with_s_.assign(ngene + 1, 0);  // at most
    std::valarray<double> s_gene(ngene);
    for (const auto& bits: all_genotypes) {
        const size_t s = bits.count();
        ++nsam_with_s_[s];
        if (s > max_sites) continue;
        genot_.push_back(bits);
        bits_t::size_type j = bits.find_first();
        while (j < bits_t::npos) {
            ++s_gene[j];
            j = bits.find_next(j);
        }
    }
    wtl::rstrip(&nsam_with_s_);
    std::cerr << "Original N_s: " << nsam_with_s_ << std::endl;
    if (max_sites + 1 < nsam_with_s_.size()) {
        nsam_with_s_.resize(max_sites + 1);
        std::cerr << "Using N_s: " << nsam_with_s_ << std::endl;
    } else {
        std::cerr << "Note: -s is too large" << std::endl;
    }
    w_gene_ = s_gene / s_gene.sum();
    std::cerr << "s_gene : " << s_gene << std::endl;
    std::cerr << "w_gene_: " << w_gene_ << std::endl;

    mle_params_.resize(annot_.size());
    mle_params_ = 1.2;
}

void ExactModel::run() {HERE;
    const std::string outfile = init_meta();
    std::cerr << "mle_params_: " << mle_params_ << std::endl;
    if (outfile.empty()) return;
    const auto axes = make_vicinity(mle_params_, BREAKS_.at(stage_), 2.0 * STEPS_.at(stage_));
    for (size_t j=0; j<names_.size(); ++j) {
        std::cerr << names_[j] << ": " << axes[j] << std::endl;
    }
    {
        wtl::ozfstream fout(outfile, std::ios::out | std::ios::app);
        std::cerr << "Writing: " << outfile << std::endl;
        run_impl(fout, wtl::itertools::product(axes));
    }
    run();
}

void ExactModel::go() {HERE;
    std::cerr << "mle_params_: " << mle_params_ << std::endl;
    if (stage_ >= STEPS_.size()) return;
    const auto axes = make_vicinity(mle_params_, BREAKS_.at(stage_), 2.0 * STEPS_.at(stage_));
    for (size_t j=0; j<names_.size(); ++j) {
        std::cerr << names_[j] << ": " << axes[j] << std::endl;
    }
    {
        std::stringstream sst;
        run_impl(sst, wtl::itertools::product(axes));
        std::cout << sst.str();
        read_results(sst);
    }
    ++stage_;
    go();
}

void ExactModel::run_impl(std::ostream& ost, wtl::itertools::Generator<std::valarray<double>>&& gen) const {HERE;
    std::cerr << skip_ << " to " << gen.max_count() << std::endl;
    if (skip_ == 0) {
        ost << "##max_count=" << gen.max_count() << "\n";
        ost << "##max_sites=" << nsam_with_s_.size() - 1 << "\n";
        ost << "##step=" << STEPS_.at(stage_) << "\n";
        ost << "loglik\t" << wtl::join(names_, "\t") << "\n";
    }
    auto buffer = wtl::make_oss();
    for (const auto& th_path: gen(skip_)) {
        buffer << calc_loglik(th_path) << "\t"
               << wtl::str_join(th_path, "\t") << "\n";
        if (gen.count() % 100 == 0) {  // snapshot for long run
            std::cerr << "*" << std::flush;
            ost << buffer.str();
            ost.flush();
            buffer.str("");
        }
        if (SIGINT_RAISED_) {throw wtl::KeyboardInterrupt();}
    }
    std::cerr << "-\n";
    ost << buffer.str();
}

std::vector<bits_t> breakdown(const bits_t& bits) {
    std::vector<bits_t> singles;
    bits_t::size_type length = bits.size();
    singles.reserve(length);
    bits_t::size_type i = bits.find_first();
    while (i < bits_t::npos) {
        bits_t x(length);
        x.set(i);
        singles.push_back(x);
        i = bits.find_next(i);
    }
    return singles;
}

class Denoms {
  public:
    Denoms() = delete;
    Denoms(const std::valarray<double>& w_gene,
        const std::valarray<double>& th_path,
        const std::vector<bits_t>& annot,
        const size_t max_sites):
        w_gene_(w_gene),
        th_path_(th_path),
        annot_(annot),
        max_sites_(max_sites),
        denoms_(max_sites + 1)
    {
        const size_t ngene = w_gene.size();
        effects_.reserve(ngene);
        for (bits_t mut_gene(ngene, 1); mut_gene.any(); mut_gene <<= 1) {
            effects_.emplace_back(translate(mut_gene));
        }
        // std::cerr << "effects_: " << effects_ << std::endl;
        mutate(bits_t(ngene, 0), bits_t(annot.size(), 0), 1.0);
        // std::cerr << "denoms_: " << denoms_ << std::endl;
    }
    const std::valarray<double> log() const {return std::log(denoms_);}

    double lnp_sample(const bits_t& genotype) const {
        double p = 0.0;
        const double p_basic = prod_w(genotype);
        auto mutations = breakdown(genotype);
        do {
            p += p_basic * discount(mutations);
        } while (std::next_permutation(mutations.begin(), mutations.end()));
        return std::log(p);
    }

  private:
    void mutate(const bits_t& genotype, const bits_t& pathtype, const double anc_p) {
        const size_t s = genotype.count() + 1;
        for (bits_t mut_gene(genotype.size(), 1); mut_gene.any(); mut_gene <<= 1) {
            if ((genotype & mut_gene).any()) continue;
            const size_t pos = mut_gene.find_first();
            const bits_t& mut_path = effects_[pos];
            double p = anc_p;
            p *= w_gene_[pos];
            p *= discount(pathtype, mut_path);
            // std::cout << (genotype | mut_gene) << " " << (pathtype | mut_path) << " " << p << std::endl;
            denoms_[s] += p;
            if (s < max_sites_) {
                mutate(genotype | mut_gene, pathtype | mut_path, p);
            }
        }
    }

    double discount(const bits_t& pathtype, const bits_t& mut_path) const {
        if (mut_path.is_subset_of(pathtype)) {
            double p = 1.0;
            const bits_t recurrent = mut_path & pathtype;
            bits_t::size_type j = recurrent.find_first();
            while (j < bits_t::npos) {
                p *= th_path_[j];
                j = recurrent.find_next(j);
            }
            return p;
        } else {return 1.0;}
    }

    double discount(const std::vector<bits_t>& mutations) const {
        double p = 1.0;
        const auto npath = annot_.size();
        bits_t pathtype(npath, 0);
        for (const auto& mut_gene: mutations) {
            const auto mut_path = translate(mut_gene);
            p *= discount(pathtype, mut_path);
            pathtype |= mut_path;
        }
        return p;
    }

    bits_t translate(const bits_t& mut_gene) const {
        bits_t mut_path(annot_.size(), 0);
        for (size_t j=0; j<annot_.size(); ++j) {
            mut_path.set(j, (annot_[j] & mut_gene).any());
        }
        return mut_path;
    }

    double prod_w(const bits_t& genotype) const {
        double p = 1.0;
        bits_t::size_type j = genotype.find_first();
        while (j < bits_t::npos) {
            p *= w_gene_[j];
            j = genotype.find_next(j);
        }
        return p;
    }
    const std::valarray<double>& w_gene_;
    const std::valarray<double>& th_path_;
    const std::vector<bits_t>& annot_;
    const size_t max_sites_;
    std::valarray<double> denoms_;
    std::vector<bits_t> effects_;
};

double ExactModel::calc_loglik(const std::valarray<double>& th_path) const {
    const size_t max_sites = nsam_with_s_.size() - 1;
    Denoms subcalc(w_gene_, th_path, annot_, max_sites);
    double loglik = 0.0;
    for (const auto& genotype: genot_) {
        loglik += subcalc.lnp_sample(genotype);
    }
    const auto lnD = subcalc.log();
    // std::cout << "lnD: " << lnD << std::endl;
    // -inf, 0, D2, D3, ...
    for (size_t s=2; s<=max_sites; ++s) {
        loglik -= nsam_with_s_[s] * lnD[s];
    }
    return loglik;
}

std::string ExactModel::init_meta() {HERE;
    if (stage_ >= STEPS_.size()) return "";
    auto oss = wtl::make_oss(2, std::ios::fixed);
    oss << "grid-" << STEPS_.at(stage_) << ".tsv.gz";
    std::string outfile = oss.str();
    try {
        wtl::izfstream ist(outfile);
        std::cerr << "Reading: " << outfile << std::endl;
        read_results(ist);
        if (skip_ == 0) {
            ++stage_;
            outfile = init_meta();
        }
    } catch (std::ios::failure& e) {
        if (errno != 2) throw;
    }
    return outfile;
}

void ExactModel::read_results(std::istream& ist) {HERE;
    size_t max_count;
    double step;
    std::tie(max_count, std::ignore, step) = read_metadata(ist);
    stage_ = guess_stage(STEPS_, step);
    std::vector<std::string> colnames;
    std::valarray<double> mle_params;
    std::tie(skip_, colnames, mle_params) = read_body(ist);
    if (skip_ == max_count) {  // is complete file
        skip_ = 0;
        mle_params_.swap(mle_params);
    }
    if (names_ != colnames) {
        std::ostringstream oss;
        oss << "Contradiction in column names:\n"
            << "genotype file: " << names_ << "\n"
            << "result file:" << colnames;
        throw std::runtime_error(oss.str());
    }
}

void ExactModel::unit_test() {HERE;
    std::stringstream sst;
    sst <<
R"({
  "pathway": ["A", "B"],
  "annotation": ["0011", "1100"],
  "sample": ["0011", "0101", "1001", "0110", "1010", "1100"]
})";
    ExactModel model(std::move(sst), 2);
    model.go();
}

} // namespace likeligrid
