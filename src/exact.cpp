// -*- mode: c++; coding: utf-8 -*-
/*! @file exact.cpp
    @brief Inplementation of ExactModel class
*/
#include "exact.hpp"

#include <functional>

#include <json.hpp>

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/itertools.hpp>
#include <wtl/algorithm.hpp>
#include <wtl/math.hpp>

namespace likeligrid {

bool ExactModel::SIGINT_RAISED_ = false;

ExactModel::ExactModel(const std::string& infile, const size_t max_sites) {HERE;
    nlohmann::json jso;
    {
        wtl::izfstream ist(infile);
        ist >> jso;
    }
    names_ = jso["pathway"].get<std::vector<std::string>>();
    const size_t npath = names_.size();
    annot_.reserve(npath);
    for (const std::string& s: jso["annotation"]) {
        annot_.emplace_back(s);
    }
    std::cerr << annot_ << std::endl;

    const size_t nsam = jso["sample"].size();
    genot_.reserve(nsam);
    for (const std::string& s: jso["sample"]) {
        genot_.emplace_back(s);
    }
    // std::cerr << genot_ << std::endl;

    const size_t ngene = genot_[0].size();
    nsam_with_s_.assign(ngene, 0);
    std::valarray<double> s_gene(ngene);
    for (const auto& bits: genot_) {
        const size_t s = bits.count();
        ++nsam_with_s_[s];
        if (s > max_sites) continue;
        for (size_t j=0; j<ngene; ++j) {
            if (bits.test(j)) ++s_gene[j];
        }
    }
    wtl::rstrip(&nsam_with_s_);
    std::cerr << "Original N_s: " << nsam_with_s_ << std::endl;
    if (max_sites + 1 < nsam_with_s_.size()) {
        nsam_with_s_.resize(max_sites + 1);
        std::cerr << "Filtered N_s: " << nsam_with_s_ << std::endl;
    } else {
        std::cerr << "Note: -s is too large" << std::endl;
    }
    const auto final_max_s = nsam_with_s_.size() - 1;
    for (size_t s=2; s<=final_max_s; ++s) {
        lnp_const_ += nsam_with_s_[s] * std::log(wtl::factorial(s));
    }
    w_gene_ = s_gene / s_gene.sum();
    std::cerr << "s_gene : " << s_gene << std::endl;
    std::cerr << "w_gene_: " << w_gene_ << std::endl;
    for (size_t j=0; j<s_gene.size(); ++j) {
        if (s_gene[j] > 0) {
            lnp_const_ += s_gene[j] * std::log(w_gene_[j]);
        }
    }
    std::cerr << "lnp_const_: " << lnp_const_ << std::endl;

    // TODO
    a_pathway_.resize(npath);
    for (size_t j=0; j<npath; ++j) {
        for (size_t i=0; i<nsam; ++i) {
            size_t s = (genot_[i] & annot_[j]).count();
            if (s > 0) {
                a_pathway_[j] += --s;
            }
        }
    }
    std::cerr << "a_pathway_: " << a_pathway_ << std::endl;
}

void ExactModel::run() {HERE;
    calc_loglik({0.5, 0.5, 0.5, 0.5});
}

double ExactModel::calc_loglik(const std::valarray<double>& params) const {
    double loglik = (a_pathway_ * std::log(params)).sum();
    calc_denom(params, 2);
    return loglik += lnp_const_;
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
        denoms_(max_sites_ + 1, 0.0) {

        mutate(bits_t(w_gene.size(), 0), bits_t(annot.size(), 0), 1.0);
    }
    const std::vector<double>& get() const {return denoms_;}

  private:
    void mutate(const bits_t& genotype, const bits_t& pathtype, const double anc_p) {
        const size_t s = genotype.count() + 1;
        for (bits_t mut_gene(genotype.size(), 1); mut_gene.any(); mut_gene <<= 1) {
            if ((genotype & mut_gene).any()) continue;
            const bits_t mut_path = translate(mut_gene);
            double p = anc_p;
            p *= w_gene_[mut_gene.find_first()];
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
            size_t j = 0;
            while ((j = recurrent.find_next(j)) < bits_t::npos) {
                p *= th_path_[j];
            }
            return p;
        } else {return 1.0;}
    }

    // TODO: memoize
    bits_t translate(const bits_t& mut_gene) const {
        bits_t mut_path(annot_.size(), 0);
        for (size_t j=0; j<annot_.size(); ++j) {
            mut_path.set(j, (annot_[j] & mut_gene).any());
        }
        return mut_path;
    }

    const std::valarray<double>& w_gene_;
    const std::valarray<double>& th_path_;
    const std::vector<bits_t>& annot_;
    const size_t max_sites_;
    std::vector<double> denoms_;
};

double ExactModel::calc_denom(const std::valarray<double>& th_path, const size_t s) const {
    std::cout << th_path << std::endl;
    double sum_prob = 0.0;
    std::cout << "D_s: " << Denoms(w_gene_, th_path, annot_, s).get() << std::endl;
    return sum_prob;
}

void ExactModel::unit_test() {HERE;
    ExactModel model("test.json.gz", 2);
    model.run();
}

} // namespace likeligrid
