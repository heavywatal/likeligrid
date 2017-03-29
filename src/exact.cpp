// -*- mode: c++; coding: utf-8 -*-
/*! @file exact.cpp
    @brief Inplementation of ExactModel class
*/
#include "exact.hpp"

#include <functional>

#include <json.hpp>

#include <wtl/debug.hpp>
#include <wtl/exception.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/eigen.hpp>

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

    ArrayXXu pathtypes(nsam, npath);
    for (size_t i=0; i<nsam; ++i) {
        std::vector<uint> row;
        row.reserve(npath);
        const auto& g = genot_[i];
        for (const auto& a: annot_) {
            row.push_back((g & a).count());
        }
        pathtypes.row(i) = wtl::eigen::ArrayX(row);
    }
    a_pathway_ = pathtypes.unaryExpr([](uint x){
        if (x > 0) {return --x;} else {return x;}
    }).colwise().sum().cast<double>();
    std::cerr << "s_pathway : " << pathtypes.colwise().sum() << std::endl;
    std::cerr << "a_pathway_: " << a_pathway_.transpose() << std::endl;

    std::vector<size_t> raw_s_sample;
    raw_s_sample.reserve(nsam);
    const size_t ngene = genot_[0].size();
    nsam_with_s_.assign(ngene, 0);
    for (const auto& bits: genot_) {
        const size_t s = bits.count();
        raw_s_sample.emplace_back(s);
        ++nsam_with_s_[s];
    }

    std::cerr << "Original N_s: " << nsam_with_s_ << std::endl;
    if (max_sites + 1 < nsam_with_s_.size()) {
        nsam_with_s_.resize(max_sites + 1);
        std::cerr << "Filtered N_s: " << nsam_with_s_ << std::endl;
    } else {
        std::cerr << "Note: -s is too large" << std::endl;
    }
    while (nsam_with_s_.back() == 0) {
        nsam_with_s_.resize(nsam_with_s_.size() - 1);
    }
    const auto final_max_s = nsam_with_s_.size() - 1;
    for (size_t s=2; s<=final_max_s; ++s) {
        lnp_const_ += nsam_with_s_[s] * std::log(wtl::factorial(s));
    }

    ArrayXXu genotypes = wtl::eigen::ArrayXX<uint>(jso["sample"]);
    genotypes = wtl::eigen::filter(genotypes, wtl::eigen::ArrayX(raw_s_sample) <= final_max_s);
    const Eigen::ArrayXd s_gene = genotypes.colwise().sum().cast<double>();
    w_gene_ = s_gene / s_gene.sum();
    std::cerr << "s_gene : " << s_gene.transpose() << std::endl;
    std::cerr << "w_gene_: " << w_gene_.transpose() << std::endl;
    for (Eigen::Index i=0; i<s_gene.size(); ++i) {
        if (s_gene[i] > 0) {
            lnp_const_ += s_gene[i] * std::log(w_gene_[i]);
        }
    }
    std::cerr << "lnp_const_: " << lnp_const_ << std::endl;
}

void ExactModel::run() {HERE;
    calc_loglik(a_pathway_ / a_pathway_.sum());
}

double ExactModel::calc_loglik(const Eigen::ArrayXd& params) const {
    double loglik = (a_pathway_ * params.log()).sum();
    calc_denom(params, 2);
    return loglik += lnp_const_;
}

class Denom {
  public:
    typedef boost::dynamic_bitset<> bits_t;
    Denom() = delete;
    Denom(const std::valarray<double>& w_genes,
      const std::valarray<double>& params,
      const std::vector<bits_t>& annot,
      const size_t nsites):
        w_genes_(w_genes),
        params_(params),
        annot_(annot),
        genotype_(w_genes.size(), 0),
        pathtype_(annot.size(), 0),
        nsites_(nsites),
        denoms_(nsites + 1, 0.0) {}
    void mutate() {
        const bits_t tmp_genotype = genotype_;
        const bits_t tmp_pathtype = pathtype_;
        const double tmp_p = p_;
        ++s_;
        for (bits_t new_mut(genotype_.size(), 1); new_mut.any(); new_mut <<= 1) {
            if ((genotype_ & new_mut).any()) continue;
            calculate(new_mut);
            if (s_ < nsites_) {
                mutate();
            }
            genotype_ = tmp_genotype;
            pathtype_ = tmp_pathtype;
            p_ = tmp_p;
        }
        --s_;
    }
    void calculate(const bits_t new_mut) {
        p_ *= w_genes_[new_mut.find_first()];
        bits_t mut_path(annot_.size(), 0);
        for (size_t i=0; i<annot_.size(); ++i) {
            if ((annot_[i] & new_mut).any()) {
               if (pathtype_.test_set(i)) {
                   p_ *= params_[i];
               }
            }
        }
        denoms_[s_] += p_;
        genotype_ |= new_mut;
        std::cout << genotype_ << " " << pathtype_ << " "
                  << s_ << " " << p_ << std::endl;
    }
    const std::vector<double>& get() const {return denoms_;}
  private:
    std::valarray<double> w_genes_;
    std::valarray<double> params_;
    std::vector<bits_t> annot_;
    bits_t genotype_;
    bits_t pathtype_;
    const size_t nsites_;
    size_t s_ = 0;
    double p_ = 1.0;
    std::vector<double> denoms_;
};

double ExactModel::calc_denom(const Eigen::ArrayXd& params, const size_t s) const {
    std::cout << params << std::endl;
    double sum_prob = 0.0;
    boost::dynamic_bitset<> flags(s);
    flags.set();
    flags.resize(w_gene_.size(), 0);
    std::cout << flags << std::endl;

    // Denom d({0.2, 0.2, 0.3, 0.3}, annot_, 3);
    Denom d(wtl::eigen::valarray(w_gene_), {0.5, 1.0, 1.0, 1.0}, annot_, s);
    d.mutate();
    std::cout << "D_s: " << d.get() << std::endl;
    return sum_prob;
}

void ExactModel::unit_test() {HERE;
    ExactModel model("test.json.gz", 2);
    model.run();
}

} // namespace likeligrid
