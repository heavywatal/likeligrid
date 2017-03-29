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

double ExactModel::calc_denom(const Eigen::ArrayXd& params, const size_t s) const {
    std::cout << params << std::endl;
    double sum_prob = 0.0;
    boost::dynamic_bitset<> flags(s);
    flags.set();
    flags.resize(w_gene_.size(), 0);
    std::cout << flags << std::endl;
    return sum_prob;
}


template <class value_type>
class Bits final: public wtl::itertools::Generator<value_type> {
  public:
    using typename wtl::itertools::Generator<value_type>::coro_t;
    using typename wtl::itertools::Generator<value_type>::size_type;
    using typename wtl::itertools::Generator<value_type>::value_size_t;
    Bits() = delete;
    explicit Bits(const size_t length, const size_t sum):
        wtl::itertools::Generator<value_type>(),
        sum_(sum),
        value_(length, 0) {}
    ~Bits() = default;

    void reset() {
        value_.reset();
        s_ = 1;
        this->cnt_ = 0;
    }
    virtual size_type max_count() const override {
        return wtl::permut(value_.size(), sum_);
    }

  private:
    virtual void source(typename coro_t::push_type& yield, const size_type skip) override {
        value_type anchor(value_);
        for (value_type current(value_.size(), 1); current.any(); current <<= 1) {
            if ((value_ & current).any()) continue;
            value_ |= current;
            if (s_ < sum_) {
                ++s_;
                yield(value_type(value_)); // intermediate genotype
                source(yield, skip);
            } else {
                if (++this->cnt_ > skip) {
                    yield(value_type(value_));  // final genotype
                }
            }
            value_ = anchor;
        }
        --s_;
    }
    const value_size_t sum_;
    value_type value_;
    value_size_t s_ = 1;
};

void ExactModel::unit_test() {HERE;
    ExactModel model("test.json.gz", 3);
    model.run();
    Bits<boost::dynamic_bitset<>> gen(5, 2);
    std::cout << "max_count: " << gen.max_count() << std::endl;
    for (const auto x: gen()) {
        std::cout << x << " " << gen.count() << std::endl;
    }
}

} // namespace likeligrid
