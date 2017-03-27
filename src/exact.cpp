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

    for (const auto s: raw_s_sample) {
        if (s <= final_max_s) {
            lnp_const_ += std::log(wtl::factorial(s));
        }
    }

    genotypes_ = wtl::eigen::ArrayXX<uint>(jso["sample"]);
    genotypes_ = wtl::eigen::filter(genotypes_, wtl::eigen::ArrayX(raw_s_sample) <= final_max_s);
    const Eigen::ArrayXd s_gene = genotypes_.colwise().sum().cast<double>();
    std::cerr << "s_gene : " << s_gene.transpose() << std::endl;
    for (Eigen::Index i=0; i<s_gene.size(); ++i) {
        if (s_gene[i] > 0) {
            lnp_const_ += s_gene[i] * std::log(s_gene[i]);
        }
    }
    const auto s_total = s_gene.sum();
    lnp_const_ -= s_total * std::log(s_total);
    std::cerr << "lnp_const_: " << lnp_const_ << std::endl;
}

void ExactModel::run() {HERE;
}

double ExactModel::calc_loglik(const Eigen::ArrayXd& params) const {
    return lnp_const_;
}

void ExactModel::unit_test() {HERE;
    ExactModel model("test.json.gz", 3);
    model.run();
}

} // namespace likeligrid
