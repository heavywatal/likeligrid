// -*- mode: c++; coding: utf-8 -*-
/*! @file genotype.cpp
    @brief Implementation of GenotypeModel class
*/
#include "genotype.hpp"

#include <json.hpp>

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/algorithm.hpp>

namespace likeligrid {

GenotypeModel::GenotypeModel(std::istream& ist, const size_t max_sites) {HERE;
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
    const size_t ngene = jso["sample"].at(0).get<std::string>().size();
    nsam_with_s_.assign(ngene + 1, 0);  // at most
    std::valarray<double> s_gene(ngene);
    for (const auto& bits: all_genotypes) {
        const size_t s = bits.count();
        ++nsam_with_s_[s];
        if (s > max_sites) continue;
        genot_.push_back(bits);
        for (size_t j=0; j<ngene; ++j) {
            if (bits[j]) ++s_gene[j];
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

    max_sites_ = nsam_with_s_.size() - 1;
    denoms_.resize(max_sites_ + 1);
    effects_.reserve(ngene);
    for (size_t j=0; j<ngene; ++j) {
        effects_.emplace_back(translate(j));
    }
    // std::cerr << "effects_: " << effects_ << std::endl;
}

void GenotypeModel::unit_test() {HERE;
    std::stringstream sst;
    sst <<
R"({
  "pathway": ["A", "B"],
  "annotation": ["0011", "1100"],
  "sample": ["0011", "0101", "1001", "0110", "1010", "1100"]
})";
    GenotypeModel model(sst, 4);
    std::cerr << model.calc_loglik({1.0, 1.0}) << std::endl;
}

} // namespace likeligrid
