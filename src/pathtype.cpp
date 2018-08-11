/*! @file pathtype.cpp
    @brief Implementation of PathtypeModel class
*/
#include "pathtype.hpp"

#include <wtl/debug.hpp>
#include <wtl/exception.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zlib.hpp>
#include <wtl/itertools.hpp>
#include <wtl/numeric.hpp>
#include <wtl/math.hpp>

#include <functional>
#include <stdexcept>
#include <cstdint>

namespace likeligrid {

PathtypeModel::PathtypeModel(const std::string& infile, const size_t max_sites):
    PathtypeModel(wtl::zlib::ifstream(infile), max_sites) {HERE;}

PathtypeModel::PathtypeModel(std::istream&& ist, const size_t max_sites) {HERE;
    names_ = wtl::read_header(ist);
    auto pathtypes = wtl::read_valarrays<uint_fast32_t>(ist);
    const auto raw_s_sample = wtl::row_sums(pathtypes);
    nsam_with_s_.assign(raw_s_sample.max() + 1u, 0u);
    for (const auto s: raw_s_sample) {
        ++nsam_with_s_[s];
    }
    std::cerr << "Original N_s: " << nsam_with_s_ << std::endl;
    if (max_sites + 1u < nsam_with_s_.size()) {
        nsam_with_s_.resize(max_sites + 1u);
        std::cerr << "Filtered N_s: " << nsam_with_s_ << std::endl;
    } else {
        std::cerr << "Note: -s is too large" << std::endl;
    }
    while (nsam_with_s_.back() == 0u) {
        nsam_with_s_.pop_back();
    }
    const auto final_max_s = static_cast<uint_fast32_t>(nsam_with_s_.size()) - 1u;
    pathtypes = wtl::filter(pathtypes, raw_s_sample <= final_max_s);

    const auto s_pathway = wtl::cast<double>(wtl::col_sums(pathtypes));
    w_pathway_ = s_pathway / s_pathway.sum();
    auto duptypes = pathtypes;
    for (auto& row: duptypes) {
        for (auto& x: row) {
            if (x > 0u) --x;
        }
    }
    a_pathway_ = wtl::cast<double>(wtl::col_sums(duptypes));
    for (size_t i=0u; i<pathtypes.size(); ++i) {
        lnp_const_ += std::log(wtl::multinomial(pathtypes[i]));
    }
    lnp_const_ += (s_pathway * std::log(w_pathway_)).sum();
    std::cerr << "s_pathway_: " << s_pathway << std::endl;
    std::cerr << "w_pathway_: " << w_pathway_ << std::endl;
    std::cerr << "a_pathway_: " << a_pathway_ << std::endl;
    std::cerr << "lnp_const_: " << lnp_const_ << std::endl;
    WTL_ASSERT(!std::isnan(lnp_const_));

    index_axes_.reserve(nsam_with_s_.size());
    std::vector<size_t> indices(a_pathway_.size());
    std::iota(std::begin(indices), std::end(indices), 0u);
    for (size_t s=0u; s<=nsam_with_s_.size(); ++s) {
        index_axes_.emplace_back(s, indices);
    }
}

double PathtypeModel::calc_loglik(const std::valarray<double>& th_path) const {
    const size_t max_sites = nsam_with_s_.size() - 1u;
    double loglik = (a_pathway_ * std::log(th_path)).sum();
    // D = 1.0 when s < 2
    for (size_t s=2u; s<=max_sites; ++s) {
        loglik -= nsam_with_s_[s] * std::log(calc_denom(w_pathway_, th_path, s));
    }
    return loglik += lnp_const_;
}

double PathtypeModel::calc_denom(
    const std::valarray<double>& w_pathway,
    const std::valarray<double>& th_pathway,
    const size_t num_mutations) const {

    if (num_mutations < 2u) return 1.0;
    auto iter = wtl::itertools::product(index_axes_[num_mutations]);
    double sum_prob = 0.0;
    std::bitset<128> bits(th_pathway.size());

    for (const auto& indices: iter()) {
        double p = 1.0;
        for (const auto j: indices) {
            p *= w_pathway[j];
            if (bits[j]) p *= th_pathway[j];
            bits.set(j);
        }
        sum_prob += p;
        bits.reset();
    }
    return sum_prob;
}

} // namespace likeligrid
