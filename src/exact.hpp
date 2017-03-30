// -*- mode: c++; coding: utf-8 -*-
/*! @file exact.hpp
    @brief Interface of ExactModel class
*/
#pragma once
#ifndef LIKELIGRID_EXACT_HPP_
#define LIKELIGRID_EXACT_HPP_

#include <cstdint>
#include <string>
#include <vector>
#include <valarray>

#include <boost/dynamic_bitset.hpp>

namespace likeligrid {

typedef uint_fast32_t uint;
typedef boost::dynamic_bitset<> bits_t;

class ExactModel {
  public:
    ExactModel() = default;
    ExactModel(
        const std::string& infile,
        const size_t max_sites=-1);
    void run();

    static void raise_sigint() {SIGINT_RAISED_ = true;}
    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    double calc_loglik(const std::valarray<double>& params) const;
    double calc_denom(
        const std::valarray<double>& params,
        const size_t num_mutations) const;

    std::vector<std::string> names_;
    std::vector<bits_t> annot_;
    std::vector<bits_t> genot_;
    std::valarray<double> w_gene_;
    std::valarray<double> a_pathway_;
    std::vector<size_t> nsam_with_s_;
    double lnp_const_ = 0.0;
    std::valarray<double> mle_params_;

    static bool SIGINT_RAISED_;
};

} // namespace likeligrid

#endif // LIKELIGRID_EXACT_HPP
