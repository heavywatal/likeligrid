// -*- mode: c++; coding: utf-8 -*-
/*! @file exact.hpp
    @brief Interface of ExactModel class
*/
#pragma once
#ifndef LIKELIGRID_EXACT_HPP_
#define LIKELIGRID_EXACT_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cstdint>

#include <boost/dynamic_bitset.hpp>
#include <Eigen/Core>
#include <wtl/itertools.hpp>

namespace likeligrid {

typedef uint_fast32_t uint;

class ExactModel {
  public:
    typedef Eigen::Array<uint, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ArrayXXu;
    typedef Eigen::Array<uint, Eigen::Dynamic, 1> ArrayXu;

    ExactModel() = default;
    ExactModel(
        const std::string& infile,
        const size_t max_sites=-1);
    void run();

    static void raise_sigint() {SIGINT_RAISED_ = true;}
    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    double calc_loglik(const Eigen::ArrayXd& params) const;
    double calc_denom(
        const Eigen::ArrayXd& params,
        const size_t num_mutations) const;

    std::vector<std::string> names_;
    std::vector<boost::dynamic_bitset<>> annot_;
    std::vector<boost::dynamic_bitset<>> genot_;
    Eigen::ArrayXd w_gene_;
    Eigen::ArrayXd a_pathway_;
    std::vector<size_t> nsam_with_s_;
    double lnp_const_ = 0.0;
    Eigen::ArrayXd mle_params_;

    static bool SIGINT_RAISED_;
};

} // namespace likeligrid

#endif // LIKELIGRID_EXACT_HPP
