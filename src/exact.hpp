// -*- mode: c++; coding: utf-8 -*-
/*! @file exact.hpp
    @brief Interface of ExactModel class
*/
#pragma once
#ifndef LIKELIGRID_EXACT_HPP_
#define LIKELIGRID_EXACT_HPP_

#include "typedef.hpp"

#include <string>
#include <vector>
#include <valarray>

#include <wtl/itertools.hpp>

namespace likeligrid {

class ExactModel {
  public:
    static const std::vector<double> STEPS_;
    static const std::vector<size_t> BREAKS_;

    ExactModel() = default;
    ExactModel(std::istream&&, const size_t max_sites);
    ExactModel(
        const std::string& infile,
        const size_t max_sites=255);
    void run(const std::string& infile="");

    static void raise_sigint() {SIGINT_RAISED_ = true;}
    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    void run_impl(std::ostream&, wtl::itertools::Generator<std::valarray<double>>&&) const;
    double calc_loglik(const std::valarray<double>& th_path) const;
    std::string init_meta(const std::string& infile);
    bool read_results(const std::string&);

    std::vector<std::string> names_;
    std::vector<bits_t> annot_;
    std::vector<bits_t> genot_;
    std::valarray<double> w_gene_;
    std::valarray<double> a_pathway_;
    std::vector<size_t> nsam_with_s_;
    double lnp_const_ = 0.0;

    std::valarray<double> mle_params_;
    size_t skip_ = 0;
    size_t stage_ = 0;

    static bool SIGINT_RAISED_;
};

} // namespace likeligrid

#endif // LIKELIGRID_EXACT_HPP
