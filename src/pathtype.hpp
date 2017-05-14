// -*- mode: c++; coding: utf-8 -*-
/*! @file pathtype.hpp
    @brief Interface of PathtypeModel class
*/
#pragma once
#ifndef LIKELIGRID_PATHTYPE_HPP_
#define LIKELIGRID_PATHTYPE_HPP_

#include "typedef.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>

#include <wtl/itertools.hpp>

namespace likeligrid {

class PathtypeModel {
  public:
    static const std::vector<double> STEPS_;
    static const std::vector<size_t> BREAKS_;

    PathtypeModel() = default;
    PathtypeModel(std::istream&&, const size_t max_sites);
    PathtypeModel(
        const std::string& infile,
        const size_t max_sites=255);
    void run(const std::string& infile="");

    double calc_loglik(const std::valarray<double>& th_path) const;
    const std::valarray<double>& mle_params() const {return mle_params_;}
    const std::vector<std::string>& names() const {return names_;}

    static void raise_sigint() {SIGINT_RAISED_ = true;}
    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    void run_impl(std::ostream&, wtl::itertools::Generator<std::valarray<double>>&&) const;
    double calc_denom(
        const std::valarray<double>& w_pathway,
        const std::valarray<double>& th_pathway,
        const size_t num_mutations) const;
    void search_limits() const;
    std::string init_meta(const std::string& infile);
    bool read_results(const std::string&);

    std::vector<std::string> names_;
    std::valarray<double> w_pathway_;
    std::valarray<double> a_pathway_;
    std::vector<size_t> nsam_with_s_;
    double lnp_const_ = 0.0;
    std::valarray<double> mle_params_;
    size_t skip_ = 0;
    size_t stage_ = 0;
    std::vector<std::vector<std::vector<size_t>>> index_axes_;

    static bool SIGINT_RAISED_;
};


class lnpnan_error: public std::runtime_error {
  public:
    lnpnan_error():
    std::runtime_error("lnp is nan; maybe some pathways have no mutation") {}
};


} // namespace likeligrid

#endif /* LIKELIGRID_PATHTYPE_HPP_ */
