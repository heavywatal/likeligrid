// -*- mode: c++; coding: utf-8 -*-
/*! @file exclusivity.hpp
    @brief Interface of Exclusivity class
*/
#pragma once
#ifndef LIKELIGRID_EXCLUSIVITY_HPP_
#define LIKELIGRID_EXCLUSIVITY_HPP_

#include "typedef.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>

#include <wtl/itertools.hpp>

namespace likeligrid {

class ExclusivityModel {
  public:
    static const std::vector<double> STEPS_;
    static const std::vector<size_t> BREAKS_;

    ExclusivityModel() = default;
    ExclusivityModel(std::istream&&, const size_t max_sites);
    ExclusivityModel(
        const std::string& infile,
        const size_t max_sites=255);
    void run(const std::string& infile="");

    static void raise_sigint() {SIGINT_RAISED_ = true;}
    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    void run_impl(std::ostream&, wtl::itertools::Generator<std::valarray<double>>&&) const;
    double calc_loglik(const std::valarray<double>& th_path) const;
    double calc_denom(
        const std::valarray<double>& w_pathway,
        const std::valarray<double>& th_pathway,
        const size_t num_mutations) const;
    void search_limits() const;
    std::unordered_map<std::string, std::valarray<double>> find_intersections() const;
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

#endif /* LIKELIGRID_EXCLUSIVITY_HPP_ */
