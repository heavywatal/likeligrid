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

#include <Eigen/Core>
#include <wtl/itertools.hpp>

namespace likeligrid {

class ExclusivityModel {
  public:
    typedef Eigen::Array<uint, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ArrayXXu;
    typedef Eigen::Array<uint, Eigen::Dynamic, 1> ArrayXu;
    static const std::vector<double> STEPS_;
    static const std::vector<size_t> BREAKS_;

    ExclusivityModel() = default;
    ExclusivityModel(
        const std::string& infile,
        const size_t max_sites=255);
    void run(const std::string& infile="");

    static void raise_sigint() {SIGINT_RAISED_ = true;}
    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    void init_axes(const std::string&);
    void run_impl(const std::string&, wtl::itertools::Generator<Eigen::ArrayXd>&&) const;
    double calc_loglik(const Eigen::ArrayXd& th_path) const;
    double calc_denom(
        const Eigen::ArrayXd& weights,
        const Eigen::ArrayXd& exclusi,
        const size_t num_mutations) const;
    void search_limits() const;
    std::unordered_map<std::string, Eigen::ArrayXd> find_intersections() const;
    bool read_results(const std::string&);
    size_t read_metadata(std::istream&);
    size_t read_body(std::istream&);

    std::vector<std::string> names_;
    Eigen::ArrayXd w_pathway_;
    Eigen::ArrayXd a_pathway_;
    std::vector<size_t> nsam_with_s_;
    double lnp_const_ = 0.0;
    Eigen::ArrayXd mle_params_;
    std::vector<Eigen::ArrayXd> axes_;
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
