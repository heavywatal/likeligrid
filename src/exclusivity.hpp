// -*- mode: c++; coding: utf-8 -*-
/*! @file exclusivity.hpp
    @brief Interface of Exclusivity class
*/
#pragma once
#ifndef LIKELIGRID_EXCLUSIVITY_HPP_
#define LIKELIGRID_EXCLUSIVITY_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include <Eigen/Core>
#include <cxxwtils/itertools.hpp>

namespace likeligrid {

class ExclusivityModel {
  public:
    typedef Eigen::Array<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ArrayXXu;
    typedef Eigen::Array<size_t, Eigen::Dynamic, 1> ArrayXu;
    static const std::vector<double> STEPS_;
    static const std::vector<size_t> BREAKS_;

    ExclusivityModel(std::istream& genotypes, const size_t max_sites=-1);

    void run(const std::string& infile="");
    void search_limits() const;

    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    void init_axes(const std::string&);
    std::string name_outfile(const std::string&) const;
    void run_impl(const std::string&, wtl::itertools::Generator<Eigen::ArrayXd>&&) const;
    double calc_loglik(const Eigen::ArrayXd& params) const;
    double calc_denom(
        const Eigen::ArrayXd& weights,
        const Eigen::ArrayXd& exclusi,
        const size_t num_mutations) const;
    std::unordered_map<std::string, Eigen::ArrayXd> find_intersections() const;
    std::ostream& write_genotypes(std::ostream&, const bool header=true) const;
    bool read_results(const std::string&);
    size_t read_metadata(std::istream&);
    size_t read_body(std::istream&);

    const std::vector<std::string> names_;
    ArrayXXu genotypes_;
    Eigen::ArrayXd w_pathway_;
    Eigen::ArrayXd a_pathway_;
    std::vector<size_t> nsam_with_s_;
    double lnp_const_ = 0.0;
    Eigen::ArrayXd mle_params_;
    std::vector<Eigen::ArrayXd> axes_;
    size_t start_ = 0;
    size_t stage_ = 0;
    std::vector<std::vector<std::vector<size_t>>> index_axes_;
};

} // namespace likeligrid

#endif /* LIKELIGRID_EXCLUSIVITY_HPP_ */
