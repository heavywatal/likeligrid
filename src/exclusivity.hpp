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
#include <map>

#include <Eigen/Core>
#include <cxxwtils/itertools.hpp>

namespace likeligrid {

class ExclusivityModel {
  public:
    typedef Eigen::Array<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ArrayXXu;
    typedef Eigen::Array<size_t, Eigen::Dynamic, 1> ArrayXu;
    static const std::vector<double> STEPS_;
    static const std::vector<size_t> BREAKS_;

    ExclusivityModel(std::istream& genotypes,
        const size_t max_sites=65535,
        const size_t max_results=65535);

    void run(const std::string& infile="");

    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    void init_axes(const std::string&);
    std::string name_outfile(const std::string&) const;
    void run_impl(const std::string&, wtl::itertools::Generator<Eigen::ArrayXd>&&);
    double calc_denom(
        const Eigen::ArrayXd& weights,
        const Eigen::ArrayXd& exclusi,
        const size_t num_mutations);
    std::ostream& write_genotypes(std::ostream&, const bool header=true) const;
    std::ostream& write_results(std::ostream&, const size_t max_rows=-1) const;
    bool read_results(const std::string&, const size_t max_rows=-1);
    void read_metadata(std::istream&);
    void read_body(std::istream&, const size_t max_rows);

    const std::vector<std::string> names_;
    ArrayXXu genotypes_;
    size_t max_results_;
    Eigen::ArrayXd w_pathway_;
    Eigen::ArrayXd a_pathway_;
    std::vector<size_t> nsam_with_s_;
    double lnp_const_ = 0.0;
    std::multimap<double, Eigen::ArrayXd, std::greater<double>> results_;
    Eigen::ArrayXd best_;
    std::vector<Eigen::ArrayXd> axes_;
    size_t start_ = 0;
    size_t stage_ = 0;
    std::vector<wtl::itertools::Product<std::vector<size_t>>> index_iters_;
};

} // namespace likeligrid

#endif /* LIKELIGRID_EXCLUSIVITY_HPP_ */
