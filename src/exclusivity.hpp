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

namespace likeligrid {

class Exclusivity {
  public:
    Exclusivity(std::istream& infile, const size_t g, const size_t n=65535);
    Exclusivity(const std::vector<std::string>& names,
          const Eigen::MatrixXd& genotypes,
          const size_t grid_density,
          const size_t max_results=65535):
          names_(names),
          genotypes_(genotypes),
          grid_density_(grid_density),
          max_results_(max_results) {}

    void run(const std::string& outfile="/dev/stdout");

    const std::vector<std::string>& names() const {return names_;}
    const Eigen::MatrixXd& genotypes() const {return genotypes_;}
    std::ostream& write_genotypes(std::ostream&, const bool header=true) const;
    std::ostream& write_results(std::ostream&, const bool header=true) const;

    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    const std::vector<std::string> names_;
    const Eigen::MatrixXd genotypes_;
    const size_t grid_density_;
    const size_t max_results_;
    std::multimap<double, std::vector<double>> results_;
};

} // namespace likeligrid

#endif /* LIKELIGRID_EXCLUSIVITY_HPP_ */
