// -*- mode: c++; coding: utf-8 -*-
/*! @file model.hpp
    @brief Interface of Model class
*/
#pragma once
#ifndef LMPP_MODEL_HPP_
#define LMPP_MODEL_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <Eigen/Core>

namespace lmpp {

class Model {
  public:
    Model(std::istream& infile, const size_t g, const size_t n=65535);
    Model(const std::vector<std::string>& names,
          const Eigen::MatrixXd& genotypes,
          const size_t grid_density,
          const size_t max_results=65535):
          names_(names),
          genotypes_(genotypes),
          grid_density_(grid_density),
          max_results_(max_results) {}

    void run(const double threshold, const double epsilon,
             const std::string& outfile="/dev/stdout");

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
    std::vector<Eigen::VectorXd> columns_;
    std::multimap<double, std::vector<double>> results_;
};

} // namespace lmpp

#endif /* LMPP_MODEL_HPP_ */
