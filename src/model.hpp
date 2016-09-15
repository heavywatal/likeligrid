// -*- mode: c++; coding: utf-8 -*-
/*! @file model.hpp
    @brief Interface of Model class
*/
#pragma once
#ifndef LMPP_MODEL_HPP_
#define LMPP_MODEL_HPP_

#include <iostream>
#include <vector>
#include <map>

#include <Eigen/Core>

namespace lmpp {

/*! @brief Represents single run
*/
class Model {
  public:
    Model(std::istream& infile, const size_t g, const size_t n=65535);

    std::multimap<double, std::vector<double>> run(const double threshold, const double epsilon);

    const Eigen::MatrixXd& genotypes() const {return genotypes_;}

    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    const Eigen::MatrixXd genotypes_;
    const size_t grid_density_;
    const size_t max_results_;
    std::vector<Eigen::VectorXd> columns_;
};

extern std::ostream&
operator<<(std::ostream&, const std::multimap<double, std::vector<double>>&);

} // namespace lmpp

#endif /* LMPP_MODEL_HPP_ */
