// -*- mode: c++; coding: utf-8 -*-
/*! @file model.hpp
    @brief Interface of Model class
*/
#pragma once
#ifndef LMPP_MODEL_HPP_
#define LMPP_MODEL_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <random>

namespace lmpp {

/*! @brief Represents single run
*/
class Model {
  public:
    Model(const double t, const double e):
        threshold_(t), epsilon_(e) {}

    //! Top level function that should be called once from main()
    void run();

    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    double likelihood(double score);
    void genotype2score();

    double threshold_ = 0.5;
    double epsilon_ = 0.1;
};

} // namespace lmpp

#endif /* LMPP_MODEL_HPP_ */
