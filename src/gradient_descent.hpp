// -*- mode: c++; coding: utf-8 -*-
/*! @file gradient_descent.hpp
    @brief Interface of GradientDescent class
*/
#pragma once
#ifndef LIKELIGRID_GRADIENT_DESCENT_HPP_
#define LIKELIGRID_GRADIENT_DESCENT_HPP_

#include "genotype.hpp"

#include <string>
#include <vector>
#include <valarray>

namespace likeligrid {

class GradientDescent {
  public:
    GradientDescent() = default;
    GradientDescent(std::istream&,
        const size_t max_sites,
        const unsigned int concurrency=1);
    GradientDescent(std::istream&& ist,
        const size_t max_sites,
        const unsigned int concurrency=1)
        : GradientDescent(ist, max_sites, concurrency){}
    GradientDescent(
        const std::string& infile,
        const size_t max_sites=255,
        const unsigned int concurrency=1);

    void run();

    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    GenotypeModel model_;

    const unsigned int concurrency_;
};

} // namespace likeligrid

#endif // LIKELIGRID_GRADIENT_DESCENT_HPP_
