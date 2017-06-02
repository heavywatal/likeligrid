// -*- mode: c++; coding: utf-8 -*-
/*! @file gradient_descent.hpp
    @brief Interface of GradientDescent class
*/
#pragma once
#ifndef LIKELIGRID_GRADIENT_DESCENT_HPP_
#define LIKELIGRID_GRADIENT_DESCENT_HPP_

#include "genotype.hpp"

#include <iosfwd>
#include <string>
#include <vector>
#include <valarray>
#include <map>

namespace likeligrid {

class lexicographical_less {
  public:
    bool operator() (const std::valarray<double>& x, const std::valarray<double>&y) const {
        return std::lexicographical_compare(std::begin(x), std::end(x), std::begin(y), std::end(y));
    }
};

typedef std::map<std::valarray<double>, double, lexicographical_less> MapGrid;

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
        const size_t max_sites,
        const unsigned int concurrency=1);

    void run(const std::pair<size_t, size_t>& epistasis_pair={0,0});
    void run(std::ostream&, const std::pair<size_t, size_t>& epistasis_pair={0,0});

    MapGrid::iterator max_iterator();
    MapGrid::const_iterator const_max_iterator() const;

    void write(std::ostream&, const std::pair<size_t, size_t>& epistasis_pair={0,0});
    std::string read_results(const std::string&, const size_t max_sites);

    static void test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    MapGrid::iterator find_better(const MapGrid::iterator&);
    std::vector<std::valarray<double>> empty_neighbors_of(const std::valarray<double>&);

    MapGrid history_;
    GenotypeModel model_;

    const unsigned int concurrency_;
};

} // namespace likeligrid

#endif // LIKELIGRID_GRADIENT_DESCENT_HPP_
