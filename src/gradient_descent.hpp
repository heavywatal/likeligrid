/*! @file gradient_descent.hpp
    @brief Interface of GradientDescent class
*/
#pragma once
#ifndef LIKELIGRID_GRADIENT_DESCENT_HPP_
#define LIKELIGRID_GRADIENT_DESCENT_HPP_

#include <iosfwd>
#include <string>
#include <vector>
#include <valarray>
#include <map>
#include <memory>

namespace likeligrid {

class GenotypeModel;

class lexicographical_less {
  public:
    bool operator() (const std::valarray<double>& x, const std::valarray<double>&y) const {
        return std::lexicographical_compare(std::begin(x), std::end(x), std::begin(y), std::end(y));
    }
};

//! existing keys are checked in empty_neighbors_of()
using MapGrid = std::map<std::valarray<double>, double, lexicographical_less>;

class GradientDescent {
  public:
    GradientDescent() = delete;
    GradientDescent(std::istream& ist,
        const size_t max_sites,
        const std::pair<size_t, size_t>& epistasis_pair={0u,0u},
        const bool pleiotropy=false,
        const unsigned int concurrency=1u);
    GradientDescent(
        const std::string& infile,
        const size_t max_sites,
        const std::pair<size_t, size_t>& epistasis_pair={0u,0u},
        const bool pleiotropy=false,
        const unsigned int concurrency=1u);
    ~GradientDescent();

    void run(std::ostream&);

    std::string outfile() const {return outfile_;};

    static void test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    MapGrid::iterator find_better(const MapGrid::iterator&);
    std::vector<std::valarray<double>> empty_neighbors_of(const std::valarray<double>&);

    void write(std::ostream&);
    std::tuple<std::string, size_t, std::string> read_results(const std::string&);

    MapGrid::iterator max_iterator();
    MapGrid::const_iterator const_max_iterator() const;

    std::unique_ptr<GenotypeModel> model_;
    std::valarray<double> starting_point_;
    MapGrid history_;
    std::string outfile_;

    const unsigned int concurrency_;
};

} // namespace likeligrid

#endif // LIKELIGRID_GRADIENT_DESCENT_HPP_
