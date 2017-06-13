// -*- mode: c++; coding: utf-8 -*-
/*! @file gridsearch.hpp
    @brief Interface of GridSearch class
*/
#pragma once
#ifndef LIKELIGRID_GRIDSEARCH_HPP_
#define LIKELIGRID_GRIDSEARCH_HPP_

#include "genotype.hpp"

#include <string>
#include <vector>
#include <valarray>

#include <wtl/itertools.hpp>

namespace likeligrid {

class GridSearch {
  public:
    GridSearch() = delete;
    GridSearch(std::istream&,
        const size_t max_sites,
        const std::pair<size_t, size_t>& epistasis_pair={0,0},
        const unsigned int concurrency=1);
    GridSearch(std::istream&& ist,
        const size_t max_sites,
        const std::pair<size_t, size_t>& epistasis_pair={0,0},
        const unsigned int concurrency=1)
        : GridSearch(ist, max_sites, epistasis_pair, concurrency){}
    GridSearch(
        const std::string& infile,
        const size_t max_sites,
        const std::pair<size_t, size_t>& epistasis_pair={0,0},
        const unsigned int concurrency=1);

    void run(const bool writing=true);

    void read_results(const std::string&);

    const std::valarray<double>& mle_params() const {return mle_params_;}

    static void test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    void run_fout();
    void run_cout();
    void run_impl(std::ostream&, wtl::itertools::Generator<std::valarray<double>>&&);
    void search_limits();
    std::string init_meta();
    void read_results(std::istream&);

    GenotypeModel model_;
    std::valarray<double> mle_params_;
    size_t skip_ = 0;
    size_t stage_ = 0;
    const unsigned int concurrency_;
};

} // namespace likeligrid

#endif // LIKELIGRID_GRIDSEARCH_HPP_
