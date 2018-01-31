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

namespace wtl {namespace itertools {
  template <class T> class Generator;
}}

namespace likeligrid {

class GridSearch {
  public:
    GridSearch() = delete;
    GridSearch(std::istream& ist,
        const size_t max_sites,
        const std::pair<size_t, size_t>& epistasis_pair={0u,0u},
        const bool pleiotropy=false,
        const unsigned int concurrency=1u)
        : model_(ist, max_sites),
          concurrency_(concurrency) {
        init(epistasis_pair, pleiotropy);
    }
    GridSearch(std::istream&& ist,
        const size_t max_sites,
        const std::pair<size_t, size_t>& epistasis_pair={0u,0u},
        const bool pleiotropy=false,
        const unsigned int concurrency=1u)
        : GridSearch(ist, max_sites, epistasis_pair, pleiotropy, concurrency){}
    GridSearch(
        const std::string& infile,
        const size_t max_sites,
        const std::pair<size_t, size_t>& epistasis_pair={0u,0u},
        const bool pleiotropy=false,
        const unsigned int concurrency=1u)
        : model_(infile, max_sites),
          concurrency_(concurrency) {
        init(epistasis_pair, pleiotropy);
    }

    void run(const bool writing=true);
    void run_cout();

    void read_results(const std::string&);

    const std::valarray<double>& mle_params() const {return mle_params_;}

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    void init(const std::pair<size_t, size_t>&, const bool pleiotropy);
    void run_fout();
    void run_impl(std::ostream&, wtl::itertools::Generator<std::valarray<double>>&&);
    void search_limits();
    std::string init_meta();
    void read_results(std::istream&);
    void write_header(std::ostream&, const size_t max_count) const;

    GenotypeModel model_;
    std::valarray<double> mle_params_;
    size_t skip_ = 0u;
    size_t stage_ = 0u;
    const unsigned int concurrency_;
};

} // namespace likeligrid

#endif // LIKELIGRID_GRIDSEARCH_HPP_
