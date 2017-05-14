// -*- mode: c++; coding: utf-8 -*-
/*! @file exact.hpp
    @brief Interface of ExactModel class
*/
#pragma once
#ifndef LIKELIGRID_EXACT_HPP_
#define LIKELIGRID_EXACT_HPP_

#include "typedef.hpp"

#include <string>
#include <vector>
#include <valarray>

#include <wtl/itertools.hpp>

namespace likeligrid {

class ExactModel {
  public:
    static const std::vector<double> STEPS_;
    static const std::vector<size_t> BREAKS_;

    ExactModel() = default;
    ExactModel(std::istream&,
        const size_t max_sites,
        const unsigned int concurrency=1);
    ExactModel(std::istream&& ist,
        const size_t max_sites,
        const unsigned int concurrency=1)
        : ExactModel(ist, max_sites, concurrency){}
    ExactModel(
        const std::string& infile,
        const size_t max_sites=255,
        const unsigned int concurrency=1);
    void run(const bool writing=true);

    double calc_loglik(const std::valarray<double>& th_path) const;
    const std::valarray<double>& mle_params() const {return mle_params_;}
    const std::vector<std::string>& names() const {return names_;}
    size_t num_genes() const {return w_gene_.size();}

    static void raise_sigint() {SIGINT_RAISED_ = true;}
    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    void run_fout();
    void run_cout();
    void run_impl(std::ostream&, wtl::itertools::Generator<std::valarray<double>>&&) const;
    void search_limits() const;
    std::string init_meta();
    void read_results(std::istream&);

    size_t breaks() const {return BREAKS_.at(stage_);}
    double radius() const {return (breaks() - 1) * STEPS_.at(stage_) * 0.5;}

    std::vector<std::string> names_;
    std::vector<bits_t> annot_;
    std::vector<bits_t> genot_;
    std::valarray<double> w_gene_;
    std::vector<size_t> nsam_with_s_;

    std::valarray<double> mle_params_;
    size_t skip_ = 0;
    size_t stage_ = 0;
    const unsigned int concurrency_;

    static bool SIGINT_RAISED_;
};

} // namespace likeligrid

#endif // LIKELIGRID_EXACT_HPP
