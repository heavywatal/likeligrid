/*! @file pathtype.hpp
    @brief Interface of PathtypeModel class
*/
#pragma once
#ifndef LIKELIGRID_PATHTYPE_HPP_
#define LIKELIGRID_PATHTYPE_HPP_

#include <iosfwd>
#include <string>
#include <vector>
#include <valarray>
#include <stdexcept>

namespace likeligrid {

class PathtypeModel {
  public:
    PathtypeModel() = default;
    PathtypeModel(std::istream&&, const size_t max_sites);
    PathtypeModel(
        const std::string& infile,
        const size_t max_sites=255u);

    double calc_loglik(const std::valarray<double>& th_path) const;
    const std::vector<std::string>& names() const {return names_;}

    static void test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    double calc_denom(
        const std::valarray<double>& w_pathway,
        const std::valarray<double>& th_pathway,
        const size_t num_mutations) const;

    std::vector<std::string> names_;
    std::valarray<double> w_pathway_;
    std::valarray<double> a_pathway_;
    std::vector<size_t> nsam_with_s_;
    double lnp_const_ = 0.0;
    std::vector<std::vector<std::vector<size_t>>> index_axes_;
};


class lnpnan_error: public std::runtime_error {
  public:
    lnpnan_error():
    std::runtime_error("lnp is nan; maybe some pathways have no mutation") {}
};


} // namespace likeligrid

#endif /* LIKELIGRID_PATHTYPE_HPP_ */
