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

namespace likeligrid {

class PathtypeModel {
  public:
    PathtypeModel() = default;
    PathtypeModel(std::istream&&, size_t max_sites);
    PathtypeModel(
        const std::string& infile,
        size_t max_sites=255u);

    double calc_loglik(const std::valarray<double>& th_path) const;
    const std::vector<std::string>& names() const {return names_;}

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    double calc_denom(
        const std::valarray<double>& w_pathway,
        const std::valarray<double>& th_pathway,
        size_t num_mutations) const;

    std::vector<std::string> names_;
    std::valarray<double> w_pathway_;
    std::valarray<double> a_pathway_;
    std::vector<size_t> nsam_with_s_;
    double lnp_const_ = 0.0;
    std::vector<std::vector<std::vector<size_t>>> index_axes_;
};

} // namespace likeligrid

#endif /* LIKELIGRID_PATHTYPE_HPP_ */
