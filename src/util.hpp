// -*- mode: c++; coding: utf-8 -*-
/*! @file util.hpp
    @brief Common functions
*/
#pragma once
#ifndef LIKELIGRID_UTIL_HPP_
#define LIKELIGRID_UTIL_HPP_

#include <wtl/iostr.hpp>

#include <string>
#include <vector>
#include <valarray>
#include <iterator>

namespace likeligrid {

inline std::tuple<size_t, size_t, double>
read_metadata(std::istream& ist) {
    std::string buffer;
    std::getline(ist, buffer, '=');
    std::getline(ist, buffer);
    const size_t max_count = std::stoul(buffer);
    std::getline(ist, buffer, '=');
    std::getline(ist, buffer);
    const size_t max_sites = std::stoul(buffer);
    std::getline(ist, buffer, '=');
    std::getline(ist, buffer);
    const double step = std::stod(buffer);
    return std::make_tuple(max_count, max_sites, step);
}

inline std::tuple<size_t, std::vector<std::string>, std::valarray<double>>
read_body(std::istream& ist) {
    std::string buffer;
    ist >> buffer; // loglik
    std::getline(ist, buffer); // header
    buffer.erase(0, 1); // \t
    const std::vector<std::string> colnames = wtl::split(buffer, "\t");
    size_t nrow = 0;
    double max_ll = std::numeric_limits<double>::lowest();
    std::vector<double> mle;
    while (std::getline(ist, buffer)) {
        ++nrow;
        std::istringstream iss(buffer);
        std::istream_iterator<double> it(iss);
        if (*it > max_ll) {
            max_ll = *it;
            mle.assign(++it, std::istream_iterator<double>());
        }
    }
    std::valarray<double> mle_params(mle.data(), mle.size());
    return std::make_tuple(nrow, colnames, mle_params);
}

} // namespace likeligrid

#endif // LIKELIGRID_UTIL_HPP_
