// -*- mode: c++; coding: utf-8 -*-
/*! @file util.hpp
    @brief Common functions
*/
#pragma once
#ifndef LIKELIGRID_UTIL_HPP_
#define LIKELIGRID_UTIL_HPP_

#include <boost/math/distributions/chi_squared.hpp>

#include <wtl/iostr.hpp>
#include <wtl/numeric.hpp>

#include <string>
#include <vector>
#include <valarray>
#include <map>
#include <iterator>

namespace likeligrid {

inline std::vector<std::valarray<double>>
make_vicinity(const std::valarray<double>& center, const size_t breaks, const double radius, const double max=2.001) {
    std::vector<std::valarray<double>> axes;
    axes.reserve(center.size());
    for (const double x: center) {
        auto axis = wtl::lin_spaced(breaks, x + radius, x - radius);
        // grid precision = 0.01
        axis = (axis * 100.0).apply(std::round) / 100.0;
        const std::valarray<bool> positive = axis > 0.0;
        axes.emplace_back(axis[positive & (axis < max)]);
    }
    return axes;
}

inline size_t guess_stage(const std::vector<double>& STEPS, const double step) {
    const auto it = std::find_if(STEPS.begin(), STEPS.end(), wtl::approx(step));
    if (it == STEPS.end()) throw std::runtime_error("invalid step size");
    return it - STEPS.begin();
}

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

inline std::valarray<double>
read_loglik(std::istream& ist, const size_t nrow) {
    std::valarray<double> values(nrow);
    std::string buffer;
    std::getline(ist, buffer);
    std::getline(ist, buffer);
    std::getline(ist, buffer);
    std::getline(ist, buffer); // header
    for (size_t l=0; l<nrow; ++l) {
        ist >> values[l];
        std::getline(ist, buffer);
    }
    return values;
}

} // namespace likeligrid

#endif // LIKELIGRID_UTIL_HPP_
