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
#include <array>
#include <map>
#include <iterator>

namespace likeligrid {

constexpr std::array<double, 6> STEPS = {{0.32, 0.16, 0.08, 0.04, 0.02, 0.01}};
constexpr std::array<size_t, 6> BREAKS = {{5, 5, 5, 5, 5, 5}};

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

inline size_t guess_stage(const double step) {
    const auto it = std::find_if(STEPS.begin(), STEPS.end(), wtl::approx(step));
    if (it == STEPS.end()) throw std::runtime_error("invalid step size");
    return static_cast<size_t>(it - STEPS.begin());
}

inline double radius(const size_t stage) {
    return (BREAKS.at(stage) - 1u) * STEPS.at(stage) * 0.5;
}

inline std::tuple<std::string, size_t, size_t, double>
read_metadata(std::istream& ist) {
    std::string buffer;
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    std::getline(ist, buffer);
    const std::string genotype_file = buffer;
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    std::getline(ist, buffer);
    const size_t max_sites = std::stoul(buffer);
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    std::getline(ist, buffer);
    const size_t max_count = std::stoul(buffer);
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    std::getline(ist, buffer);
    const double step = std::stod(buffer);
    return std::make_tuple(genotype_file, max_sites, max_count, step);
}

inline double d2_from_neutral(std::valarray<double> v) {
    v -= 1.0;
    v *= v;
    return v.sum();
}

inline double d2_from_neutral(const std::vector<double>& v) {
    double d = 0.0;
    for (auto x: v) {
        x -= 1.0;
        x *= x;
        d += x;
    }
    return d;
}

inline std::tuple<size_t, std::vector<std::string>, std::valarray<double>>
read_body(std::istream& ist) {
    std::string buffer;
    ist >> buffer; // loglik
    std::getline(ist, buffer); // header
    buffer.erase(0u, 1u); // \t
    const std::vector<std::string> colnames = wtl::split(buffer, "\t");
    size_t nrow = 0u;
    double max_ll = std::numeric_limits<double>::lowest();
    std::vector<double> mle;
    while (std::getline(ist, buffer)) {
        ++nrow;
        std::istringstream iss(buffer);
        std::istream_iterator<double> it(iss);
        if (*it > max_ll) {
            max_ll = *it;
            mle.assign(++it, std::istream_iterator<double>());
        } else if (*it == max_ll) {
            std::vector<double> chalenger(++it, std::istream_iterator<double>());
            if (d2_from_neutral(chalenger) < d2_from_neutral(mle)) {
                mle.swap(chalenger);
            }
        }
    }
    std::valarray<double> mle_params(mle.data(), mle.size());
    return std::make_tuple(nrow, colnames, mle_params);
}

inline std::valarray<double>
read_loglik(std::istream& ist, const size_t nrow) {
    std::valarray<double> values(nrow);
    for (size_t i=0u; i<5u; ++i) {
        ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    for (size_t i=0u; i<nrow; ++i) {
        ist >> values[i];
        ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return values;
}

} // namespace likeligrid

#endif // LIKELIGRID_UTIL_HPP_
