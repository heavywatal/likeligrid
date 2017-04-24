// -*- mode: c++; coding: utf-8 -*-
/*! @file typedef.hpp
    @brief Type definitions
*/
#pragma once
#ifndef LIKELIGRID_TYPEDEF_HPP_
#define LIKELIGRID_TYPEDEF_HPP_

#include <cstdint>
#include <bitset>

namespace likeligrid {

typedef uint_fast32_t uint;

typedef std::bitset<128> bits_t;
// TODO: Try boost::multiprecision::uint256_t

} // namespace likeligrid

#endif // LIKELIGRID_TYPEDEF_HPP_
