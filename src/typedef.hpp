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

using uint =  uint_fast32_t;

using bits_t = std::bitset<128>;
// TODO: Try boost::multiprecision::uint256_t

} // namespace likeligrid

#endif // LIKELIGRID_TYPEDEF_HPP_
