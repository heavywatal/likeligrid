// -*- mode: c++; coding: utf-8 -*-
/*! @file program.hpp
    @brief Interface of Program class
*/
#pragma once
#ifndef LIKELIGRID_PROGRAM_HPP_
#define LIKELIGRID_PROGRAM_HPP_

#include <iostream>
#include <sstream>
#include <vector>

#include <cxxwtils/exception.hpp>

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace likeligrid {

/*! @brief Represents single run
*/
class Program {
  public:
    //! Parse command arguments
    Program(const std::vector<std::string>& args);

    //! Top level function that should be called once from main()
    void run();

    std::string conf() const {return CONFIG_STRING;}

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    boost::program_options::options_description options_desc();
    boost::program_options::options_description positional_desc();
    void help_and_exit();

    size_t GRID_DENSITY = 10;
    size_t MAX_SITES = 65535;
    size_t MAX_RESULTS = 16;
    std::string INFILE = "-";
    std::string OUTFILE = "-";

    std::string COMMAND_ARGS;
    std::string CONFIG_STRING;
};

} // namespace likeligrid

#endif /* LIKELIGRID_PROGRAM_HPP_ */
