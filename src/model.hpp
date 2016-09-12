// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.h
    @brief Interface of Simulation class
*/
#pragma once
#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <random>

#include <cxxwtils/exception.hpp>

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace lmpp {

/*! @brief Represents single run
*/
class Model {
  public:
    //! Parse command arguments
    Model(const std::vector<std::string>& args);

    //! Top level function that should be called once from main()
    void run();

    void write() const;

    std::string conf() const {return CONFIG_STRING;}

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    boost::program_options::options_description options_desc();
    boost::program_options::options_description positional_desc();
    void help_and_exit();

    double normal_lccdf(double mean, const double sd);
    void genotype2score();

    double epsilon_ = 0.1;
    double threshold_ = 0.5;

    bool WRITE_TO_FILES = false;

    //! Seed for random number generator
    unsigned int SEED = std::random_device{}();

    //! Target directory to which the contents in WORK_DIR are moved
    std::string OUT_DIR;

    std::string COMMAND_ARGS;
    std::string CONFIG_STRING;
};

} // namespace tumopp

#endif /* MODEL_HPP_ */
