// -*- mode: c++; coding: utf-8 -*-
/*! @file program.cpp
    @brief Inplementation of Program class
*/
#include "program.hpp"

#include <cstdlib>

#include <cxxwtils/debug.hpp>
#include <cxxwtils/iostr.hpp>
#include <cxxwtils/getopt.hpp>
#include <cxxwtils/eigen.hpp>

#include "exclusivity.hpp"

namespace likeligrid {

namespace po = boost::program_options;

inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "print this help")
        ("verbose,v", po::value<bool>()->default_value(false)->implicit_value(true), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1));
    return description;
}

po::options_description Program::options_desc() {HERE;
    po::options_description description("Program");
    description.add_options()
        ("grid,g", po::value(&GRID_DENSITY)->default_value(GRID_DENSITY))
        ("results,n", po::value(&MAX_RESULTS)->default_value(MAX_RESULTS))
        ("outfile,o", po::value<std::string>(&OUTFILE)->default_value(OUTFILE));
    return description;
}

po::options_description Program::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("infile", po::value(&INFILE)->default_value(INFILE));
    return description;
}

void Program::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "Usage: likeligrid [options] infile\n" << std::endl;
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

//! Unit test for each class
inline void test(const int flg) {HERE;
    switch (flg) {
      case 0:
        break;
      case 1:
        ExclusivityModel::unit_test();
        throw wtl::ExitSuccess();
      default:
        throw std::runtime_error("Unknown argument for --test");
    }
}

Program::Program(const std::vector<std::string>& arguments) {HERE;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);
    COMMAND_ARGS = wtl::str_join(arguments, " ");

    auto description = general_desc();
    description.add(options_desc());
    description.add(positional_desc());
    po::positional_options_description positional;
    positional.add("infile", 1);
    po::variables_map vm;
    po::store(po::command_line_parser(arguments).
              options(description).
              positional(positional).run(), vm);
    if (vm["help"].as<bool>()) {help_and_exit();}
    po::notify(vm);

    CONFIG_STRING = wtl::flags_into_string(vm);
    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << CONFIG_STRING << std::endl;
    }
    test(vm["test"].as<int>());
    if (INFILE == "-") {
        INFILE = "/dev/stdin";
    }
    if (OUTFILE == "-") {
        OUTFILE = "/dev/stdout";
    }
}

void Program::run() {HERE;
    wtl::Fin fin(INFILE);
    ExclusivityModel model(fin, GRID_DENSITY, MAX_RESULTS);
    model.run(OUTFILE);
}

} // namespace likeligrid
