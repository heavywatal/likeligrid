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
        ("max-sites,s", po::value(&MAX_SITES)->default_value(MAX_SITES))
        ("limits,l", po::value(&SEARCH_LIMITS)->default_value(SEARCH_LIMITS)->implicit_value(true))
        ("prefix,p", po::value(&PREFIX)->default_value(PREFIX), "prefix or previous result");
    return description;
}

po::options_description Program::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("genotypes", po::value(&GENOTYPES_FILE)->default_value(GENOTYPES_FILE));
    return description;
}

void Program::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "Usage: likeligrid [options] genotypes\n" << std::endl;
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
    positional.add("genotypes", 1);
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
    if (GENOTYPES_FILE == "-") {
        GENOTYPES_FILE = "/dev/stdin";
    }
}

void Program::run() {HERE;
    wtl::Fin fin(GENOTYPES_FILE);
    ExclusivityModel model(fin, MAX_SITES);
    model.run(PREFIX);
    if (SEARCH_LIMITS) model.search_limits();
}

} // namespace likeligrid
