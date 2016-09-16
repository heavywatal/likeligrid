// -*- mode: c++; coding: utf-8 -*-
/*! @file program.cpp
    @brief Inplementation of Program class
*/
#include "program.hpp"

#include <cstdlib>

#include <boost/filesystem.hpp>

#include <cxxwtils/debug.hpp>
#include <cxxwtils/iostr.hpp>
#include <cxxwtils/getopt.hpp>
#include <cxxwtils/os.hpp>
#include <cxxwtils/gz.hpp>

#include "model.hpp"

namespace lmpp {

namespace fs = boost::filesystem;
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
        ("max,n", po::value(&MAX_RESULTS)->default_value(MAX_RESULTS))
        ("epsilon,e", po::value(&EPSILON)->default_value(EPSILON))
        ("threshold,c", po::value(&THRESHOLD)->default_value(THRESHOLD))
        ("outfile,o", po::value<std::string>(&OUTFILE)->default_value(OUTFILE));
    return description;
}

po::options_description Program::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("genotype", po::value(&GENOTYPE_FILE)->default_value(GENOTYPE_FILE));
    return description;
}

void Program::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "Usage: lmpp [options] genotype\n" << std::endl;
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

//! Unit test for each class
inline void test(const int flg) {HERE;
    switch (flg) {
      case 0:
        break;
      case 1:
        Model::unit_test();
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
    positional.add("genotype", 1);
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
}

void Program::run() {HERE;
    wtl::Fin fin(GENOTYPE_FILE);
    Model model(fin, GRID_DENSITY, MAX_RESULTS);
    const auto result = model.run(THRESHOLD, EPSILON);
    if (OUTFILE == "-") {
        print(std::cout, result);
    } else {
        print(wtl::Fout(OUTFILE), result);
    }
}

} // namespace lmpp
