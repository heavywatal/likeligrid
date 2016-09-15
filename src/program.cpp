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
#include <cxxwtils/prandom.hpp>
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
        ("grid,g", po::value<size_t>(&GRID_DENSITY)->default_value(GRID_DENSITY))
        ("max,n", po::value<size_t>(&MAX_RESULTS)->default_value(MAX_RESULTS))
        ("out_dir,o", po::value<std::string>(&OUT_DIR)->default_value(OUT_DIR))
        ("seed", po::value<unsigned int>(&SEED)->default_value(SEED));
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
    OUT_DIR = wtl::strftime("lmpp_%Y%m%d_%H%M_") + std::to_string(::getpid());

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
    wtl::sfmt().seed(SEED);

    // CONFIG_STRING = wtl::flags_into_string(vm);
    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << CONFIG_STRING << std::endl;
    }
    test(vm["test"].as<int>());
    OUT_DIR = fs::system_complete(OUT_DIR).string();
}

void Program::run() {HERE;
    wtl::Fin fin(GENOTYPE_FILE);
    Model model(fin, GRID_DENSITY, MAX_RESULTS);
    std::vector<double> v_thr{0.5};
    std::vector<double> v_eps{0.0, 0.1};
    for (const auto threshold: v_thr) {
        for (const auto epsilon: v_eps) {
            std::cout << model.run(threshold, epsilon) << std::endl;
        }
    }
}

void Program::write() const {HERE;
    if (false) {
        derr("mkdir && cd to " << OUT_DIR << std::endl);
        fs::create_directory(OUT_DIR);
        wtl::cd(OUT_DIR);
        wtl::Fout{"program_options.conf"} << CONFIG_STRING;
        std::cerr << wtl::iso8601datetime() << std::endl;
    }
}

} // namespace lmpp
