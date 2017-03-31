// -*- mode: c++; coding: utf-8 -*-
/*! @file program.cpp
    @brief Inplementation of Program class
*/
#include "program.hpp"

#include <csignal>
#include <cstdlib>
#include <memory>
#include <regex>

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/os.hpp>
#include <wtl/getopt.hpp>
#include <wtl/eigen.hpp>

#include "exclusivity.hpp"
#include "exact.hpp"

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
        ("max-sites,s", po::value(&MAX_SITES)->default_value(MAX_SITES));
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
        ExactModel::unit_test();
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
    std::signal(SIGINT, [](int signum){
        if (signum == SIGINT) {
            ExclusivityModel::raise_sigint();
            ExactModel::raise_sigint();
        }
    });

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

    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << wtl::flags_into_string(vm) << std::endl;
    }
    test(vm["test"].as<int>());
}

void Program::run() {HERE;
    std::string prefix = "stdin";
    if (GENOTYPES_FILE == "-") {
        GENOTYPES_FILE = "/dev/stdin";
    } else {
        std::smatch mobj;
        std::regex_search(GENOTYPES_FILE, mobj, std::regex("([^/]+?)\\.[^/]+$"));
        prefix = mobj.str(1);
    }
    std::ostringstream oss;
    oss << prefix << "-s" << MAX_SITES;
    const std::string outdir = oss.str();
    try {
        ExclusivityModel model(GENOTYPES_FILE, MAX_SITES);
        wtl::mkdir(outdir);  // after constructor success
        wtl::Pushd cd(outdir);
        model.run();
    } catch (const lnpnan_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

} // namespace likeligrid
