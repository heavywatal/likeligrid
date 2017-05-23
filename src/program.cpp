// -*- mode: c++; coding: utf-8 -*-
/*! @file program.cpp
    @brief Implementation of Program class
*/
#include "program.hpp"
#include "pathtype.hpp"
#include "genotype.hpp"
#include "gridsearch.hpp"
#include "gradient_descent.hpp"

#include <csignal>
#include <cstdlib>
#include <memory>
#include <regex>

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/os.hpp>
#include <wtl/getopt.hpp>
#include <wtl/zfstream.hpp>

namespace likeligrid {

bool SIGINT_RAISED = false;

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
        ("parallel,j", po::value(&CONCURRENCY)->default_value(CONCURRENCY))
        ("max-sites,s", po::value(&MAX_SITES)->default_value(MAX_SITES))
        ("previous,p", po::value(&PREVIOUS_RESULT)->default_value(PREVIOUS_RESULT));
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
void Program::test(const int flag) {HERE;
    switch (flag) {
      case 0:
        break;
      case 1:
        GridSearch::unit_test();
        GradientDescent::unit_test();
        throw wtl::ExitSuccess();
      case 2:
        GenotypeModel::unit_test();
        PathtypeModel::unit_test();
        throw wtl::ExitSuccess();
      case 3: {
        wtl::izfstream ist(GENOTYPES_FILE);
        GenotypeModel model(ist, MAX_SITES);
        model.benchmark(CONCURRENCY);
        throw wtl::ExitSuccess();
      }
      default:
        throw std::runtime_error("Unknown argument for --test");
    }
}

Program::Program(const std::vector<std::string>& arguments) {HERE;
    std::cout << wtl::join(arguments, " ") << std::endl;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);
    std::signal(SIGINT, [](int signum){
        if (signum == SIGINT) {
            SIGINT_RAISED = true;
        }
    });

    auto description = general_desc();
    description.add(options_desc());
    description.add(positional_desc());
    po::positional_options_description positional;
    positional.add("genotypes", 1);
    po::variables_map vm;
    po::store(po::command_line_parser({arguments.begin() + 1, arguments.end()}).
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
    try {
        if (PREVIOUS_RESULT != "") {
            GradientDescent gradient_descent(GENOTYPES_FILE, MAX_SITES, CONCURRENCY);
            GridSearch grid_search(GENOTYPES_FILE);
            grid_search.read_results(PREVIOUS_RESULT);
            const std::string outdir = make_outdir();
            wtl::Pushd cd(outdir);
            wtl::ozfstream ost("gradient_descent.tsv.gz");
            gradient_descent.run(ost, grid_search.mle_params());
            std::cout << *gradient_descent.mle_params() << std::endl;
        } else if (GENOTYPES_FILE == "-") {
            GridSearch searcher(std::cin, MAX_SITES, CONCURRENCY);
            searcher.run(false);
        } else {
            GridSearch searcher(GENOTYPES_FILE, MAX_SITES, CONCURRENCY);
            // after constructor success
            const std::string outdir = make_outdir();
            wtl::Pushd cd(outdir);
            searcher.run(true);
        }
    } catch (const wtl::KeyboardInterrupt& e) {
        std::cerr << e.what() << std::endl;
    } catch (const lnpnan_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

std::string Program::make_outdir() const {
    std::smatch mobj;
    std::regex_search(GENOTYPES_FILE, mobj, std::regex("([^/]+?)\\.[^/]+$"));
    std::ostringstream oss;
    oss << mobj.str(1) << "-s" << MAX_SITES;
    const std::string outdir = oss.str();
    wtl::mkdir(outdir);
    return outdir;
}

} // namespace likeligrid
