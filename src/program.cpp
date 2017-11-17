/*! @file program.cpp
    @brief Implementation of Program class
*/
#include "version.hpp"
#include "program.hpp"
#include "pathtype.hpp"
#include "genotype.hpp"
#include "gridsearch.hpp"
#include "gradient_descent.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/getopt.hpp>
#include <wtl/zfstream.hpp>

#include <boost/filesystem.hpp>

#include <regex>

namespace likeligrid {

namespace fs = boost::filesystem;
namespace po = boost::program_options;

inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::bool_switch(), "print this help")
        ("verbose,v", po::bool_switch(), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1));
    return description;
}

po::options_description Program::options_desc() {HERE;
    po::options_description description("Program");
    description.add_options()
        ("parallel,j", po::value(&CONCURRENCY)->default_value(CONCURRENCY))
        ("max-sites,s", po::value(&MAX_SITES)->default_value(MAX_SITES))
        ("gradient,g", po::bool_switch(&GRADIENT_MODE))
        ("epistasis,e", po::value(&EPISTASIS_PAIR)->default_value(EPISTASIS_PAIR)->multitoken())
        ("pleiotropy,p", po::bool_switch(&PLEIOTROPY));
    return description;
}

po::options_description Program::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("infile", po::value(&INFILE)->default_value(INFILE));
    return description;
}

[[noreturn]] void Program::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "commit " << GIT_COMMIT_HASH
              << " [" << GIT_BRANCH << "]\n"
              << "Date:  " << GIT_COMMIT_TIME << std::endl;    std::cout << "Usage: tek [options]\n" << std::endl;
    std::cout << "Usage: likeligrid [options] infile\n" << std::endl;
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

//! Unit test for each class
void Program::test(const int flag) {HERE;
    switch (flag) {
      case 0:
        break;
      case 1:
        GridSearch::test();
        GradientDescent::test();
        throw wtl::ExitSuccess();
      case 2:
        GenotypeModel::test();
        PathtypeModel::test();
        throw wtl::ExitSuccess();
      case 3: {
        wtl::izfstream ist(INFILE);
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
    wtl::set_SIGINT_handler();

    auto description = general_desc();
    description.add(options_desc());
    description.add(positional_desc());
    po::positional_options_description positional;
    positional.add("infile", 1);
    po::variables_map vm;
    po::store(po::command_line_parser({arguments.begin() + 1, arguments.end()}).
              options(description).
              positional(positional).run(), vm);
    if (vm["help"].as<bool>()) {help_and_exit();}
    po::notify(vm);
    if (EPISTASIS_PAIR.size() != 2u) {
        throw std::runtime_error("EPISTASIS_PAIR.size() != 2U");
    }

    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << wtl::flags_into_string(vm) << std::endl;
    }
    test(vm["test"].as<int>());
}

void Program::run() {HERE;
    std::pair<size_t, size_t> epistasis{EPISTASIS_PAIR[0u], EPISTASIS_PAIR[1u]};
    if (PLEIOTROPY && (epistasis.first == epistasis.second)) {
        throw std::runtime_error("PLEIOTROPY && (epistasis.first == epistasis.second)");
    }
    try {
        if (GRADIENT_MODE) {
            std::valarray<double> mle_params;
            std::ostringstream outfile;
            if (wtl::endswith(INFILE, ".tsv.gz")) {// previous result
                std::string previous = fs::path(INFILE).filename().string();
                wtl::izfstream ist(INFILE);
                size_t max_sites;
                std::tie(INFILE, max_sites, std::ignore, std::ignore) = read_metadata(ist);
                std::tie(std::ignore, std::ignore, mle_params) = read_body(ist);
                outfile << "grad-from-s" << max_sites << "-" << previous;
            } else {
                outfile << "grad-from-center.tsv.gz";
            }
            GradientDescent gradient_descent(INFILE, MAX_SITES, epistasis, PLEIOTROPY, CONCURRENCY);
            const std::string outdir = make_outdir();
            fs::current_path(outdir);
            wtl::ozfstream ost(outfile.str());
            ost.precision(std::cout.precision());
            std::cerr << "outfile: " << ost.path() << std::endl;
            gradient_descent.run(ost, mle_params);
        } else if (INFILE == "-") {
            GridSearch searcher(std::cin, MAX_SITES, epistasis, PLEIOTROPY, CONCURRENCY);
            searcher.run(false);
        } else {
            GridSearch searcher(INFILE, MAX_SITES, epistasis, PLEIOTROPY, CONCURRENCY);
            // after constructor success
            const std::string outdir = make_outdir();
            fs::current_path(outdir);
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
    std::regex_search(INFILE, mobj, std::regex("([^/]+?)\\.[^/]+$"));
    std::ostringstream oss;
    oss << mobj.str(1) << "-s" << MAX_SITES;
    if (GRADIENT_MODE) {oss << "-g";}
    if (EPISTASIS_PAIR[0u] != EPISTASIS_PAIR[1u]) {
        oss << "-e" << EPISTASIS_PAIR[0u]
            <<  "x" << EPISTASIS_PAIR[1u];
    }
    const std::string outdir = oss.str();
    fs::create_directory(outdir);
    return outdir;
}

} // namespace likeligrid
