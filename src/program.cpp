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
#include <wtl/chrono.hpp>
#include <wtl/getopt.hpp>
#include <wtl/zlib.hpp>
#include <wtl/filesystem.hpp>

#include <regex>

namespace likeligrid {

namespace fs = wtl::filesystem;
namespace po = boost::program_options;

inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::bool_switch(), "print this help")
        ("version", po::bool_switch(), "print version")
        ("verbose,v", po::bool_switch(), "verbose output")
        ("test", po::bool_switch());
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
    std::cout << "Usage: " << PROJECT_NAME << " [options] infile\n\n";
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

Program::Program(const std::vector<std::string>& arguments) {HERE;
    wtl::join(arguments, std::cout, " ") << std::endl;
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
    if (vm["version"].as<bool>()) {
        std::cout << PROJECT_VERSION << "\n";
        throw wtl::ExitSuccess();
    }
    po::notify(vm);
    WTL_ASSERT(EPISTASIS_PAIR.size() == 2u);

    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << wtl::flags_into_string(vm) << std::endl;
    }
    if (vm["test"].as<bool>()) {
        wtl::zlib::ifstream ist(INFILE);
        GenotypeModel model(ist, MAX_SITES);
        model.benchmark(CONCURRENCY);
        throw wtl::ExitSuccess();
    }
}

inline std::string extract_prefix(const std::string& infile) {
    fs::path inpath(infile);
    for (fs::path p=inpath.filename(); !p.extension().empty(); p=p.stem()) {
        if (p.extension() == ".json") {
            return p.stem().string();
        }
    }
    std::smatch mobj;
    std::string dir = inpath.parent_path();
    if (std::regex_search(dir, mobj, std::regex("([\\w-]+)-s\\d"))) {
        return mobj[1];
    }
    throw std::runtime_error("Cannot extract prefix: " + infile);
}

void Program::run() {HERE;
    std::pair<size_t, size_t> epistasis{EPISTASIS_PAIR[0u], EPISTASIS_PAIR[1u]};
    WTL_ASSERT(!PLEIOTROPY || (epistasis.first != epistasis.second));
    try {
        if (GRADIENT_MODE) {
            if (INFILE == "-") {
                GradientDescent searcher(std::cin, MAX_SITES, epistasis, PLEIOTROPY, CONCURRENCY);
                searcher.run(std::cout);
                return;
            }
            GradientDescent searcher(INFILE, MAX_SITES, epistasis, PLEIOTROPY, CONCURRENCY);
            const auto outdir = make_outdir(extract_prefix(INFILE));
            const auto outfile = fs::path(outdir) / searcher.outfile();
            std::cerr << "outfile: " << outfile << std::endl;
            wtl::zlib::ofstream ost(outfile.native());
            ost.precision(std::cout.precision());
            searcher.run(ost);
        } else if (INFILE == "-") {
            GridSearch searcher(std::cin, MAX_SITES, epistasis, PLEIOTROPY, CONCURRENCY);
            searcher.run(false);
        } else {
            GridSearch searcher(INFILE, MAX_SITES, epistasis, PLEIOTROPY, CONCURRENCY);
            // after constructor success
            const std::string outdir = make_outdir(extract_prefix(INFILE));
            fs::current_path(outdir);
            searcher.run(true);
        }
    } catch (const wtl::KeyboardInterrupt& e) {
        std::cerr << e.what() << std::endl;
    }
}

std::string Program::make_outdir(const std::string& prefix) const {
    std::ostringstream oss;
    oss << prefix << "-s" << MAX_SITES;
    if (GRADIENT_MODE) {oss << "-g";}
    if (EPISTASIS_PAIR[0u] != EPISTASIS_PAIR[1u]) {
        oss << "-e" << EPISTASIS_PAIR[0u]
            <<  "x" << EPISTASIS_PAIR[1u];
    }
    if (PLEIOTROPY) {oss << "-p";}
    const std::string outdir = oss.str();
    fs::create_directory(outdir);
    return outdir;
}

} // namespace likeligrid
