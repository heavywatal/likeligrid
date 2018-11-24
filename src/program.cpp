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
#include <wtl/zlib.hpp>
#include <wtl/filesystem.hpp>
#include <clippson/clippson.hpp>

#include <regex>

namespace likeligrid {

namespace fs = wtl::filesystem;

//! Global variables mapper of commane-line arguments
nlohmann::json VM;

//! Options description for general purpose
inline clipp::group general_options(nlohmann::json* vm) {
    return (
      wtl::option(vm, {"h", "help"}, false, "print this help"),
      wtl::option(vm, {"version"}, false, "print version"),
      wtl::option(vm, {"v", "verbose"}, false, "verbose output"),
      wtl::option(vm, {"test"}, false, "run tests")
    ).doc("General:");
}

inline clipp::group program_options(nlohmann::json* vm, Program* prog) {
    return (
      wtl::option(vm, {"j", "parallel"}, &prog->CONCURRENCY),
      wtl::option(vm, {"s", "max-sites"}, &prog->MAX_SITES),
      wtl::option(vm, {"g", "gradient"}, &prog->GRADIENT_MODE),
      wtl::option(vm, {"e", "epistasis"}, &prog->EPISTASIS_PAIR),
      wtl::option(vm, {"p", "pleiotropy"}, &prog->PLEIOTROPY)
    ).doc("Program:");
}

Program::Program(const std::vector<std::string>& arguments) {HERE;
    wtl::join(arguments, std::cout, " ") << std::endl;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);
    wtl::set_SIGINT_handler();

    VM.clear();
    nlohmann::json vm_local;
    auto cli = (
      general_options(&vm_local),
      program_options(&VM, this),
      wtl::value<std::string>(&VM, "infile", &INFILE)
    );
    wtl::parse(cli, arguments);
    if (vm_local["help"]) {
        auto fmt = wtl::doc_format();
        std::cout << "Usage: " << PROJECT_NAME << " [options]\n\n";
        std::cout << clipp::documentation(cli, fmt) << "\n";
        throw wtl::ExitSuccess();
    }
    if (vm_local["version"]) {
        std::cout << PROJECT_VERSION << "\n";
        throw wtl::ExitSuccess();
    }
    WTL_ASSERT(EPISTASIS_PAIR.size() == 2u);
    if (vm_local["verbose"]) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << VM.dump(2) << std::endl;
    }
    if (vm_local["test"]) {
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
