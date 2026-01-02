/*! @file program.cpp
    @brief Implementation of Program class
*/
#include "version.hpp"
#include "program.hpp"
#include "genotype.hpp"
#include "gridsearch.hpp"
#include "gradient_descent.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <wtl/zlib.hpp>
#include <clippson/clippson.hpp>

#include <filesystem>
#include <regex>

namespace likeligrid {

namespace fs = std::filesystem;

//! Global variables mapper of command-line arguments
nlohmann::json VM;

//! Options description for general purpose
inline clipp::group general_options(nlohmann::json* vm) {
    return (
      clippson::option(vm, {"h", "help"}, false, "print this help"),
      clippson::option(vm, {"version"}, false, "print version"),
      clippson::option(vm, {"v", "verbose"}, false, "verbose output"),
      clippson::option(vm, {"test"}, false, "run tests")
    ).doc("General:");
}

inline clipp::group program_options(nlohmann::json* vm) {
    std::vector<size_t> EPISTASIS_PAIR = {0u, 0u};
    return (
      clippson::option(vm, {"j", "parallel"}, 1),
      clippson::option(vm, {"s", "max-sites"}, 3u),
      clippson::option(vm, {"g", "gradient"}, false),
      clippson::option(vm, {"e", "epistasis"}, EPISTASIS_PAIR),
      clippson::option(vm, {"p", "pleiotropy"}, false)
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
      program_options(&VM),
      clippson::value<std::string>(&VM, "infile")
    );
    clippson::parse(cli, arguments);
    if (vm_local["help"]) {
        auto fmt = clippson::doc_format();
        std::cout << "Usage: " << PROJECT_NAME << " [options] infile\n\n";
        std::cout << clipp::documentation(cli, fmt) << "\n";
        throw wtl::ExitSuccess();
    }
    if (vm_local["version"]) {
        std::cout << PROJECT_VERSION << "\n";
        throw wtl::ExitSuccess();
    }
    WTL_ASSERT(VM.at("epistasis").size() == 2u);
    if (vm_local["verbose"]) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << VM.dump(2) << std::endl;
    }
    if (vm_local["test"]) {
        std::string infile = VM.at("--")[0u];
        wtl::zlib::ifstream ist(infile);
        GenotypeModel model(ist, VM.at("max-sites"));
        model.benchmark(VM.at("parallel"));
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

inline std::string make_outdir(const std::string& prefix) {
    const std::vector<size_t> epistasis_pair = VM.at("epistasis");
    std::ostringstream oss;
    oss << prefix << "-s" << VM.at("max-sites");
    if (VM.at("gradient")) {oss << "-g";}
    if (epistasis_pair[0u] != epistasis_pair[1u]) {
        oss << "-e" << epistasis_pair[0u]
            <<  "x" << epistasis_pair[1u];
    }
    if (VM.at("pleiotropy")) {oss << "-p";}
    const std::string outdir = oss.str();
    fs::create_directory(outdir);
    return outdir;
}

void Program::run() {HERE;
    const int concurrency = VM.at("parallel");
    const unsigned max_sites = VM.at("max-sites");
    const bool pleiotropy = VM.at("pleiotropy");
    const std::string infile = VM.at("--")[0u];
    const std::pair<size_t, size_t> epistasis{VM.at("epistasis")[0u], VM.at("epistasis")[1u]};
    WTL_ASSERT(!pleiotropy || (epistasis.first != epistasis.second));
    try {
        if (VM.at("gradient")) {
            if (infile == "-") {
                GradientDescent searcher(std::cin, max_sites, epistasis, pleiotropy, concurrency);
                searcher.run(std::cout);
                return;
            }
            GradientDescent searcher(infile, max_sites, epistasis, pleiotropy, concurrency);
            const auto outdir = make_outdir(extract_prefix(infile));
            const auto outfile = fs::path(outdir) / searcher.outfile();
            std::cerr << "outfile: " << outfile << std::endl;
            wtl::zlib::ofstream ost(outfile.native());
            ost.precision(std::cout.precision());
            searcher.run(ost);
        } else if (infile == "-") {
            GridSearch searcher(std::cin, max_sites, epistasis, pleiotropy, concurrency);
            searcher.run(false);
        } else {
            GridSearch searcher(infile, max_sites, epistasis, pleiotropy, concurrency);
            // after constructor success
            const std::string outdir = make_outdir(extract_prefix(infile));
            fs::current_path(outdir);
            searcher.run(true);
        }
    } catch (const wtl::KeyboardInterrupt& e) {
        std::cerr << e.what() << std::endl;
    }
}

} // namespace likeligrid
