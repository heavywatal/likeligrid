// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.cpp
    @brief Inplementation of Simulation class
*/
#include "model.hpp"

#include <cstdlib>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/math/distributions/normal.hpp>

#include "Eigen/Core"

#include <cxxwtils/iostr.hpp>
#include <cxxwtils/getopt.hpp>
#include <cxxwtils/prandom.hpp>
#include <cxxwtils/os.hpp>
#include <cxxwtils/gz.hpp>

namespace lmpp {

namespace bmath = boost::math;
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

po::options_description Model::options_desc() {HERE;
    po::options_description description("Simulation");
    description.add_options()
        ("write,w", po::value<bool>(&WRITE_TO_FILES)->default_value(WRITE_TO_FILES)->implicit_value(true))
        ("out_dir,o", po::value<std::string>(&OUT_DIR)->default_value(OUT_DIR))
        ("seed", po::value<unsigned int>(&SEED)->default_value(SEED));
    return description;
}

po::options_description Model::positional_desc() {HERE;
    po::options_description description("Positional");
    description.add_options()
        ("epsilon", po::value(&epsilon_)->default_value(epsilon_))
        ("threshold", po::value(&threshold_)->default_value(threshold_));
    return description;
}

void Model::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "Usage: tumopp [options] nsam howmany\n" << std::endl;
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

//! Unit test for each class
inline void test(const int flg) {HERE;
    switch (flg) {
      case 0:
        break;
      case 1:
        throw wtl::ExitSuccess();
      default:
        throw std::runtime_error("Unknown argument for --test");
    }
}

Model::Model(const std::vector<std::string>& arguments) {HERE;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);
    COMMAND_ARGS = wtl::str_join(arguments, " ");
    OUT_DIR = wtl::strftime("tumopp_%Y%m%d_%H%M_") + std::to_string(::getpid());

    auto description = general_desc();
    description.add(options_desc());
    description.add(positional_desc());
    po::positional_options_description positional;
    positional.add("nsam", 1)
              .add("howmany", 1);
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

void Model::genotype2score() {
    size_t dimensions = 2;
    size_t num_samples = 4;
    Eigen::VectorXd parameter(dimensions);
    Eigen::MatrixXd genotype(num_samples, dimensions);
    parameter << 0.3, 0.4;
    genotype << 0, 0,
        1, 0,
        0, 1,
        1, 1;
    std::cout << genotype * parameter << std::endl;
}

double Model::normal_lccdf(double score, const double sd) {
    bmath::normal bm_normal(score += epsilon_, sd);
    return std::log(1.0 - bmath::cdf(bm_normal, threshold_));
}

Eigen::MatrixXd read_matrix(const std::string& path) {
    std::ifstream ifs(path.c_str());
    std::string buffer;
    std::vector<std::string> lines;
    while (std::getline(ifs, buffer, '\t')) {
        lines.push_back(buffer);
    }
}

void Model::run() {HERE;
    std::cout << "lmpp" << std::endl;
    genotype2score();
    std::cout << normal_lccdf(0.1, 0.5) << std::endl;
    std::cout << normal_lccdf(0.3, 0.5) << std::endl;
    std::cout << normal_lccdf(0.4, 0.5) << std::endl;
    std::cout << normal_lccdf(0.5, 0.5) << std::endl;
    std::cout << normal_lccdf(0.7, 0.5) << std::endl;
}

void Model::write() const {HERE;
    if (WRITE_TO_FILES) {
        derr("mkdir && cd to " << OUT_DIR << std::endl);
        fs::create_directory(OUT_DIR);
        wtl::cd(OUT_DIR);
        wtl::Fout{"program_options.conf"} << CONFIG_STRING;
        std::cerr << wtl::iso8601datetime() << std::endl;
    }
}

} // namespace tumopp
