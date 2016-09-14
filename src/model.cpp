// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.cpp
    @brief Inplementation of Simulation class
*/
#include "model.hpp"

#include <cstdlib>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/coroutine2/coroutine.hpp>

#define EIGEN_NO_DEBUG
#include "Eigen/Core"

#include <cxxwtils/iostr.hpp>
#include <cxxwtils/getopt.hpp>
#include <cxxwtils/prandom.hpp>
#include <cxxwtils/algorithm.hpp>
#include <cxxwtils/os.hpp>
#include <cxxwtils/gz.hpp>
#include <cxxwtils/eigen.hpp>

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

inline double normal_cdf(const double x, const double mean, const double sd, const bool complement=false) {
    if (complement) {
        return bmath::cdf(bmath::complement(bmath::normal(mean, sd), x));
    } else {
        return bmath::cdf(bmath::normal(mean, sd), x);
    }
}

inline double Model::normal_ccdf(double score, const double sd) {
    return normal_cdf(threshold_, score += epsilon_, sd, true);
}

inline bool equals_one(double x) {
    return std::fabs(x - 1.0) < std::numeric_limits<double>::epsilon();
}

template <class value_type>
class BruteForce {
  public:
    typedef boost::coroutines2::coroutine<value_type> coro_t;
    BruteForce() = delete;
    BruteForce(const value_type& axis, const size_t dimensions):
        dimensions_(dimensions),
        level_(dimensions),
        current_(dimensions),
        axis_(axis) {}

    typename coro_t::pull_type generate() {
        return typename coro_t::pull_type([this](typename coro_t::push_type& yield){source(yield);});
    }

  private:
    void source(typename coro_t::push_type& yield) {
        if (--level_ > 0) {
            for (const auto x: axis_) {
                current_[level_] = x;
                source(yield);
            }
            ++level_;
        } else {
            for (const auto x: axis_) {
                current_[level_] = x;
                yield(value_type(current_));
            }
            ++level_;
        }
    }

    const double dimensions_;
    size_t level_;
    value_type current_;
    value_type axis_;
};

template <class T>
BruteForce<T> brute_force(const T& axis, const size_t dimensions) {
    return BruteForce<T>(axis, dimensions);
}

class SimplexIterator;

class Simplex {
  public:
    typedef std::vector<double> value_type;
    typedef SimplexIterator iterator;
    typedef boost::coroutines2::coroutine<value_type> coro_t;
    Simplex() = delete;
    Simplex(const value_type& axis, const size_t dimensions): bf_(axis, dimensions) {}

    iterator begin();
    iterator end();
    Simplex& operator++() {
        return *this;
    }
    typename coro_t::pull_type generate() {
        return typename coro_t::pull_type([this](typename coro_t::push_type& yield){source(yield);});
    }

  private:
    void source(typename coro_t::push_type& yield) {
        const double e = std::numeric_limits<double>::epsilon();
        for (const auto& v: bf_.generate()) {
            if (std::fabs(wtl::sum(v) - 1.0) < e) {yield(v);}
        }
    }
    BruteForce<value_type> bf_;
};

class SimplexIterator: public std::iterator<std::forward_iterator_tag, double> {
  public:
    SimplexIterator(Simplex* x): ptr_(x) {}
    Simplex& operator*() const {
        return *ptr_;
    }
    Simplex* operator->() const {
        return ptr_;
    }
    SimplexIterator& operator++() {
        ptr_->operator++();
        return *this;
    }
    bool operator!= (const SimplexIterator& it) const {
        return (ptr_ != it.ptr_);
    }
  private:
    Simplex* ptr_;
};

SimplexIterator Simplex::begin() {return SimplexIterator(this);}
SimplexIterator Simplex::end() {return SimplexIterator(this);}

template <class T>
inline void nested_loop(const T& range, size_t nest=3, T current={}) {
    if (current.empty()) current.assign(nest, 0);
    if (--nest > 0) {
        for (const auto& x: range) {
            current[nest] = x;
            nested_loop(range, nest, current);
        }
    } else {
        for (const auto& x: range) {
            current[nest] = x;
            std::cout << current << std::endl;
        }
    }
}

template <>
inline void nested_loop<Eigen::VectorXd>(const Eigen::VectorXd& range, size_t nest, Eigen::VectorXd current) {
    if (current.size() == 0) current.resize(nest);
    auto end = range.data() + range.size();
    if (--nest > 0) {
        for (auto p=range.data(); p<end; ++p) {
            current[nest] = *p;
            nested_loop(range, nest, current);
        }
    } else {
        for (auto p=range.data(); p<end; ++p) {
            current[nest] = *p;
            std::cout << wtl::eigen::as_vector(current) << std::endl;
        }
    }
}

void Model::genotype2score() {HERE;
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
    std::cout << normal_ccdf(0.1, 0.5) << std::endl;
    std::cout << normal_ccdf(0.3, 0.5) << std::endl;
    std::cout << normal_ccdf(0.4, 0.5) << std::endl;
    std::cout << normal_ccdf(0.5, 0.5) << std::endl;
    std::cout << normal_ccdf(0.7, 0.5) << std::endl;
}

void Model::run() {HERE;
    // genotype2score();
    Eigen::VectorXd vxd = Eigen::VectorXd::LinSpaced(3, 0.0, 1.0);
    auto vd = wtl::eigen::as_vector(vxd);
    auto bf = brute_force(vd, 3);
    for (const auto& x: bf.generate()) {
        std::cout << x << std::endl;
    }

    Simplex sim(vd, 3);
    for (const auto& x: sim.generate()) {
        std::cout << x << std::endl;
    }
}

void Model::write() const {HERE;
    auto mat = Eigen::Matrix3d::Random();
    std::ofstream("test.tsv") << mat.format(wtl::eigen::tsv());
    std::cout << wtl::eigen::read_matrix<double>("test.tsv").format(wtl::eigen::tsv()) << std::endl;
    if (WRITE_TO_FILES) {
        derr("mkdir && cd to " << OUT_DIR << std::endl);
        fs::create_directory(OUT_DIR);
        wtl::cd(OUT_DIR);
        wtl::Fout{"program_options.conf"} << CONFIG_STRING;
        std::cerr << wtl::iso8601datetime() << std::endl;
    }
}

} // namespace tumopp
