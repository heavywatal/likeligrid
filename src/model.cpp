// -*- mode: c++; coding: utf-8 -*-
/*! @file model.cpp
    @brief Inplementation of Model class
*/
#include "model.hpp"

#include <boost/math/distributions/normal.hpp>

#include <cxxwtils/debug.hpp>
#include <cxxwtils/iostr.hpp>
#include <cxxwtils/eigen.hpp>
#include <cxxwtils/itertools.hpp>

namespace lmpp {

namespace bmath = boost::math;

inline double normal_cdf(const double x, const double mean, const double sd) {
    return bmath::cdf(bmath::normal(mean, sd), x);
}

inline double normal_ccdf(const double x, const double mean, const double sd) {
    if (sd == 0.0) {return 0.0;}
    return bmath::cdf(bmath::complement(bmath::normal(mean, sd), x));
}

class IsoSigma {
public:
    IsoSigma(const double t, const double e):
        threshold_(t), epsilon_(e) {}
    double operator()(double x) const {
        x += epsilon_;
        return normal_ccdf(threshold_, x, x);
    }
private:
    const double threshold_;
    const double epsilon_;
};

class ConstSigma {
public:
    ConstSigma(const double t, const double e, const double s):
        threshold_(t), epsilon_(e), sd_(s) {}
    double operator()(double x) const {
        return normal_ccdf(threshold_, x += epsilon_, sd_);
    }
private:
    const double threshold_;
    const double epsilon_;
    const double sd_;
};

/////////1/////////2/////////3/////////4

Model::Model(std::istream& infile, const size_t g, const size_t n):
    genotypes_(wtl::eigen::read_matrix<double>(infile)),
    grid_density_(g), max_results_(n) {}

std::multimap<double, std::vector<double>>
Model::run(const double threshold, const double epsilon) {HERE;
    Eigen::VectorXd axis = Eigen::VectorXd::LinSpaced(grid_density_, 0.0, 1.0 - epsilon);
    columns_ = std::vector<Eigen::VectorXd>(genotypes_.cols(), axis);

    auto sim = wtl::itertools::simplex(columns_, 1.0 - epsilon);
    std::multimap<double, std::vector<double>> leaders;
    std::function<double(double)> calc_lik;
    if (true) {
        calc_lik = IsoSigma(threshold, epsilon);
    } else {
        calc_lik = ConstSigma(threshold, epsilon, 0.1);
    }
    for (const auto& coefs: sim()) {
        const double loglik = (genotypes_ * coefs).array().unaryExpr(calc_lik).log().sum();
        leaders.emplace(loglik, wtl::eigen::as_vector(coefs));
        while (leaders.size() > max_results_) {
            leaders.erase(leaders.begin());
        }
    }
    return leaders;
}

std::ostream& operator<<(std::ostream& ost, const std::multimap<double, std::vector<double>>&m) {
    for (const auto& p: m) {
        ost << p.first << "\t" << wtl::join(p.second, "\t") << "\n";
    }
    return ost;
}

void Model::unit_test() {HERE;
    std::string geno = "0\t0\n0\t1\n1\t0\n1\t1\n";
    std::istringstream iss(geno);
    Model model(iss, 10);
    std::cout << model.genotypes() << std::endl;
    std::cout << model.run(0.5, 0.1) << std::endl;

    auto calc_lik = IsoSigma(0.5, 0.1);
    const std::vector<double> scores{0.1, 0.2, 0.4, 0.5, 0.9};
    for (const auto x: scores) {
        std::cout << calc_lik(x) << std::endl;
    }
}

} // namespace lmpp
