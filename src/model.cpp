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

class IsoVar {
public:
    IsoVar(const double t, const double e):
        threshold_(t), epsilon_(e) {}
    double operator()(double x) const {
        x += epsilon_;
        return normal_ccdf(threshold_, x, std::sqrt(x));
    }
private:
    const double threshold_;
    const double epsilon_;
};

class ConstVar {
public:
    ConstVar(const double t, const double e, const double s):
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
    names_(wtl::read_header(infile)),
    genotypes_(wtl::eigen::read_matrix<double>(infile, names_.size())),
    grid_density_(g), max_results_(n) {}

void Model::run(const double threshold, const double epsilon) {HERE;
    results_.clear();
    Eigen::VectorXd axis = Eigen::VectorXd::LinSpaced(grid_density_, 0.0, 1.0 - epsilon);
    columns_ = std::vector<Eigen::VectorXd>(genotypes_.cols(), axis);

    auto sim = wtl::itertools::simplex(columns_, 1.0 - epsilon);
    std::function<double(double)> calc_lik;
    if (true) {
        calc_lik = IsoVar(threshold, epsilon);
    } else {
        calc_lik = ConstVar(threshold, epsilon, 0.1);
    }
    for (const auto& coefs: sim()) {
        const double loglik = (genotypes_ * coefs).array().unaryExpr(calc_lik).log().sum();
        results_.emplace(loglik, wtl::eigen::as_vector(coefs));
        while (results_.size() > max_results_) {
            results_.erase(results_.begin());
        }
    }
}

std::ostream& Model::write_genotypes(std::ostream& ost, const bool header) const {
    if (header) {ost << wtl::join(names_, "\t") << "\n";}
    return ost << genotypes_.format(wtl::eigen::tsv());
}

std::ostream& Model::write_results(std::ostream& ost, const bool header) const {
    if (header) {ost << "loglik\t" << wtl::join(names_, "\t") << "\n";}
    for (const auto& p: results_) {
        ost << p.first << "\t" << wtl::join(p.second, "\t") << "\n";
    }
    return ost;
}

void Model::unit_test() {HERE;
    std::string geno = "a\tb\n0\t0\n0\t1\n1\t0\n1\t1\n";
    std::istringstream iss(geno);
    Model model(iss, 10);
    model.write_genotypes(std::cout);
    model.run(0.5, 0.1);
    model.write_results(std::cout);

    auto calc_lik = IsoVar(0.5, 0.1);
    const std::vector<double> scores{0.1, 0.2, 0.4, 0.5, 0.9};
    for (const auto x: scores) {
        std::cout << calc_lik(x) << std::endl;
    }
}

} // namespace lmpp
