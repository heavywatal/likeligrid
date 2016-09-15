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
    return bmath::cdf(bmath::complement(bmath::normal(mean, sd), x));
}

inline double Model::likelihood(double score) {
    score += epsilon_;
    return normal_ccdf(threshold_, score, score);
}

void Model::genotype2score() {HERE;
}

void Model::run() {HERE;
    wtl::Fin fin("dummy_genotype_cancer.tsv");
    auto genotype = wtl::eigen::read_matrix<double>(fin);
    std::cout << genotype.format(wtl::eigen::tsv()) << std::endl;
    const size_t ncol = genotype.cols();

    Eigen::VectorXd vxd = Eigen::VectorXd::LinSpaced(5, 0.0, 1.0);
    std::vector<Eigen::VectorXd> columns(ncol, vxd);
    auto sim = wtl::itertools::simplex(columns);

    for (const auto& coefs: sim()) {
        auto scores = (genotype * coefs).array();
        auto loglik = scores.unaryExpr([this](double x) {
            return likelihood(x);
        }).log().sum();
        std::cout << loglik << std::endl;
    }
    genotype2score();
}

void Model::unit_test() {HERE;
    Model model(0.5, 0.1);
    const std::vector<double> scores{0.1, 0.2, 0.4, 0.5, 0.9};
    for (const auto x: scores) {
        std::cout << model.likelihood(x) << std::endl;
    }
}

} // namespace lmpp
