// -*- mode: c++; coding: utf-8 -*-
/*! @file exclusivity.cpp
    @brief Inplementation of Exclusivity class
*/
#include "exclusivity.hpp"

#include <unordered_map>

#include <cxxwtils/debug.hpp>
#include <cxxwtils/iostr.hpp>
#include <cxxwtils/numeric.hpp>
#include <cxxwtils/math.hpp>
#include <cxxwtils/eigen.hpp>
#include <cxxwtils/itertools.hpp>

namespace likeligrid {

class FuncObj {
public:
    FuncObj() = default;
    double operator()(double x) const {
        return x;
    }
};

/////////1/////////2/////////3/////////4

Exclusivity::Exclusivity(std::istream& infile, const size_t g):
    names_(wtl::read_header(infile)),
    genotypes_(wtl::eigen::read_matrix<double>(infile, names_.size())),
    grid_density_(g) {}


inline double calc_denom(
    const Eigen::VectorXd& weights,
    const Eigen::VectorXd& exclusi,
    const size_t num_mutations) {

    const size_t ncol = weights.size();
    std::vector<size_t> indices(ncol);
    std::iota(std::begin(indices), std::end(indices), 0);
    const std::vector<std::vector<size_t>> columns(num_mutations, indices);
    auto gen = wtl::itertools::product(columns);
    std::vector<double> probs;
    probs.reserve(static_cast<size_t>(gen.count_max()));
    for (const auto v: gen()) {
        std::unordered_set<size_t> mutated;
        double p = 1.0;
        for (const auto x: v) {
            p *= weights[x];
            if (!mutated.insert(x).second) {
                p *= exclusi[x];
            }
        }
        probs.push_back(p);
    }
    return wtl::sum(probs);
}


void Exclusivity::run(const std::string& outfile) {HERE;
    // Calculate denominator with all posibble genotypes for each S
    Eigen::VectorXd freqs = genotypes_.colwise().sum() / genotypes_.sum();
    Eigen::VectorXd alphas = Eigen::VectorXd::Ones(freqs.size()) * 0.6;
    for (size_t i=1; i<=5; ++i) {
        std::cout << i << ": " << calc_denom(freqs, alphas, i) << std::endl;;
    }

    // deprecated
    results_.clear();
    Eigen::VectorXd axis = Eigen::VectorXd::LinSpaced(grid_density_, 0.0, 1.0);
    std::vector<Eigen::VectorXd> columns_ = std::vector<Eigen::VectorXd>(genotypes_.cols(), axis);

    auto sim = wtl::itertools::simplex(columns_, 1.0);
    const auto num_gridpoints = sim.count_max();
    std::function<double(double)> calc_lik;
    calc_lik = FuncObj();
    for (const auto& coefs: sim()) {
        if (sim.count() % 1000 == 0) {
            wtl::Fout fout(outfile);
            fout << "# " << sim.count() << " in "
                 << static_cast<double>(sim.count_all()) / static_cast<double>(num_gridpoints)
                 << " (" << sim.count_all() << " / " << num_gridpoints << ")\n";
            write_results(fout);
        }
        const double loglik = (genotypes_ * coefs).array().unaryExpr(calc_lik).log().sum();
        results_.emplace(loglik, wtl::eigen::as_vector(coefs));
    }
    wtl::Fout fout(outfile);
    write_results(fout);
}

std::ostream& Exclusivity::write_genotypes(std::ostream& ost, const bool header) const {
    if (header) {ost << wtl::join(names_, "\t") << "\n";}
    return ost << genotypes_.format(wtl::eigen::tsv());
}

std::ostream& Exclusivity::write_results(std::ostream& ost, const bool header) const {
    if (header) {ost << "loglik\t" << wtl::join(names_, "\t") << "\n";}
    for (const auto& p: results_) {
        ost << p.first << "\t" << wtl::join(p.second, "\t") << "\n";
    }
    return ost;
}

void Exclusivity::unit_test() {HERE;
    std::string geno = "a\tb\n0\t0\n0\t1\n1\t0\n1\t1\n";
    std::istringstream iss(geno);
    Exclusivity exclusivity(iss, 10);
    exclusivity.write_genotypes(std::cout);
    exclusivity.run();
    exclusivity.write_results(std::cout);

    auto calc_lik = FuncObj();
    const std::vector<double> scores{0.1, 0.2, 0.4, 0.5, 0.9};
    for (const auto x: scores) {
        std::cout << calc_lik(x) << std::endl;
    }
}

} // namespace likeligrid
