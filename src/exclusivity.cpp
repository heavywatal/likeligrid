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

ExclusivityModel::ExclusivityModel(std::istream& infile, const size_t g, const size_t n):
    names_(wtl::read_header(infile)),
    genotypes_(wtl::eigen::read_array<size_t>(infile, names_.size())),
    grid_density_(g),
    max_results_(n) {}


inline double calc_denom(
    const Eigen::ArrayXd& weights,
    const Eigen::ArrayXd& exclusi,
    const size_t num_mutations) {

    if (num_mutations < 2) return 1.0;
    std::vector<size_t> indices(weights.size());
    std::iota(std::begin(indices), std::end(indices), 0);
    const std::vector<std::vector<size_t>> columns(num_mutations, indices);
    auto iter = wtl::itertools::product(columns);
    std::vector<double> probs;
    probs.reserve(static_cast<size_t>(iter.count_max()));
    for (const auto& v: iter()) {
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

void ExclusivityModel::run(const std::string& outfile) {HERE;
    results_.clear();
    const size_t nsam = genotypes_.rows();
    const ArrayXu vs = genotypes_.rowwise().sum();
    const size_t max_sites = vs.maxCoeff();
    const ArrayXu freqs = genotypes_.colwise().sum();
    const Eigen::ArrayXd weights = freqs.cast<double>() / freqs.sum();
    double lnp_const = (freqs.cast<double>() * weights.log()).sum();

    const Eigen::ArrayXd dups = genotypes_.unaryExpr([](size_t x){
        if (x > 0) {return --x;} else {return x;}
    }).colwise().sum().cast<double>();

    std::vector<size_t> s_counts(max_sites + 1, 0);
    for (size_t i=0; i<nsam; ++i) {
        ++s_counts[vs[i]];
        auto v = wtl::eigen::as_valarray(genotypes_.row(i));
        lnp_const += std::log(wtl::polynomial(v));
    }

    const double step = 1.0 / grid_density_;
    const Eigen::ArrayXd axis = Eigen::VectorXd::LinSpaced(grid_density_, 1.0, step).array();
    const auto columns = std::vector<Eigen::ArrayXd>(genotypes_.cols(), axis);
    auto iter = wtl::itertools::product(columns);
    const auto num_gridpoints = iter.count_max();
    for (const auto& params: iter()) {
        if (iter.count() % 1000 == 0) {  // snapshot for long run
            wtl::Fout fout(outfile);
            fout << "# " << iter.count() << " in " << num_gridpoints << "\n";
            write_results(fout);
        }

        double loglik = lnp_const;
        loglik += (dups * params.log()).sum();
        for (size_t s=0; s<=max_sites; ++s) {
            loglik -= s_counts[s] * std::log(calc_denom(weights, params, s));
        }

        results_.emplace(loglik, wtl::eigen::as_vector(params));
        while (results_.size() > max_results_) {
            results_.erase(results_.begin());
        }
    }
    wtl::Fout fout(outfile);
    write_results(fout);
}

std::ostream& ExclusivityModel::write_genotypes(std::ostream& ost, const bool header) const {
    if (header) {ost << wtl::join(names_, "\t") << "\n";}
    return ost << genotypes_.format(wtl::eigen::tsv());
}

std::ostream& ExclusivityModel::write_results(std::ostream& ost, const bool header) const {
    if (header) {ost << "loglik\t" << wtl::join(names_, "\t") << "\n";}
    for (const auto& p: results_) {
        ost << p.first << "\t" << wtl::join(p.second, "\t") << "\n";
    }
    return ost;
}

void ExclusivityModel::unit_test() {HERE;
    std::string geno = "a\tb\n0\t0\n0\t1\n1\t0\n1\t1\n";
    std::istringstream iss(geno);
    ExclusivityModel model(iss, 10);
    model.write_genotypes(std::cout);
    model.run();
}

} // namespace likeligrid
