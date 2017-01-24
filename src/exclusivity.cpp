// -*- mode: c++; coding: utf-8 -*-
/*! @file exclusivity.cpp
    @brief Inplementation of Exclusivity class
*/
#include "exclusivity.hpp"

#include <unordered_map>

#include <boost/dynamic_bitset.hpp>

#include <cxxwtils/debug.hpp>
#include <cxxwtils/exception.hpp>
#include <cxxwtils/iostr.hpp>
#include <cxxwtils/numeric.hpp>
#include <cxxwtils/math.hpp>
#include <cxxwtils/eigen.hpp>
#include <cxxwtils/itertools.hpp>

namespace likeligrid {

ExclusivityModel::ExclusivityModel(std::istream& infile,
    const size_t grid_density,
    const size_t max_sites,
    const size_t max_results):
    names_(wtl::read_header(infile)),
    genotypes_(wtl::eigen::read_array<size_t>(infile, names_.size())),
    grid_density_(grid_density),
    max_results_(max_results) {
        const auto pred = genotypes_.rowwise().sum().array() < max_sites;
        genotypes_ = wtl::eigen::filter(genotypes_, pred);

        std::vector<size_t> indices(genotypes_.cols());
        std::iota(std::begin(indices), std::end(indices), 0);
        const size_t max_sites_real = genotypes_.rowwise().sum().maxCoeff();
        index_iters_.reserve(max_sites_real);
        for (size_t i=0; i<=max_sites_real; ++i) {
            index_iters_.emplace_back(std::vector<std::vector<size_t>>(i, indices));
        }
    }

double ExclusivityModel::calc_denom(
    const Eigen::ArrayXd& weights,
    const Eigen::ArrayXd& exclusi,
    const size_t num_mutations) {

    if (num_mutations < 2) return 1.0;
    auto& iter = index_iters_[num_mutations];
    double sum_prob = 0.0;
    boost::dynamic_bitset<> mutated(exclusi.size());
    for (const auto& v: iter()) {
        double p = 1.0;
        for (const auto x: v) {
            p *= weights[x];
            if (mutated.test(x)) p*= exclusi[x];
            mutated.set(x);
        }
        sum_prob += p;
        mutated.reset();
    }
    iter.reset();
    return sum_prob;
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

    if (outfile != "/dev/stdout") {
        std::ifstream fin(outfile);
        read_results(fin);
    }
    const double step = 1.0 / grid_density_;
    const Eigen::ArrayXd axis = Eigen::VectorXd::LinSpaced(grid_density_, 1.0, step).array();
    const auto columns = std::vector<Eigen::ArrayXd>(genotypes_.cols(), axis);
    auto iter = wtl::itertools::product(columns);
    const auto num_gridpoints = iter.count_max();
    for (const auto& params: iter(start_)) {
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

void ExclusivityModel::read_results(std::istream& ist) {
    if (ist.fail()) return;
    std::string buffer;
    ist >> buffer;
    if (buffer != "#") throw wtl::ExitSuccess("Completed already");
    ist >> start_;
    std::getline(ist, buffer); // in count_max()
    std::getline(ist, buffer); // header
    while (std::getline(ist, buffer)) {
        std::istringstream iss(buffer);
        std::istream_iterator<double> it(iss);
        results_.emplace(double(*it), std::vector<double>(++it, std::istream_iterator<double>()));
    }
}

void ExclusivityModel::unit_test() {HERE;
    std::string geno = "a\tb\n0\t0\n0\t1\n1\t0\n1\t1\n";
    std::istringstream iss(geno);
    ExclusivityModel model(iss, 5);
    model.write_genotypes(std::cout);
    model.run();
}

} // namespace likeligrid
