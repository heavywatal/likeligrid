// -*- mode: c++; coding: utf-8 -*-
/*! @file exclusivity.cpp
    @brief Inplementation of Exclusivity class
*/
#include "exclusivity.hpp"

#include <functional>
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

const std::vector<double> ExclusivityModel::STEPS_ = {0.2, 0.1, 0.05, 0.02, 0.01};
const std::vector<size_t> ExclusivityModel::BREAKS_ = {5, 5, 6, 5, 5};

ExclusivityModel::ExclusivityModel(std::istream& genotypes,
    const size_t max_sites,
    const size_t max_results):
    names_(wtl::read_header(genotypes)),
    genotypes_(wtl::eigen::read_array<size_t>(genotypes, names_.size())),
    max_results_(max_results) {

    const auto pred = genotypes_.rowwise().sum().array() <= max_sites;
    genotypes_ = wtl::eigen::filter(genotypes_, pred);

    max_sites_ = genotypes_.rowwise().sum().maxCoeff();
    index_iters_.reserve(max_sites_);
    std::vector<size_t> indices(genotypes_.cols());
    std::iota(std::begin(indices), std::end(indices), 0);
    for (size_t i=0; i<=max_sites_; ++i) {
        index_iters_.emplace_back(std::vector<std::vector<size_t>>(i, indices));
    }
}

void ExclusivityModel::run(const std::string& infile) {HERE;
    init_axes(infile);
    if (stage_ == STEPS_.size()) {
        std::cerr << "Done: step size = " << STEPS_.back() << std::endl;
        return;
    }
    const std::string outfile = name_outfile(infile);
    if (read_results(outfile) && start_ == 0) {
        std::cerr << outfile << " is already completed:" << std::endl;
        write_results(std::cout);
        run(outfile);
        return;
    }
    std::cerr << "Start: " << start_ << std::endl;
    for (size_t i=0; i<names_.size(); ++i) {
        std::cerr << names_[i] << ": " << axes_[i].transpose() << std::endl;
    }
    run_impl(outfile);
    if (outfile != "/dev/null") {
        run(outfile);
    }
}

void ExclusivityModel::init_axes(const std::string& infile) {HERE;
    if (read_results(infile)) {
        if (start_ > 0) {
            throw std::runtime_error("infile must be a complete result");
        }
        max_results_ = results_.size();
        const double step = STEPS_.at(stage_);
        const size_t breaks = BREAKS_.at(stage_);
        axes_.clear();
        axes_.reserve(names_.size());
        for (const double x: best_result()) {
            Eigen::ArrayXd axis = Eigen::ArrayXd::LinSpaced(breaks, x + step, x - step);
            axes_.push_back(wtl::eigen::filter(axis, axis > 0.0));
        }
        results_.clear();
        ++stage_;
    } else {
        const double step = STEPS_[0];
        const size_t breaks = BREAKS_[0];
        axes_.assign(names_.size(), Eigen::VectorXd::LinSpaced(breaks, 1.0, step).array());
    }
}

std::string ExclusivityModel::name_outfile(const std::string& infile) const {HERE;
    std::string prefix = "output";
    if (infile == "/dev/null") {
        return infile;
    } else if (!infile.empty()) {
        prefix = wtl::split(infile, "-")[0];
    }
    std::ostringstream oss(prefix, std::ios::ate);
    oss << "-s" << max_sites_
        << "-step" << STEPS_.at(stage_)
        << ".tsv";
    return oss.str();
}

void ExclusivityModel::run_impl(const std::string& outfile) {HERE;
    const ArrayXu freqs = genotypes_.colwise().sum();
    const Eigen::ArrayXd weights = freqs.cast<double>() / freqs.sum();
    double lnp_const = (freqs.cast<double>() * weights.log()).sum();

    const Eigen::ArrayXd dups = genotypes_.unaryExpr([](size_t x){
        if (x > 0) {return --x;} else {return x;}
    }).colwise().sum().cast<double>();

    const ArrayXu s_samples = genotypes_.rowwise().sum();
    std::vector<size_t> s_counts(max_sites_ + 1, 0);
    for (Eigen::Index i=0; i<genotypes_.rows(); ++i) {
        ++s_counts[s_samples[i]];
        auto v = wtl::eigen::valarray(genotypes_.row(i));
        lnp_const += std::log(wtl::polynomial(v));
    }

    auto iter = wtl::itertools::product(axes_);
    const auto max_count = iter.max_count();
    for (const auto& params: iter(start_)) {
        if (iter.count() % 1000 == 0) {  // snapshot for long run
            std::cerr << "\r" << iter.count() << " in " << max_count << std::flush;
            wtl::Fout fout(outfile);
            fout << "# " << iter.count() << " in " << max_count << "\n";
            write_results(fout);
        }

        double loglik = lnp_const;
        loglik += (dups * params.log()).sum();
        for (size_t s=0; s<=max_sites_; ++s) {
            loglik -= s_counts[s] * std::log(calc_denom(weights, params, s));
        }

        results_.emplace(loglik, wtl::eigen::vector(params));
        while (results_.size() > max_results_) {
            results_.erase(results_.begin());
        }
    }
    write_results(std::cout);
    std::cerr << "\nWriting to " << outfile << std::endl;
    wtl::Fout fout(outfile);
    write_results(fout);
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

std::ostream& ExclusivityModel::write_genotypes(std::ostream& ost, const bool header) const {HERE;
    if (header) {ost << wtl::join(names_, "\t") << "\n";}
    return ost << genotypes_.format(wtl::eigen::tsv());
}

std::ostream& ExclusivityModel::write_results(std::ostream& ost) const {
    ost << "##max_sites=" << max_sites_ << "\n";
    ost << "##step=" << STEPS_.at(stage_) << "\n";
    ost << "loglik\t" << wtl::join(names_, "\t") << "\n";
    for (const auto& p: results_) {
        ost << p.first << "\t" << wtl::join(p.second, "\t") << "\n";
    }
    return ost;
}

bool ExclusivityModel::read_results(const std::string& infile) {HERE;
    results_.clear();
    std::ifstream ist(infile);
    if (ist.fail() || ist.bad() || infile == "/dev/null") return false;
    std::cerr << "Reading: " << infile << std::endl;
    read_metadata(ist);
    read_body(ist);
    return true;
}

void ExclusivityModel::read_metadata(std::istream& ist) {HERE;
    std::string buffer;
    ist >> buffer;
    if (buffer == "#") { // incomplete
        ist >> start_;
        std::getline(ist, buffer); // in count_max()
        ist >> buffer;
    } else {
        start_ = 0;
    }
    max_sites_ = std::stoul(wtl::split(buffer, "=")[1]);
    ist >> buffer;
    const double step = std::stod(wtl::split(buffer, "=")[1]);
    auto pred = std::bind(wtl::equal<double>, std::placeholders::_1, step);
    const auto it = std::find_if(STEPS_.begin(), STEPS_.end(), pred);
    if (it == STEPS_.end()) throw std::runtime_error("invalid step size");
    stage_ = it - STEPS_.begin();
}

void ExclusivityModel::read_body(std::istream& ist) {HERE;
    std::string buffer;
    ist >> buffer; // loglik
    std::getline(ist, buffer); // header
    buffer.erase(0, 1); // \t
    if (names_ != wtl::split(buffer, "\t")) {
        std::cerr << names_ << std::endl;
        std::cerr << wtl::split(buffer, "\t") << std::endl;
        throw std::runtime_error("Column names are wrong");
    }
    while (std::getline(ist, buffer)) {
        std::istringstream iss(buffer);
        std::istream_iterator<double> it(iss);
        results_.emplace(double(*it), std::vector<double>(++it, std::istream_iterator<double>()));
    }
}

void ExclusivityModel::unit_test() {HERE;
    std::string geno = "A\tB\n0\t0\n0\t1\n1\t0\n1\t1\n";
    std::istringstream iss(geno);
    ExclusivityModel model(iss);
    model.write_genotypes(std::cerr);
    model.run("/dev/null");
}

} // namespace likeligrid
