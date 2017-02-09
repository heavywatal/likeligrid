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
#include <cxxwtils/gz.hpp>

namespace likeligrid {

const std::vector<double> ExclusivityModel::STEPS_ = {0.2, 0.1, 0.05, 0.02, 0.01};
const std::vector<size_t> ExclusivityModel::BREAKS_ = {5, 5, 6, 5, 5};

ExclusivityModel::ExclusivityModel(std::istream& genotypes,
    const size_t max_sites,
    const size_t max_results):
    names_(wtl::read_header(genotypes)),
    genotypes_(wtl::eigen::read_array<size_t>(genotypes, names_.size())),
    max_results_(max_results)
    {HERE;

    const ArrayXu raw_s_sample = genotypes_.rowwise().sum().array();
    nsam_with_s_.assign(raw_s_sample.maxCoeff() + 1, 0);
    for (Eigen::Index i=0; i<raw_s_sample.size(); ++i) {
        ++nsam_with_s_[raw_s_sample[i]];
    }
    std::cerr << "Original N_s: " << nsam_with_s_ << std::endl;
    if (nsam_with_s_.size() - 1 < max_sites) {
        std::cerr << "Note: -s is too large" << std::endl;
    } else {
        nsam_with_s_.resize(max_sites + 1);
        std::cerr << "Filtered N_s: " << nsam_with_s_ << std::endl;
    }
    const auto final_max_s = nsam_with_s_.size() - 1;
    genotypes_ = wtl::eigen::filter(genotypes_, raw_s_sample <= final_max_s);

    const Eigen::ArrayXd s_pathway = genotypes_.colwise().sum().cast<double>();
    w_pathway_ = s_pathway / s_pathway.sum();
    a_pathway_ = genotypes_.unaryExpr([](size_t x){
        if (x > 0) {return --x;} else {return x;}
    }).colwise().sum().cast<double>();
    lnp_const_ = (s_pathway * w_pathway_.log()).sum();
    for (Eigen::Index i=0; i<genotypes_.rows(); ++i) {
        auto v = wtl::eigen::valarray(genotypes_.row(i));
        lnp_const_ += std::log(wtl::polynomial(v));
    }
    std::cerr << "s_pathway_: " << s_pathway.transpose() << std::endl;
    std::cerr << "w_pathway_: " << w_pathway_.transpose() << std::endl;
    std::cerr << "a_pathway_: " << a_pathway_.transpose() << std::endl;
    std::cerr << "lnp_const_: " << lnp_const_ << std::endl;

    index_iters_.reserve(nsam_with_s_.size());
    std::vector<size_t> indices(genotypes_.cols());
    std::iota(std::begin(indices), std::end(indices), 0);
    for (size_t s=0; s<=nsam_with_s_.size(); ++s) {
        index_iters_.emplace_back(std::vector<std::vector<size_t>>(s, indices));
    }
}

inline std::vector<Eigen::ArrayXd>
make_vicinity(const Eigen::ArrayXd& center, const size_t breaks, const double radius, const double max=1.0) {
    std::vector<Eigen::ArrayXd> axes;
    axes.reserve(center.size());
    for (const double x: wtl::eigen::vector(center)){
        Eigen::ArrayXd axis = Eigen::ArrayXd::LinSpaced(breaks, x + radius, x - radius);
        axes.push_back(wtl::eigen::filter(axis, (0.0 < axis) * (axis < max)));
    }
    return axes;
}

void ExclusivityModel::run(const std::string& infile) {HERE;
    init_axes(infile);
    if (stage_ == STEPS_.size()) {
        std::cerr << "Done: step size = " << STEPS_.back() << std::endl;
        --stage_;
        max_results_ = -1;
        if (true) {
            axes_.assign(names_.size(), Eigen::ArrayXd::LinSpaced(100, 1.0, 0.01));
            run_impl(name_outfile("uniaxis"), wtl::itertools::uniaxis(axes_, best_));
        } else {
            //TODO
            const auto axes = make_vicinity(best_, 21, 1.0);
            const auto vicinity = make_vicinity(best_, 3, 0.01, 1.2);
            run_impl(name_outfile("prototype"), wtl::itertools::prototype(axes, vicinity));
        }
        return;
    }
    const std::string outfile = name_outfile(infile);
    if (read_results(outfile) && start_ == 0) {
        std::cerr << outfile << " is already completed:" << std::endl;
        write_results(std::cout, 10);
        run(outfile);
        return;
    }
    std::cerr << "Start: " << start_ << std::endl;
    for (size_t j=0; j<names_.size(); ++j) {
        std::cerr << names_[j] << ": " << axes_[j].transpose() << std::endl;
    }
    run_impl(outfile, wtl::itertools::product(axes_));
    if (outfile != "/dev/null") {
        run(outfile);
    }
}

void ExclusivityModel::init_axes(const std::string& infile) {HERE;
    if (!results_.empty() || read_results(infile, 1)) {
        write_results(std::cerr, 1);
        if (start_ > 0) {
            throw std::runtime_error("infile must be a complete result");
        }
        best_ = results_.begin()->second;
        results_.clear();
        axes_ = make_vicinity(best_, BREAKS_.at(stage_), STEPS_.at(stage_), 2.0);
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
    oss << "-s" << nsam_with_s_.size() - 1
        << "-step" << STEPS_.at(stage_)
        << ".tsv.gz";
    return oss.str();
}

void ExclusivityModel::run_impl(const std::string& outfile, wtl::itertools::Generator<Eigen::ArrayXd>&& gen) {HERE;
    std::cerr << "\nWriting to " << outfile << std::endl;
    const auto max_count = gen.max_count();
    for (const auto& params: gen(start_)) {
        if (gen.count() % 10000 == 0) {  // snapshot for long run
            std::cerr << "\r" << gen.count() << " in " << max_count << std::flush;
            auto oss = wtl::make_oss();
            oss << "# " << gen.count() << " in " << max_count << "\n";
            write_results(oss);
            wtl::gzip{wtl::Fout(outfile)} << oss.str();
        }

        double loglik = 0.0;
        loglik += (a_pathway_ * params.log()).sum();
        for (size_t s=0; s<nsam_with_s_.size(); ++s) {
            loglik -= nsam_with_s_[s] * std::log(calc_denom(w_pathway_, params, s));
        }

        results_.emplace(loglik += lnp_const_, params);
        while (results_.size() > max_results_) {
            results_.erase(--results_.end());
        }
    }
    write_results(std::cout, 20);
    auto oss = wtl::make_oss();
    write_results(oss);
    wtl::gzip{wtl::Fout(outfile)} << oss.str();
}

double ExclusivityModel::calc_denom(
    const Eigen::ArrayXd& w_pathway,
    const Eigen::ArrayXd& x_pathway,
    const size_t num_mutations) {

    if (num_mutations < 2) return 1.0;
    auto& iter = index_iters_[num_mutations];
    double sum_prob = 0.0;
    boost::dynamic_bitset<> bits(x_pathway.size());
    for (const auto& indices: iter()) {
        double p = 1.0;
        for (const auto j: indices) {
            p *= w_pathway[j];
            if (bits[j]) p *= x_pathway[j];
            bits.set(j);
        }
        sum_prob += p;
        bits.reset();
    }
    iter.reset();
    return sum_prob;
}

std::ostream& ExclusivityModel::write_genotypes(std::ostream& ost, const bool header) const {HERE;
    if (header) {ost << wtl::join(names_, "\t") << "\n";}
    return ost << genotypes_.format(wtl::eigen::tsv());
}

std::ostream& ExclusivityModel::write_results(std::ostream& ost, const size_t max_rows) const {
    ost << "##max_sites=" << nsam_with_s_.size() - 1 << "\n";
    ost << "##step=" << STEPS_.at(stage_) << "\n";
    ost << "loglik\t" << wtl::join(names_, "\t") << "\n";
    size_t i = 0;
    for (const auto& p: results_) {
        if (++i > max_rows) {
            ost << "# ... with " << results_.size() - max_rows << " more rows\n";
            break;
        }
        ost << p.first << "\t" << wtl::join(wtl::eigen::vector(p.second), "\t") << "\n";
    }
    return ost;
}

bool ExclusivityModel::read_results(const std::string& infile, const size_t max_rows) {HERE;
    results_.clear();
    std::ifstream ist(infile);
    if (ist.fail() || ist.bad() || infile == "/dev/null") return false;
    wtl::gunzip zist(ist);
    std::cerr << "Reading: " << infile << std::endl;
    read_metadata(zist);
    read_body(zist, max_rows);
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
    const size_t max_sites = std::stoul(wtl::split(buffer, "=")[1]);
    if (nsam_with_s_.size() - 1 != max_sites) {
        std::cerr << "Note: -s is different between arg and file\n";
    }
    ist >> buffer;
    const double step = std::stod(wtl::split(buffer, "=")[1]);
    auto pred = std::bind(wtl::equal<double>, std::placeholders::_1, step);
    const auto it = std::find_if(STEPS_.begin(), STEPS_.end(), pred);
    if (it == STEPS_.end()) throw std::runtime_error("invalid step size");
    stage_ = it - STEPS_.begin();
}

void ExclusivityModel::read_body(std::istream& ist, const size_t max_rows) {HERE;
    std::string buffer;
    ist >> buffer; // loglik
    std::getline(ist, buffer); // header
    buffer.erase(0, 1); // \t
    if (names_ != wtl::split(buffer, "\t")) {
        std::cerr << names_ << std::endl;
        std::cerr << wtl::split(buffer, "\t") << std::endl;
        throw std::runtime_error("Column names are wrong");
    }
    size_t i = 0;
    while (std::getline(ist, buffer)) {
        if (++i > max_rows) break;
        std::istringstream iss(buffer);
        std::istream_iterator<double> it(iss);
        results_.emplace(double(*it), wtl::eigen::ArrayX(std::vector<double>(++it, std::istream_iterator<double>())));
    }
}

void ExclusivityModel::unit_test() {HERE;
    std::string geno = "A\tB\n0\t1\n1\t0\n1\t1\n0\t2\n";
    std::istringstream iss(geno);
    ExclusivityModel model(iss);
    model.write_genotypes(std::cerr);
    model.run("/dev/null");
}

} // namespace likeligrid
