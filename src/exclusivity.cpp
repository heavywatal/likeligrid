// -*- mode: c++; coding: utf-8 -*-
/*! @file exclusivity.cpp
    @brief Inplementation of Exclusivity class
*/
#include "exclusivity.hpp"

#include <functional>
#include <unordered_map>
#include <csignal>

#include <boost/dynamic_bitset.hpp>

#include <cxxwtils/debug.hpp>
#include <cxxwtils/exception.hpp>
#include <cxxwtils/iostr.hpp>
#include <cxxwtils/numeric.hpp>
#include <cxxwtils/math.hpp>
#include <cxxwtils/eigen.hpp>
#include <cxxwtils/itertools.hpp>
#include <cxxwtils/gz.hpp>

namespace {
    bool SIGINT_RAISED = false;
}

namespace likeligrid {

const std::vector<double> ExclusivityModel::STEPS_ = {0.4, 0.2, 0.1, 0.05, 0.02, 0.01};
const std::vector<size_t> ExclusivityModel::BREAKS_ = {5, 5, 5, 6, 5, 5};

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
    std::signal(SIGINT, [](int signum){
        if (signum == SIGINT) SIGINT_RAISED = true;
        // TODO: how to handle other signum?
    });
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
            run_impl(name_outfile("uniaxis"), wtl::itertools::uniaxis(axes_, mle_params_));
        } else {
            //TODO
            const auto axes = make_vicinity(mle_params_, 21, 1.0);
            const auto vicinity = make_vicinity(mle_params_, 3, 0.01, 1.2);
            run_impl(name_outfile("prototype"), wtl::itertools::prototype(axes, vicinity));
        }
        return;
    }
    const std::string outfile = name_outfile(infile);
    if (read_results(outfile) && start_ == 0) {
        std::cerr << outfile << " is already completed:" << std::endl;
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
    if (mle_params_.size() > 0 || read_results(infile)) {
        std::cerr << mle_params_.transpose() << std::endl;
        if (start_ > 0) {
            throw std::runtime_error("infile must be a complete result");
        }
        axes_ = make_vicinity(mle_params_, BREAKS_.at(stage_), STEPS_.at(stage_), 2.0);
        ++stage_;
    } else {
        const double step = STEPS_[0];
        const size_t breaks = BREAKS_[0];
        axes_.assign(names_.size(), Eigen::VectorXd::LinSpaced(breaks, 2.0, step).array());
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
    const auto max_count = gen.max_count();
    auto buffer = wtl::make_oss();
    std::ios::openmode mode = std::ios::out;
    if (start_ == 0) {
        buffer << "##max_count=" << max_count << "\n";
        buffer << "##max_sites=" << nsam_with_s_.size() - 1 << "\n";
        buffer << "##step=" << STEPS_.at(stage_) << "\n";
        buffer << "loglik\t" << wtl::join(names_, "\t") << "\n";
    } else {
        mode |= std::ios::app;
    }
    std::cerr << "\nWriting to " << outfile << std::endl;
    wtl::ogzstream fout(outfile, mode);
    double max_ll = std::numeric_limits<double>::lowest();
    for (const auto& params: gen(start_)) {
        if (gen.count() % 10000 == 0) {  // snapshot for long run
            std::cerr << "\r" << gen.count() << " in " << max_count << std::flush;
            fout << buffer.str();
            fout.strict_sync();
            buffer.str("");
        }

        double loglik = 0.0;
        loglik += (a_pathway_ * params.log()).sum();
        for (size_t s=0; s<nsam_with_s_.size(); ++s) {
            loglik -= nsam_with_s_[s] * std::log(calc_denom(w_pathway_, params, s));
        }

        buffer << (loglik += lnp_const_) << "\t"
               << wtl::join(wtl::eigen::vector(params), "\t") << "\n";
        if (loglik > max_ll) {
            max_ll = loglik;
            mle_params_ = params;
        }
        if (SIGINT_RAISED) {throw wtl::ExitSuccess("KeyboardInterrupt");}
    }
    fout << buffer.str();
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

bool ExclusivityModel::read_results(const std::string& infile) {HERE;
    wtl::igzstream ist(infile);
    if (!ist || infile == "/dev/null") return false;
    std::cerr << "Reading: " << infile << std::endl;
    const size_t max_count = read_metadata(ist);
    start_ = read_body(ist);
    if (start_ == max_count) start_ = 0;
    return true;
}

size_t ExclusivityModel::read_metadata(std::istream& ist) {HERE;
    std::string buffer;
    std::getline(ist, buffer);
    const size_t max_count = std::stoul(wtl::split(buffer, "=")[1]);
    std::getline(ist, buffer);
    const size_t max_sites = std::stoul(wtl::split(buffer, "=")[1]);
    if (nsam_with_s_.size() - 1 != max_sites) {
        std::cerr << "Note: -s is different between arg and file\n";
    }
    std::getline(ist, buffer);
    const double step = std::stod(wtl::split(buffer, "=")[1]);
    auto pred = std::bind(wtl::equal<double>, std::placeholders::_1, step);
    const auto it = std::find_if(STEPS_.begin(), STEPS_.end(), pred);
    if (it == STEPS_.end()) throw std::runtime_error("invalid step size");
    stage_ = it - STEPS_.begin();
    return max_count;
}

size_t ExclusivityModel::read_body(std::istream& ist) {HERE;
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
    double max_ll = std::numeric_limits<double>::lowest();
    while (std::getline(ist, buffer)) {
        ++i;
        std::istringstream iss(buffer);
        std::istream_iterator<double> it(iss);
        if (*it > max_ll) {
            max_ll = *it;
            mle_params_ = wtl::eigen::ArrayX(std::vector<double>(++it, std::istream_iterator<double>()));
        }
    }
    return i;
}

void ExclusivityModel::unit_test() {HERE;
    std::string geno = "A\tB\n0\t1\n1\t0\n1\t1\n0\t2\n";
    std::istringstream iss(geno);
    ExclusivityModel model(iss);
    model.write_genotypes(std::cerr);
    model.run("/dev/null");
}

} // namespace likeligrid
