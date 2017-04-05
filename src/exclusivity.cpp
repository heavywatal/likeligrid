// -*- mode: c++; coding: utf-8 -*-
/*! @file exclusivity.cpp
    @brief Inplementation of Exclusivity class
*/
#include "exclusivity.hpp"
#include "util.hpp"

#include <functional>
#include <unordered_map>

#include <boost/math/distributions/chi_squared.hpp>

#include <wtl/debug.hpp>
#include <wtl/exception.hpp>
#include <wtl/iostr.hpp>
#include <wtl/numeric.hpp>
#include <wtl/math.hpp>
#include <wtl/eigen.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/os.hpp>

namespace likeligrid {

const std::vector<double> ExclusivityModel::STEPS_ = {0.4, 0.2, 0.1, 0.05, 0.02, 0.01};
const std::vector<size_t> ExclusivityModel::BREAKS_ = {5, 5, 5, 5, 6, 5};
bool ExclusivityModel::SIGINT_RAISED_ = false;

ExclusivityModel::ExclusivityModel(const std::string& infile, const size_t max_sites) {HERE;
    wtl::izfstream ifs(infile);
    names_ = wtl::read_header(ifs);
    auto pathtypes = wtl::eigen::read_array<ArrayXXu::value_type>(ifs, names_.size());
    ifs.close();
    const ArrayXu raw_s_sample = pathtypes.rowwise().sum().array();
    nsam_with_s_.assign(raw_s_sample.maxCoeff() + 1, 0);
    for (Eigen::Index i=0; i<raw_s_sample.size(); ++i) {
        ++nsam_with_s_[raw_s_sample[i]];
    }
    std::cerr << "Original N_s: " << nsam_with_s_ << std::endl;
    if (max_sites + 1 < nsam_with_s_.size()) {
        nsam_with_s_.resize(max_sites + 1);
        std::cerr << "Filtered N_s: " << nsam_with_s_ << std::endl;
    } else {
        std::cerr << "Note: -s is too large" << std::endl;
    }
    while (nsam_with_s_.back() == 0) {
        nsam_with_s_.resize(nsam_with_s_.size() - 1);
    }
    const auto final_max_s = nsam_with_s_.size() - 1;
    pathtypes = wtl::eigen::filter(pathtypes, raw_s_sample <= final_max_s);

    const Eigen::ArrayXd s_pathway = pathtypes.colwise().sum().cast<double>();
    w_pathway_ = s_pathway / s_pathway.sum();
    a_pathway_ = pathtypes.unaryExpr([](size_t x){
        if (x > 0) {return --x;} else {return x;}
    }).colwise().sum().cast<double>();
    for (Eigen::Index i=0; i<pathtypes.rows(); ++i) {
        auto v = wtl::eigen::valarray(pathtypes.row(i));
        lnp_const_ += std::log(wtl::multinomial(v));
    }
    lnp_const_ += (s_pathway * w_pathway_.log()).sum();
    std::cerr << "s_pathway_: " << s_pathway.transpose() << std::endl;
    std::cerr << "w_pathway_: " << w_pathway_.transpose() << std::endl;
    std::cerr << "a_pathway_: " << a_pathway_.transpose() << std::endl;
    std::cerr << "lnp_const_: " << lnp_const_ << std::endl;
    if (std::isnan(lnp_const_)) throw lnpnan_error();

    index_axes_.reserve(nsam_with_s_.size());
    std::vector<size_t> indices(pathtypes.cols());
    std::iota(std::begin(indices), std::end(indices), 0);
    for (size_t s=0; s<=nsam_with_s_.size(); ++s) {
        index_axes_.emplace_back(s, indices);
    }
}

inline std::vector<Eigen::ArrayXd>
make_vicinity(const Eigen::ArrayXd& center, const size_t breaks, const double radius, const double max=2.0) {
    std::vector<Eigen::ArrayXd> axes;
    axes.reserve(center.size());
    for (const double x: wtl::eigen::vector(center)){
        Eigen::ArrayXd axis = Eigen::ArrayXd::LinSpaced(breaks, x + radius, x - radius);
        axis = (axis * 100.0).round() / 100.0;  // grid precision = 0.01
        axes.push_back(wtl::eigen::filter(axis, (0.0 < axis) * (axis < max)));
    }
    return axes;
}

void ExclusivityModel::run(const std::string& infile) {HERE;
    init_axes(infile);
    if (stage_ == STEPS_.size()) {
        std::cerr << "Done: step size = " << STEPS_.back() << std::endl;
        --stage_;
        axes_.assign(names_.size(), Eigen::ArrayXd::LinSpaced(200, 2.0, 0.01));
        run_impl("uniaxis.tsv.gz", wtl::itertools::uniaxis(axes_, mle_params_));
        search_limits();
        return;
    }
    std::string outfile = "/dev/stdout";
    if (infile != "/dev/null") {
        auto oss = wtl::make_oss(2, std::ios::fixed);
        oss << "grid-" << STEPS_.at(stage_) << ".tsv.gz";
        outfile = oss.str();
        read_results(outfile);
    }
    for (size_t j=0; j<names_.size(); ++j) {
        std::cerr << names_[j] << ": " << axes_[j].transpose() << std::endl;
    }
    run_impl(outfile, wtl::itertools::product(axes_));
    if (outfile != "/dev/stdout") {
        run(outfile);
    }
}

void ExclusivityModel::search_limits() const {HERE;
    for (const auto& p: find_intersections()) {
        std::cerr << p.first << ": " << p.second.transpose() << std::endl;
        const auto axes = make_vicinity(p.second, 5, 0.02, 2.0);
        run_impl("limit-" + p.first + ".tsv.gz", wtl::itertools::product(axes));
    }
}

std::unordered_map<std::string, Eigen::ArrayXd> ExclusivityModel::find_intersections() const {HERE;
    namespace bmath = boost::math;
    bmath::chi_squared_distribution<> chisq(1.0);
    const double step = 0.01;
    const double max_ll = calc_loglik(mle_params_);
    const double threshold = max_ll - 0.5 * bmath::quantile(bmath::complement(chisq, 0.05));
    std::unordered_map<std::string, Eigen::ArrayXd> intersections;
    for (Eigen::Index j=0; j<mle_params_.size(); ++j) {
        auto th_path = mle_params_;
        for (size_t i=0; i<200; ++i) {
            th_path[j] -= step;
            if (th_path[j] < step || calc_loglik(th_path) < threshold) {
                intersections.emplace(names_[j] + "_L", th_path);
                break;
            }
        }
        for (size_t i=0; i<200; ++i) {
            th_path[j] += step;
            if (th_path[j] >= 2.0 || calc_loglik(th_path) < threshold) {
                intersections.emplace(names_[j] + "_U", th_path);
                break;
            }
        }
    }
    return intersections;
}

void ExclusivityModel::init_axes(const std::string& infile) {HERE;
    if (read_results(infile)) {
        if (skip_ > 0) {
            throw std::runtime_error("infile must be a complete result");
        }
        std::cerr << "mle_params_: " << mle_params_.transpose() << std::endl;
        axes_ = make_vicinity(mle_params_, BREAKS_.at(stage_), 2.0 * STEPS_.at(stage_), 2.0);
        ++stage_;
    } else {
        const double step = STEPS_[0];
        const size_t breaks = BREAKS_[0];
        axes_.assign(names_.size(), Eigen::VectorXd::LinSpaced(breaks, 2.0, step).array());
    }
}

void ExclusivityModel::run_impl(const std::string& outfile, wtl::itertools::Generator<Eigen::ArrayXd>&& gen) const {HERE;
    auto buffer = wtl::make_oss();
    std::ios::openmode mode = std::ios::out;
    if (wtl::exists(outfile)) {
        if (skip_ == 0 && outfile != "/dev/stdout") {
            std::cerr << outfile << " exists." << std::endl;
            return;
        }
        mode |= std::ios::app;
    } else {
        buffer << "##max_count=" << gen.max_count() << "\n";
        buffer << "##max_sites=" << nsam_with_s_.size() - 1 << "\n";
        buffer << "##step=" << STEPS_.at(stage_) << "\n";
        buffer << "loglik\t" << wtl::join(names_, "\t") << "\n";
    }
    std::cerr << "Writing: " << outfile << std::endl;
    std::cerr << skip_ << " to " << gen.max_count() << std::endl;
    wtl::ozfstream fout(outfile, mode);
    for (const auto& th_path: gen(skip_)) {
        buffer << calc_loglik(th_path) << "\t"
               << th_path.transpose().format(wtl::eigen::tsv()) << "\n";
        if (gen.count() % 10000 == 0) {  // snapshot for long run
            std::cerr << "*" << std::flush;
            fout << buffer.str();
            fout.strict_sync();
            buffer.str("");
        }
        if (SIGINT_RAISED_) {throw wtl::ExitSuccess("KeyboardInterrupt");}
    }
    std::cerr << "\n";
    fout << buffer.str();
}

double ExclusivityModel::calc_loglik(const Eigen::ArrayXd& th_path) const {
    double loglik = (a_pathway_ * th_path.log()).sum();
    // D = 1.0 when s < 2
    for (size_t s=2; s<nsam_with_s_.size(); ++s) {
        loglik -= nsam_with_s_[s] * std::log(calc_denom(w_pathway_, th_path, s));
    }
    return loglik += lnp_const_;
}

double ExclusivityModel::calc_denom(
    const Eigen::ArrayXd& w_pathway,
    const Eigen::ArrayXd& x_pathway,
    const size_t num_mutations) const {

    if (num_mutations < 2) return 1.0;
    auto iter = wtl::itertools::product(index_axes_[num_mutations]);
    double sum_prob = 0.0;
    bits_t bits(x_pathway.size());

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
    return sum_prob;
}

bool ExclusivityModel::read_results(const std::string& infile) {HERE;
    if (infile == "/dev/null")
        return false;
    try {
        wtl::izfstream ist(infile);
        std::cerr << "Reading: " << infile << std::endl;
        const size_t max_count = read_metadata(ist);
        skip_ = read_body(ist);
        if (skip_ == max_count) skip_ = 0;
        return true;
    } catch (std::ios::failure& e) {
        if (errno != 2) throw;
        return false;
    }
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
    stage_ = guess_stage(STEPS_, step);
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
    ExclusivityModel model("test.tsv");
    model.run("/dev/null");
}

} // namespace likeligrid
