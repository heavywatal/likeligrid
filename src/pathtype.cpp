// -*- mode: c++; coding: utf-8 -*-
/*! @file exclusivity.cpp
    @brief Implementation of Exclusivity class
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
#include <wtl/zfstream.hpp>
#include <wtl/os.hpp>

namespace likeligrid {

const std::vector<double> ExclusivityModel::STEPS_ = {0.4, 0.2, 0.1, 0.05, 0.02, 0.01};
const std::vector<size_t> ExclusivityModel::BREAKS_ = {5, 5, 5, 5, 6, 5};
bool ExclusivityModel::SIGINT_RAISED_ = false;

ExclusivityModel::ExclusivityModel(const std::string& infile, const size_t max_sites):
    ExclusivityModel(wtl::izfstream(infile), max_sites) {HERE;}

ExclusivityModel::ExclusivityModel(std::istream&& ist, const size_t max_sites) {HERE;
    names_ = wtl::read_header(ist);
    auto pathtypes = wtl::read_valarrays<uint>(ist);
    const auto raw_s_sample = wtl::row_sums(pathtypes);
    nsam_with_s_.assign(raw_s_sample.max() + 1, 0);
    for (const auto s: raw_s_sample) {
        ++nsam_with_s_[s];
    }
    std::cerr << "Original N_s: " << nsam_with_s_ << std::endl;
    if (max_sites + 1 < nsam_with_s_.size()) {
        nsam_with_s_.resize(max_sites + 1);
        std::cerr << "Filtered N_s: " << nsam_with_s_ << std::endl;
    } else {
        std::cerr << "Note: -s is too large" << std::endl;
    }
    while (nsam_with_s_.back() == 0) {
        nsam_with_s_.pop_back();
    }
    const auto final_max_s = nsam_with_s_.size() - 1;
    pathtypes = wtl::filter(pathtypes, raw_s_sample <= final_max_s);

    const auto s_pathway = wtl::cast<double>(wtl::col_sums(pathtypes));
    w_pathway_ = s_pathway / s_pathway.sum();
    auto duptypes = pathtypes;
    for (auto& row: duptypes) {
        for (auto& x: row) {
            if (x > 0) --x;
        }
    }
    a_pathway_ = wtl::cast<double>(wtl::col_sums(duptypes));
    for (size_t i=0; i<pathtypes.size(); ++i) {
        lnp_const_ += std::log(wtl::multinomial(pathtypes[i]));
    }
    lnp_const_ += (s_pathway * std::log(w_pathway_)).sum();
    std::cerr << "s_pathway_: " << s_pathway << std::endl;
    std::cerr << "w_pathway_: " << w_pathway_ << std::endl;
    std::cerr << "a_pathway_: " << a_pathway_ << std::endl;
    std::cerr << "lnp_const_: " << lnp_const_ << std::endl;
    if (std::isnan(lnp_const_)) throw lnpnan_error();

    index_axes_.reserve(nsam_with_s_.size());
    std::vector<size_t> indices(a_pathway_.size());
    std::iota(std::begin(indices), std::end(indices), 0);
    for (size_t s=0; s<=nsam_with_s_.size(); ++s) {
        index_axes_.emplace_back(s, indices);
    }
    mle_params_.resize(a_pathway_.size());
    mle_params_ = 1.2;
}

void ExclusivityModel::run(const std::string& infile) {HERE;
    const std::string outfile = init_meta(infile);
    std::cerr << "mle_params_: " << mle_params_ << std::endl;
    if (outfile == "") {
        std::cerr << "Done: step size = " << STEPS_.at(--stage_) << std::endl;
        search_limits();
        return;
    }
    const auto axes = make_vicinity(mle_params_, BREAKS_.at(stage_), 2.0 * STEPS_.at(stage_));
    for (size_t j=0; j<names_.size(); ++j) {
        std::cerr << names_[j] << ": " << axes[j] << std::endl;
    }
    {
        wtl::ozfstream fout(outfile, std::ios::out | std::ios::app);
        std::cerr << "Writing: " << outfile << std::endl;
        run_impl(fout, wtl::itertools::product(axes));
    }
    if (outfile != "/dev/stdout") {
        run(outfile);
    }
}

void ExclusivityModel::search_limits() const {HERE;
    namespace bmath = boost::math;
    bmath::chi_squared_distribution<> chisq(1.0);
    const double diff95 = 0.5 * bmath::quantile(bmath::complement(chisq, 0.05));
    auto axis = wtl::round(wtl::lin_spaced(200, 2.0, 0.01), 100);
    axis = (axis * 100.0).apply(std::round) / 100.0;
    std::map<std::string, std::valarray<double>> intersections;
    for (size_t i=0; i<names_.size(); ++i) {
        const std::string outfile = "uniaxis-" + names_[i] + ".tsv.gz";
        std::cerr << outfile << std::endl;
        std::stringstream sst;
        run_impl(sst, wtl::itertools::uniaxis(axis, mle_params_, i));
        wtl::ozfstream(outfile) << sst.str();
        const auto logliks = read_loglik(sst, axis.size());
        const double threshold = logliks.max() - diff95;
        const std::valarray<double> range = axis[logliks > threshold];
        auto bound_params = mle_params_;
        bound_params[i] = std::max(range.min() - 0.01, 0.01);
        intersections.emplace(names_[i] + "_L", bound_params);
        bound_params[i] = std::min(range.max() + 0.01, 2.00);
        intersections.emplace(names_[i] + "_U", bound_params);
    }
    for (const auto& p: intersections) {
        const std::string outfile = "limit-" + p.first + ".tsv.gz";
        std::cerr << outfile << ": " << p.second << std::endl;
        const auto axes = make_vicinity(p.second, 5, 0.02);
        wtl::ozfstream fout(outfile);
        //TODO: if exists
        run_impl(fout, wtl::itertools::product(axes));
    }
}

void ExclusivityModel::run_impl(std::ostream& ost, wtl::itertools::Generator<std::valarray<double>>&& gen) const {HERE;
    std::cerr << skip_ << " to " << gen.max_count() << std::endl;
    if (skip_ == 0) {
        ost << "##max_count=" << gen.max_count() << "\n";
        ost << "##max_sites=" << nsam_with_s_.size() - 1 << "\n";
        ost << "##step=" << STEPS_.at(stage_) << "\n";
        ost << "loglik\t" << wtl::join(names_, "\t") << "\n";
    }
    auto buffer = wtl::make_oss();
    for (const auto& th_path: gen(skip_)) {
        buffer << calc_loglik(th_path) << "\t"
               << wtl::str_join(th_path, "\t") << "\n";
        if (gen.count() % 10000 == 0) {  // snapshot for long run
            std::cerr << "*" << std::flush;
            ost << buffer.str();
            ost.flush();
            buffer.str("");
        }
        if (SIGINT_RAISED_) {throw wtl::KeyboardInterrupt();}
    }
    std::cerr << "-\n";
    ost << buffer.str();
}

double ExclusivityModel::calc_loglik(const std::valarray<double>& th_path) const {
    const size_t max_sites = nsam_with_s_.size() - 1;
    double loglik = (a_pathway_ * std::log(th_path)).sum();
    // D = 1.0 when s < 2
    for (size_t s=2; s<=max_sites; ++s) {
        loglik -= nsam_with_s_[s] * std::log(calc_denom(w_pathway_, th_path, s));
    }
    return loglik += lnp_const_;
}

double ExclusivityModel::calc_denom(
    const std::valarray<double>& w_pathway,
    const std::valarray<double>& th_pathway,
    const size_t num_mutations) const {

    if (num_mutations < 2) return 1.0;
    auto iter = wtl::itertools::product(index_axes_[num_mutations]);
    double sum_prob = 0.0;
    bits_t bits(th_pathway.size());

    for (const auto& indices: iter()) {
        double p = 1.0;
        for (const auto j: indices) {
            p *= w_pathway[j];
            if (bits[j]) p *= th_pathway[j];
            bits.set(j);
        }
        sum_prob += p;
        bits.reset();
    }
    return sum_prob;
}

std::string ExclusivityModel::init_meta(const std::string& infile) {HERE;
    if (infile == "/dev/null") return "/dev/stdout";
    if (stage_ >= STEPS_.size()) return "";
    auto oss = wtl::make_oss(2, std::ios::fixed);
    oss << "grid-" << STEPS_.at(stage_) << ".tsv.gz";
    std::string outfile = oss.str();
    if (read_results(outfile) && skip_ == 0) {
        ++stage_;
        outfile = init_meta(outfile);
    }
    return outfile;
}

bool ExclusivityModel::read_results(const std::string& infile) {HERE;
    if (infile == "/dev/null")
        return false;
    try {
        wtl::izfstream ist(infile);
        std::cerr << "Reading: " << infile << std::endl;
        size_t max_count;
        double step;
        std::tie(max_count, std::ignore, step) = read_metadata(ist);
        stage_ = guess_stage(STEPS_, step);
        std::vector<std::string> colnames;
        std::valarray<double> mle_params;
        std::tie(skip_, colnames, mle_params) = read_body(ist);
        if (skip_ == max_count) {  // is complete file
            skip_ = 0;
            mle_params_.swap(mle_params);
        }
        if (names_ != colnames) {
            std::ostringstream oss;
            oss << "Contradiction in column names:\n"
                << "genotype file: " << names_ << "\n"
                << "result file:" << colnames;
            throw std::runtime_error(oss.str());
        }
        return true;
    } catch (std::ios::failure& e) {
        if (errno != 2) throw;
        return false;
    }
}

void ExclusivityModel::unit_test() {HERE;
    std::stringstream sst;
    sst <<
R"(A B
0 1
1 0
1 1
0 2
)";
    ExclusivityModel model(std::move(sst), 3);
    model.run("/dev/null");
}

} // namespace likeligrid
