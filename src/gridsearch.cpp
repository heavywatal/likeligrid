/*! @file gridsearch.cpp
    @brief Implementation of GridSearch class
*/
#include "gridsearch.hpp"
#include "util.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/numeric.hpp>
#include <wtl/math.hpp>
#include <wtl/concurrent.hpp>

#include <chrono>

namespace likeligrid {

void GridSearch::init(const std::pair<size_t, size_t>& epistasis_pair, const bool pleiotropy) {HERE;
    model_.set_epistasis(epistasis_pair, pleiotropy);
    mle_params_.resize(model_.names().size());
    mle_params_ = 1.0;
}

void GridSearch::run(const bool writing) {HERE;
    while (stage_ < STEPS.size()) {
        if (writing) {run_fout();} else {run_cout();}
    }
    --stage_;
    search_limits();
}

void GridSearch::run_fout() {HERE;
    const std::string outfile = init_meta();
    std::cerr << "mle_params_: " << mle_params_ << std::endl;
    if (outfile.empty()) return;
    const auto axes = make_vicinity(mle_params_, BREAKS.at(stage_), radius(stage_));
    for (size_t j=0u; j<model_.names().size(); ++j) {
        std::cerr << model_.names()[j] << ": " << axes[j] << std::endl;
    }
    {
        wtl::ozfstream fout(outfile, std::ios::out | std::ios::app);
        std::cerr << "Writing: " << outfile << std::endl;
        run_impl(fout, wtl::itertools::product(axes));
    }
}

void GridSearch::run_cout() {HERE;
    const auto axes = make_vicinity(mle_params_, BREAKS.at(stage_), radius(stage_));
    for (size_t j=0u; j<model_.names().size(); ++j) {
        std::cerr << model_.names()[j] << ": " << axes[j] << std::endl;
    }
    {
        std::stringstream sst;
        run_impl(sst, wtl::itertools::product(axes));
        std::cout << sst.str();
        read_results(sst);
    }
    std::cerr << "mle_params_: " << mle_params_ << std::endl;
    ++stage_;
}

void GridSearch::search_limits() {HERE;
    namespace bmath = boost::math;
    bmath::chi_squared_distribution<> chisq(1.0);
    const double diff95 = 0.5 * bmath::quantile(bmath::complement(chisq, 0.05));
    auto axis = wtl::round(wtl::lin_spaced(200, 2.0, 0.01), 100);
    axis = (axis * 100.0).apply(std::round) / 100.0;
    std::map<std::string, std::valarray<double>> intersections;
    for (size_t i=0u; i<model_.names().size(); ++i) {
        const std::string outfile = "uniaxis-" + model_.names()[i] + ".tsv.gz";
        std::cerr << outfile << std::endl;
        std::stringstream sst;
        run_impl(sst, wtl::itertools::uniaxis(axis, mle_params_, i));
        wtl::ozfstream(outfile) << sst.str();
        const auto logliks = read_loglik(sst, axis.size());
        const double threshold = logliks.max() - diff95;
        const std::valarray<double> range = axis[logliks > threshold];
        auto bound_params = mle_params_;
        bound_params[i] = std::max(range.min() - 0.01, 0.01);
        intersections.emplace(model_.names()[i] + "_L", bound_params);
        bound_params[i] = std::min(range.max() + 0.01, 2.00);
        intersections.emplace(model_.names()[i] + "_U", bound_params);
    }
    for (const auto& p: intersections) {
        const std::string outfile = "limit-" + p.first + ".tsv.gz";
        std::cerr << outfile << ": " << p.second << std::endl;
        const auto axes = make_vicinity(p.second, 5u, 0.02);
        wtl::ozfstream fout(outfile);
        //TODO: if exists
        run_impl(fout, wtl::itertools::product(axes));
    }
}

void GridSearch::run_impl(std::ostream& ost, wtl::itertools::Generator<std::valarray<double>>&& gen) {HERE;
    std::cerr << skip_ << " to " << gen.max_count() << std::endl;
    if (skip_ == 0u) {
        write_header(ost, gen.max_count());
    }

    auto task = [this](const std::valarray<double> th_path) {
        // argument and model are copied for each thread
        auto buffer = wtl::make_oss();
        auto model_copy = this->model_;
        buffer << model_copy.calc_loglik(th_path) << "\t"
               << wtl::str_join(th_path, "\t") << "\n";
        return buffer.str();
    };

    wtl::ThreadPool pool(concurrency_);
    std::vector<std::future<std::string>> futures;
    futures.reserve(gen.max_count() - skip_);
    for (const auto& th_path: gen(skip_)) {
        futures.push_back(pool.submit(task, th_path));
    }

    auto buffer = wtl::make_oss();
    size_t stars = 0u;
    size_t i = skip_;
    const auto min_interval = std::chrono::seconds(1);
    auto next_time = std::chrono::system_clock::now();
    for (auto& ftr: futures) {
        buffer << ftr.get();
        ++i;
        auto now = std::chrono::system_clock::now();
        if (now > next_time || &ftr == &futures.back()) {
            next_time = now + min_interval;
            ost << buffer.str();
            buffer.str("");
            buffer.clear();
            for (size_t n= static_cast<size_t>(20.0 * i / gen.max_count()); stars<n; ++stars) {
                std::cerr << "*";
            }
        }
        if (wtl::SIGINT_RAISED()) {throw wtl::KeyboardInterrupt();}
    }
    std::cerr << "\n";
}

std::string GridSearch::init_meta() {HERE;
    if (stage_ >= STEPS.size()) return "";
    auto oss = wtl::make_oss(2u, std::ios::fixed);
    oss << "grid-" << STEPS.at(stage_) << ".tsv.gz";
    std::string outfile = oss.str();
    try {
        wtl::izfstream ist(outfile);
        std::cerr << "Reading: " << outfile << std::endl;
        read_results(ist);
        if (skip_ == 0u) {
            ++stage_;
            outfile = init_meta();
        }
    } catch (std::ios::failure& e) {
        if (errno != 2u) throw;
    }
    return outfile;
}

void GridSearch::read_results(std::istream& ist) {HERE;
    size_t max_count;
    double step;
    std::tie(std::ignore, std::ignore, max_count, step) = read_metadata(ist);
    stage_ = guess_stage(step);
    std::vector<std::string> colnames;
    std::valarray<double> mle_params;
    std::tie(skip_, colnames, mle_params) = read_body(ist);
    if (skip_ == max_count) {  // is complete file
        skip_ = 0u;
        mle_params_.swap(mle_params);
    }
}

void GridSearch::read_results(const std::string& infile) {
    wtl::izfstream ist(infile);
    read_results(ist);
}

void GridSearch::write_header(std::ostream& ost, const size_t max_count) const {
    ost << "##genotype_file=" << model_.filename() << "\n";
    ost << "##max_sites=" << model_.max_sites() << "\n";
    ost << "##max_count=" << max_count << "\n";
    ost << "##step=" << STEPS.at(stage_) << "\n";
    ost << "loglik\t" << wtl::join(model_.names(), "\t") << "\n";
}

void GridSearch::test() {HERE;
    std::stringstream sst;
    sst <<
R"({
  "pathway": ["A", "B"],
  "annotation": ["0011", "1100"],
  "sample": ["0011", "0101", "1001", "0110", "1010", "1100"]
})";
    GridSearch searcher(sst, 4u, {0, 1});
    searcher.run_cout();
}

} // namespace likeligrid
