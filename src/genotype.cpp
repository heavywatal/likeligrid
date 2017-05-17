// -*- mode: c++; coding: utf-8 -*-
/*! @file genotype.cpp
    @brief Implementation of GenotypeModel class
*/
#include "genotype.hpp"
#include "util.hpp"

#include <functional>
#include <chrono>

#include <json.hpp>

#include <wtl/debug.hpp>
#include <wtl/exception.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/algorithm.hpp>
#include <wtl/numeric.hpp>
#include <wtl/math.hpp>
#include <wtl/os.hpp>
#include <wtl/concurrent.hpp>

namespace likeligrid {

bool GenotypeModel::SIGINT_RAISED_ = false;

GenotypeModel::GenotypeModel(const std::string& infile,
    const size_t max_sites,
    const unsigned int concurrency)
    : GenotypeModel(wtl::izfstream(infile), max_sites, concurrency) {HERE;}

GenotypeModel::GenotypeModel(
    std::istream& ist,
    const size_t max_sites,
    const unsigned int concurrency)
    : concurrency_(concurrency) {HERE;

    nlohmann::json jso;
    ist >> jso;
    names_ = jso["pathway"].get<std::vector<std::string>>();
    const size_t npath = names_.size();
    annot_.reserve(npath);
    for (const std::string& s: jso["annotation"]) {
        annot_.emplace_back(s);
    }
    std::cerr << "annot_: " << annot_ << std::endl;

    const size_t nsam = jso["sample"].size();
    std::vector<bits_t> all_genotypes;
    all_genotypes.reserve(nsam);
    for (const std::string& s: jso["sample"]) {
        all_genotypes.emplace_back(s);
    }

    genot_.reserve(nsam);
    const size_t ngene = jso["sample"].at(0).get<std::string>().size();
    nsam_with_s_.assign(ngene + 1, 0);  // at most
    std::valarray<double> s_gene(ngene);
    for (const auto& bits: all_genotypes) {
        const size_t s = bits.count();
        ++nsam_with_s_[s];
        if (s > max_sites) continue;
        genot_.push_back(bits);
        for (size_t j=0; j<ngene; ++j) {
            if (bits[j]) ++s_gene[j];
        }
    }
    wtl::rstrip(&nsam_with_s_);
    std::cerr << "Original N_s: " << nsam_with_s_ << std::endl;
    if (max_sites + 1 < nsam_with_s_.size()) {
        nsam_with_s_.resize(max_sites + 1);
        std::cerr << "Using N_s: " << nsam_with_s_ << std::endl;
    } else {
        std::cerr << "Note: -s is too large" << std::endl;
    }
    w_gene_ = s_gene / s_gene.sum();
    std::cerr << "s_gene : " << s_gene << std::endl;
    std::cerr << "w_gene_: " << w_gene_ << std::endl;

    mle_params_.resize(annot_.size());
    mle_params_ = 1.0;
}

void GenotypeModel::run(const bool writing) {HERE;
    while (stage_ < STEPS.size()) {
        if (writing) {run_fout();} else {run_cout();}
    }
    --stage_;
    search_limits();
}

void GenotypeModel::run_fout() {HERE;
    const std::string outfile = init_meta();
    std::cerr << "mle_params_: " << mle_params_ << std::endl;
    if (outfile.empty()) return;
    const auto axes = make_vicinity(mle_params_, BREAKS.at(stage_), radius(stage_));
    for (size_t j=0; j<names_.size(); ++j) {
        std::cerr << names_[j] << ": " << axes[j] << std::endl;
    }
    {
        wtl::ozfstream fout(outfile, std::ios::out | std::ios::app);
        std::cerr << "Writing: " << outfile << std::endl;
        run_impl(fout, wtl::itertools::product(axes));
    }
}

void GenotypeModel::run_cout() {HERE;
    const auto axes = make_vicinity(mle_params_, BREAKS.at(stage_), radius(stage_));
    for (size_t j=0; j<names_.size(); ++j) {
        std::cerr << names_[j] << ": " << axes[j] << std::endl;
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

void GenotypeModel::search_limits() const {HERE;
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

void GenotypeModel::run_impl(std::ostream& ost, wtl::itertools::Generator<std::valarray<double>>&& gen) const {HERE;
    std::cerr << skip_ << " to " << gen.max_count() << std::endl;
    if (skip_ == 0) {
        ost << "##max_count=" << gen.max_count() << "\n";
        ost << "##max_sites=" << nsam_with_s_.size() - 1 << "\n";
        ost << "##step=" << STEPS.at(stage_) << "\n";
        ost << "loglik\t" << wtl::join(names_, "\t") << "\n";
        ost.flush();
    }

    wtl::Semaphore semaphore(concurrency_);
    auto task = [this,&semaphore](const std::valarray<double> th_path) {
        // argument is copied for thread
        auto buffer = wtl::make_oss();
        buffer << calc_loglik(th_path) << "\t"
               << wtl::str_join(th_path, "\t") << "\n";
        semaphore.unlock();
        return buffer.str();
    };

    std::deque<std::future<std::string>> futures;
    const auto min_interval = std::chrono::seconds(2);
    auto next_time = std::chrono::system_clock::now() + min_interval;
    size_t stars = 0;
    for (const auto& th_path: gen(skip_)) {
        semaphore.lock();
        futures.push_back(std::async(std::launch::async, task, th_path));
        auto now = std::chrono::system_clock::now();
        if (now > next_time) {
            next_time = now + min_interval;
            size_t progress = 0;
            while (!futures.empty() && wtl::is_ready(futures.front())) {
                ost << futures.front().get();
                futures.pop_front();
                ++progress;
            }
            if (progress > 0) {
                ost.flush();
                for (size_t n= 0.2 * gen.percent(); stars<n; ++stars) {
                    std::cerr << "*";
                }
            }
        }
        if (SIGINT_RAISED_) {throw wtl::KeyboardInterrupt();}
    }
    for (auto& ftr: futures) {
        ost << ftr.get();
    }
    for (; stars<20U; ++stars) {
        std::cerr << "*";
    }
    std::cerr << std::endl;
}

inline std::valarray<size_t> to_indices(const bits_t& bits) {
    std::valarray<size_t> indices(bits.count());
    for (size_t i=0, j=0; i<indices.size(); ++j) {
        if (bits[j]) {
            indices[i] = j;
            ++i;
        }
    }
    return indices;
}

inline double slice_prod(const std::valarray<double>& coefs, const bits_t& bits) {
    double p = 1.0;
    for (size_t i=0; i<coefs.size(); ++i) {
        if (bits[i]) p *= coefs[i];
    }
    return p;
}

class Denoms {
  public:
    Denoms() = delete;
    Denoms(const std::valarray<double>& w_gene,
        const std::valarray<double>& th_path,
        const std::vector<bits_t>& annot,
        const size_t max_sites):
        w_gene_(w_gene),
        th_path_(th_path),
        annot_(annot),
        max_sites_(max_sites),
        denoms_(max_sites + 1)
    {
        const size_t ngene = w_gene.size();
        effects_.reserve(ngene);
        for (size_t j=0; j<ngene; ++j) {
            effects_.emplace_back(translate(j));
        }
        // std::cerr << "effects_: " << effects_ << std::endl;
        mutate(bits_t(), bits_t(), 1.0);
        // std::cerr << "denoms_: " << denoms_ << std::endl;
    }
    const std::valarray<double> log() const {return std::log(denoms_);}

    double lnp_sample(const bits_t& genotype) const {
        double p = 0.0;
        const double p_basic = slice_prod(w_gene_, genotype);
        auto mut_route = to_indices(genotype);
        do {
            p += p_basic * discount(mut_route);
        } while (std::next_permutation(std::begin(mut_route), std::end(mut_route)));
        return std::log(p);
    }

  private:
    void mutate(const bits_t& genotype, const bits_t& pathtype, const double anc_p) {
        const auto s = genotype.count() + 1;
        for (size_t j=0; j<w_gene_.size(); ++j) {
            if (genotype[j]) continue;
            const bits_t& mut_path = effects_[j];
            double p = anc_p;
            p *= w_gene_[j];
            p *= discount_if_subset(pathtype, mut_path);
            denoms_[s] += p;
            if (s < max_sites_) {
                mutate(bits_t(genotype).set(j), pathtype | mut_path, p);
            }
        }
    }

    double discount_if_subset(const bits_t& pathtype, const bits_t& mut_path) const {
        double p = 1.0;
        for (size_t i=0; i<th_path_.size(); ++i) {
            if (mut_path[i]) {
                if (pathtype[i]) {
                    p *= th_path_[i];
                } else {
                    return 1.0;
                }
            }
        }
        return p;
    }

    double discount(const std::valarray<size_t>& mut_route) const {
        double p = 1.0;
        bits_t pathtype;
        for (const auto j: mut_route) {
            const auto& mut_path = effects_[j];
            p *= discount_if_subset(pathtype, mut_path);
            pathtype |= mut_path;
        }
        return p;
    }

    bits_t translate(const size_t& mut_idx) const {
        bits_t mut_path;
        for (size_t j=0; j<annot_.size(); ++j) {
            mut_path.set(j, annot_[j][mut_idx]);
        }
        return mut_path;
    }
    const std::valarray<double>& w_gene_;
    const std::valarray<double>& th_path_;
    const std::vector<bits_t>& annot_;
    const size_t max_sites_;
    std::valarray<double> denoms_;
    std::vector<bits_t> effects_;
};

double GenotypeModel::calc_loglik(const std::valarray<double>& th_path) const {
    const size_t max_sites = nsam_with_s_.size() - 1;
    Denoms subcalc(w_gene_, th_path, annot_, max_sites);
    double loglik = 0.0;
    for (const auto& genotype: genot_) {
        loglik += subcalc.lnp_sample(genotype);
    }
    const auto lnD = subcalc.log();
    // std::cerr << "lnD: " << lnD << std::endl;
    // -inf, 0, D2, D3, ...
    for (size_t s=2; s<=max_sites; ++s) {
        loglik -= nsam_with_s_[s] * lnD[s];
    }
    return loglik;
}

std::string GenotypeModel::init_meta() {HERE;
    if (stage_ >= STEPS.size()) return "";
    auto oss = wtl::make_oss(2, std::ios::fixed);
    oss << "grid-" << STEPS.at(stage_) << ".tsv.gz";
    std::string outfile = oss.str();
    try {
        wtl::izfstream ist(outfile);
        std::cerr << "Reading: " << outfile << std::endl;
        read_results(ist);
        if (skip_ == 0) {
            ++stage_;
            outfile = init_meta();
        }
    } catch (std::ios::failure& e) {
        if (errno != 2) throw;
    }
    return outfile;
}

void GenotypeModel::read_results(std::istream& ist) {HERE;
    size_t max_count;
    double step;
    std::tie(max_count, std::ignore, step) = read_metadata(ist);
    stage_ = guess_stage(step);
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
}

void GenotypeModel::unit_test() {HERE;
    std::stringstream sst;
    sst <<
R"({
  "pathway": ["A", "B"],
  "annotation": ["0011", "1100"],
  "sample": ["0011", "0101", "1001", "0110", "1010", "1100"]
})";
    GenotypeModel model(sst, 4);
    model.run(false);
}

} // namespace likeligrid
