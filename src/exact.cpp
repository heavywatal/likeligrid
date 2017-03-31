// -*- mode: c++; coding: utf-8 -*-
/*! @file exact.cpp
    @brief Inplementation of ExactModel class
*/
#include "exact.hpp"

#include <functional>

#include <json.hpp>

#include <wtl/debug.hpp>
#include <wtl/exception.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/algorithm.hpp>
#include <wtl/math.hpp>
#include <wtl/os.hpp>
#include <wtl/eigen.hpp>

namespace likeligrid {

const std::vector<double> ExactModel::STEPS_ = {0.4, 0.2, 0.1, 0.05, 0.02, 0.01};
const std::vector<size_t> ExactModel::BREAKS_ = {5, 5, 5, 6, 5, 5};
bool ExactModel::SIGINT_RAISED_ = false;

ExactModel::ExactModel(const std::string& infile, const size_t max_sites) {HERE;
    nlohmann::json jso;
    {
        wtl::izfstream ist(infile);
        ist >> jso;
    }
    names_ = jso["pathway"].get<std::vector<std::string>>();
    const size_t npath = names_.size();
    annot_.reserve(npath);
    for (const std::string& s: jso["annotation"]) {
        annot_.emplace_back(s);
    }
    std::cerr << "annot_: " << annot_ << std::endl;

    const size_t nsam = jso["sample"].size();
    genot_.reserve(nsam);
    for (const std::string& s: jso["sample"]) {
        genot_.emplace_back(s);
    }
    // std::cerr << "genot_: " << genot_ << std::endl;

    const size_t ngene = genot_[0].size();
    nsam_with_s_.assign(ngene + 1, 0);  // at most
    std::valarray<double> s_gene(ngene);
    for (const auto& bits: genot_) {
        const size_t s = bits.count();
        ++nsam_with_s_[s];
        if (s > max_sites) continue;
        bits_t::size_type j = bits.find_first();
        while (j < bits_t::npos) {
            ++s_gene[j];
            j = bits.find_next(j);
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
    const auto final_max_s = nsam_with_s_.size() - 1;
    for (size_t s=2; s<=final_max_s; ++s) {
        lnp_const_ += nsam_with_s_[s] * std::log(wtl::factorial(s));
    }
    w_gene_ = s_gene / s_gene.sum();
    std::cerr << "s_gene : " << s_gene << std::endl;
    std::cerr << "w_gene_: " << w_gene_ << std::endl;
    for (size_t j=0; j<s_gene.size(); ++j) {
        if (s_gene[j] > 0) {
            lnp_const_ += s_gene[j] * std::log(w_gene_[j]);
        }
    }
    std::cerr << "lnp_const_: " << lnp_const_ << std::endl;

    // TODO
    a_pathway_.resize(npath);
    for (size_t j=0; j<npath; ++j) {
        for (size_t i=0; i<nsam; ++i) {
            size_t s = (genot_[i] & annot_[j]).count();
            if (s > 0) {
                a_pathway_[j] += --s;
            }
        }
    }
    std::cerr << "a_pathway_: " << a_pathway_ << std::endl;
}

void ExactModel::run(const std::string& infile) {HERE;
    const std::string outfile = init_meta(infile);
    if (outfile.empty()) return;
    const auto axes = make_axes();
    for (size_t j=0; j<names_.size(); ++j) {
        std::cerr << names_[j] << ": " << axes[j] << std::endl;
    }
    run_impl(outfile, wtl::itertools::product(axes));
    if (outfile != "/dev/stdout") {
        run(outfile);
    }
}

void ExactModel::run_impl(const std::string& outfile, wtl::itertools::Generator<std::valarray<double>>&& gen) const {HERE;
    auto buffer = wtl::make_oss();
    std::ios::openmode mode = std::ios::out;
    if (wtl::exists(outfile)) {
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
               << wtl::str_join(th_path, "\t") << "\n";
        if (gen.count() % 100 == 0) {  // snapshot for long run
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
        for (bits_t mut_gene(ngene, 1); mut_gene.any(); mut_gene <<= 1) {
            effects_.emplace_back(translate(mut_gene));
        }
        // std::cerr << "effects_: " << effects_ << std::endl;
        mutate(bits_t(ngene, 0), bits_t(annot.size(), 0), 1.0);
        // std::cerr << "denoms_: " << denoms_ << std::endl;
    }
    const std::valarray<double> log() const {return std::log(denoms_);}

  private:
    void mutate(const bits_t& genotype, const bits_t& pathtype, const double anc_p) {
        const size_t s = genotype.count() + 1;
        for (bits_t mut_gene(genotype.size(), 1); mut_gene.any(); mut_gene <<= 1) {
            if ((genotype & mut_gene).any()) continue;
            const size_t pos = mut_gene.find_first();
            const bits_t& mut_path = effects_[pos];
            double p = anc_p;
            p *= w_gene_[pos];
            p *= discount(pathtype, mut_path);
            // std::cout << (genotype | mut_gene) << " " << (pathtype | mut_path) << " " << p << std::endl;
            denoms_[s] += p;
            if (s < max_sites_) {
                mutate(genotype | mut_gene, pathtype | mut_path, p);
            }
        }
    }

    double discount(const bits_t& pathtype, const bits_t& mut_path) const {
        if (mut_path.is_subset_of(pathtype)) {
            double p = 1.0;
            const bits_t recurrent = mut_path & pathtype;
            bits_t::size_type j = recurrent.find_first();
            while (j < bits_t::npos) {
                p *= th_path_[j];
                j = recurrent.find_next(j);
            }
            return p;
        } else {return 1.0;}
    }

    bits_t translate(const bits_t& mut_gene) const {
        bits_t mut_path(annot_.size(), 0);
        for (size_t j=0; j<annot_.size(); ++j) {
            mut_path.set(j, (annot_[j] & mut_gene).any());
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

double ExactModel::calc_loglik(const std::valarray<double>& th_path) const {
    const size_t max_sites = nsam_with_s_.size() - 1;
    double loglik = (a_pathway_ * std::log(th_path)).sum();
    const auto lnD = Denoms(w_gene_, th_path, annot_, max_sites).log();
    // std::cout << "lnD: " << lnD << std::endl;
    // -inf, 0, D2, D3, ...
    for (size_t s=2; s<=max_sites; ++s) {
        loglik -= nsam_with_s_[s] * lnD[s];
    }
    return loglik += lnp_const_;
}

inline std::vector<std::valarray<double>>
make_vicinity(const std::valarray<double>& center, const size_t breaks, const double radius, const double max=2.0) {
    std::vector<std::valarray<double>> axes;
    axes.reserve(center.size());
    for (const double x: center) {
        Eigen::ArrayXd axis = Eigen::ArrayXd::LinSpaced(breaks, x + radius, x - radius);
        axis = (axis * 100.0).round() / 100.0;  // grid precision = 0.01
        axes.push_back(wtl::eigen::valarray(wtl::eigen::filter(axis, (0.0 < axis) * (axis < max))));
    }
    return axes;
}

std::vector<std::valarray<double>> ExactModel::make_axes() const {HERE;
    if (mle_params_.size() > 0) {
        std::cerr << "mle_params_: " << mle_params_ << std::endl;
        return make_vicinity(mle_params_, BREAKS_.at(stage_), STEPS_.at(stage_), 2.0);
    } else {
        const double step = STEPS_[0];
        const size_t breaks = BREAKS_[0];
        Eigen::ArrayXd axis = Eigen::VectorXd::LinSpaced(breaks, 2.0, step).array();
        return std::vector<std::valarray<double>>(names_.size(), wtl::eigen::valarray(axis));
    }
}

std::string ExactModel::init_meta(const std::string& infile) {HERE;
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

bool ExactModel::read_results(const std::string& infile) {HERE;
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

size_t ExactModel::read_metadata(std::istream& ist) {HERE;
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
    auto pred = std::bind(wtl::approx<double>, std::placeholders::_1, step);
    const auto it = std::find_if(STEPS_.begin(), STEPS_.end(), pred);
    if (it == STEPS_.end()) throw std::runtime_error("invalid step size");
    stage_ = it - STEPS_.begin();
    return max_count;
}

size_t ExactModel::read_body(std::istream& ist) {HERE;
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
            std::vector<double> v(++it, std::istream_iterator<double>());
            mle_params_ = std::valarray<double>(v.data(), v.size());
        }
    }
    return i;
}

void ExactModel::unit_test() {HERE;
    ExactModel model("test.json.gz", 2);
    model.run("/dev/null");
}

} // namespace likeligrid
