// -*- mode: c++; coding: utf-8 -*-
/*! @file gradient_descent.cpp
    @brief Implementation of GradientDescent class
*/
#include "gradient_descent.hpp"
#include "genotype.hpp"
#include "util.hpp"

#include <sfmt.hpp>
#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/itertools.hpp>
#include <wtl/concurrent.hpp>
#include <wtl/scope.hpp>

namespace likeligrid {

// std::unique_ptr needs to know GenotypeModel implementation
GradientDescent::~GradientDescent() = default;

GradientDescent::GradientDescent(const std::string& infile,
    const size_t max_sites,
    const std::pair<size_t, size_t>& epistasis_pair,
    const unsigned int concurrency)
    : concurrency_(concurrency) {HERE;
    std::string genotype_file;
    size_t prev_max_sites;
    std::string prev_epistasis;
    std::tie(genotype_file, prev_max_sites, prev_epistasis)
        = read_results(infile);
    std::cerr << "genotype: " << genotype_file << std::endl;
    model_ = std::make_unique<GenotypeModel>(genotype_file, max_sites);
    size_t dimensions = model_->names().size();
    bool is_restarting = (max_sites != prev_max_sites);
    if (epistasis_pair.first != epistasis_pair.second) {
        if (prev_epistasis.empty()) {
            is_restarting = true;
        } else {
            const std::string epistasis_name =
                model_->names().at(epistasis_pair.first) + ":" +
                model_->names().at(epistasis_pair.second);
            if (epistasis_name != prev_epistasis) {
                throw std::runtime_error("epistasis_name != prev_epistasis");
            }
        }
        model_->set_epistasis(epistasis_pair);
        ++dimensions;
    }
    const auto it = const_max_iterator();
    std::cerr << "old best: " << *it << std::endl;
    if (is_restarting) {
        const auto& old_best = it->first;
        std::valarray<double> new_start(1.0, dimensions);
        std::copy(std::begin(old_best), std::end(old_best), std::begin(new_start));
        history_.clear();
        const double loglik = model_->calc_loglik(new_start);
        history_.emplace(new_start, loglik);
        std::cerr << "new start: " << *history_.begin() << std::endl;
    }
}

GradientDescent::GradientDescent(
    std::istream& ist,
    const size_t max_sites,
    const std::pair<size_t, size_t>& epistasis_pair,
    const unsigned int concurrency)
    : model_(std::make_unique<GenotypeModel>(ist, max_sites)),
      concurrency_(concurrency) {HERE;
    size_t dimensions = model_->names().size();
    if (epistasis_pair.first != epistasis_pair.second) {
        model_->set_epistasis(epistasis_pair);
        ++dimensions;
    }
    const std::valarray<double> new_start(1.0, dimensions);
    const double loglik = model_->calc_loglik(new_start);
    history_.emplace(new_start, loglik);
}

void GradientDescent::run(std::ostream& ost) {HERE;
    auto at_exit = wtl::scope_exit([&ost,this](){
        std::cerr << "\n" << *const_max_iterator() << std::endl;
        write(ost);
    });
    for (auto it = max_iterator();
         it != history_.end();
         it = find_better(it)) {
    }
}

MapGrid::iterator GradientDescent::find_better(const MapGrid::iterator& prev_it) {
    auto task = [](GenotypeModel model, const std::valarray<double> theta) {
        // arguments are copied for each thread
        return std::make_pair(theta, model.calc_loglik(theta));
    };
    std::vector<std::future<std::pair<std::valarray<double>, double>>> futures;
    futures.reserve(concurrency_);

    auto surrounding = empty_neighbors_of(prev_it->first);
    std::shuffle(std::begin(surrounding), std::end(surrounding), wtl::sfmt());
    for (const auto& theta: surrounding) {
        futures.push_back(std::async(std::launch::async, task, *model_, theta));
        if (futures.size() == concurrency_) {
            auto better_it = prev_it;
            for (auto& ftr: futures) {
                auto result_it = history_.insert(ftr.get()).first;
                std::cerr << "." << std::flush;
                if (result_it->second > better_it->second) {
                    better_it = result_it;
                }
            }
            if (better_it != prev_it) {
                return better_it;
            }
            futures.clear();
        }
    }
    return history_.end();
}

std::vector<std::valarray<double>> GradientDescent::empty_neighbors_of(const std::valarray<double>& center) {
    const auto axes = make_vicinity(center, 3, 0.01);
    auto iter = wtl::itertools::product(axes);
    std::vector<std::valarray<double>> empty_neighbors;
    empty_neighbors.reserve(iter.max_count());
    for (const auto& x: iter()) {
        if (history_.find(x) == history_.end()) {
            empty_neighbors.push_back(x);
        }
    }
    return empty_neighbors;
}

std::string GradientDescent::outfile() const {HERE;
    std::ostringstream oss;
    oss << "gradient-s" << model_->max_sites();
    if (model_->epistasis_pair().first != model_->epistasis_pair().second) {
        oss << "-e" << model_->epistasis_pair().first
            <<  "x" << model_->epistasis_pair().second;
    }
    oss << ".tsv.gz";
    return oss.str();
}

void GradientDescent::write(std::ostream& ost) {HERE;
    ost << "##genotype_file=" << model_->filename() << "\n";
    ost << "##max_sites=" << model_->max_sites() << "\n";
    ost << "##max_count=" << 0U << "\n";
    ost << "##step=" << 0.01 << "\n";
    ost << "loglik\t" << wtl::join(model_->names(), "\t");
    const auto& epistasis_pair = model_->epistasis_pair();
    if (epistasis_pair.first != epistasis_pair.second) {
        ost << "\t" << model_->names().at(epistasis_pair.first)
            <<  ":" << model_->names().at(epistasis_pair.second);
    }
    ost << "\n";
    for (const auto& p: history_) {
        ost << p.second << "\t"
            << wtl::str_join(p.first, "\t") << "\n";
    }
}

std::tuple<std::string, size_t, std::string> GradientDescent::read_results(const std::string& infile) {HERE;
    wtl::izfstream ist(infile);
    std::string genotype_file;
    size_t prev_max_sites;
    std::tie(genotype_file, prev_max_sites, std::ignore, std::ignore) = read_metadata(ist);

    std::string buffer;
    ist >> buffer; // loglik
    std::getline(ist, buffer); // header
    buffer.erase(0, 1); // \t
    const std::vector<std::string> colnames = wtl::split(buffer, "\t");
    std::string epistasis_name = colnames.back();
    if (epistasis_name.find(':') == std::string::npos) epistasis_name.clear();

    while (std::getline(ist, buffer)) {
        std::istringstream iss(buffer);
        std::istream_iterator<double> it(iss);
        double loglik = *it;
        std::vector<double> vec(++it, std::istream_iterator<double>());
        history_.emplace(std::valarray<double>(vec.data(), vec.size()), loglik);
    }
    return {genotype_file, prev_max_sites, epistasis_name};
}

MapGrid::iterator GradientDescent::max_iterator() {HERE;
    return std::max_element(std::begin(history_), std::end(history_),
        [](MapGrid::value_type& x, MapGrid::value_type& y){
            return x.second < y.second;
        });
}

MapGrid::const_iterator GradientDescent::const_max_iterator() const {HERE;
    return std::max_element(std::begin(history_), std::end(history_),
        [](const MapGrid::value_type& x, const MapGrid::value_type& y){
            return x.second < y.second;
        });
}

void GradientDescent::test() {HERE;
    std::stringstream sst;
    sst <<
R"({
  "pathway": ["A", "B"],
  "annotation": ["0011", "1100"],
  "sample": ["0011", "0101", "1001", "0110", "1010", "1100"]
})";
    GradientDescent searcher(sst, 4, {0, 1});
    searcher.run(std::cout);
}

} // namespace likeligrid
