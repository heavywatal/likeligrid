// -*- mode: c++; coding: utf-8 -*-
/*! @file gradient_descent.cpp
    @brief Implementation of GradientDescent class
*/
#include "gradient_descent.hpp"
#include "util.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/prandom.hpp>
#include <wtl/itertools.hpp>
#include <wtl/concurrent.hpp>
#include <wtl/scope.hpp>

namespace likeligrid {

GradientDescent::GradientDescent(const std::string& infile,
    const size_t max_sites,
    const unsigned int concurrency)
    : GradientDescent(wtl::izfstream(infile), max_sites, concurrency) {HERE;}

GradientDescent::GradientDescent(
    std::istream& ist,
    const size_t max_sites,
    const unsigned int concurrency)
    : model_(ist, max_sites),
      concurrency_(concurrency) {HERE;

}

void GradientDescent::run(std::ostream& ost) {HERE;
    auto at_exit = wtl::scope_exit([&ost,this](){
        std::cerr << std::endl;
        write(ost);
    });
    if (history_.empty()) {
        const size_t dimensions = model_.names().size();
        const std::valarray<double> initial_values(0.90, dimensions);
        const double loglik = model_.calc_loglik(initial_values);
        history_.emplace(initial_values, loglik);
    }
    for (auto it = max_iterator();
         it != history_.end();
         it = find_better(it)) {
    }
}

MapGrid::iterator GradientDescent::find_better(const MapGrid::iterator& prev_it) {
    auto task = [](decltype(model_) model, const std::valarray<double> theta) {
        // arguments are copied for each thread
        return std::make_pair(theta, model.calc_loglik(theta));
    };
    std::vector<std::future<std::pair<std::valarray<double>, double>>> futures;
    futures.reserve(concurrency_);

    auto surrounding = empty_neighbors_of(prev_it->first);
    std::shuffle(std::begin(surrounding), std::end(surrounding), wtl::sfmt());
    for (const auto& theta: surrounding) {
        futures.push_back(std::async(std::launch::async, task, model_, theta));
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

void GradientDescent::write(std::ostream& ost) {HERE;
    ost << "##max_count=" << 0U << "\n";
    ost << "##max_sites=" << model_.max_sites() << "\n";
    ost << "##step=" << 0.01 << "\n";
    ost << "loglik\t" << wtl::join(model_.names(), "\t") << "\n";
    for (const auto& p: history_) {
        ost << p.second << "\t"
            << wtl::str_join(p.first, "\t") << "\n";
    }
}

void GradientDescent::read_results(std::istream& ist) {HERE;
    size_t max_sites;
    std::tie(std::ignore, max_sites, std::ignore) = read_metadata(ist);

    std::string buffer;
    ist >> buffer; // loglik
    std::getline(ist, buffer); // header
    buffer.erase(0, 1); // \t
    const std::vector<std::string> colnames = wtl::split(buffer, "\t");

    bool is_resuming = false;
    if (model_.names() != colnames) {
        std::ostringstream oss;
        oss << "Contradiction in column names:\n"
            << "genotype file: " << model_.names() << "\n"
            << "result file:" << colnames;
        throw std::runtime_error(oss.str());
    } else {
        if (max_sites == model_.max_sites()) {
            is_resuming = true;
        }
    }

    while (std::getline(ist, buffer)) {
        std::istringstream iss(buffer);
        std::istream_iterator<double> it(iss);
        double loglik = *it;
        std::vector<double> vec(++it, std::istream_iterator<double>());
        if (is_resuming) {
            history_.emplace(std::valarray<double>(vec.data(), vec.size()), loglik);
        }
    }
}

void GradientDescent::unit_test() {HERE;
    std::stringstream sst;
    sst <<
R"({
  "pathway": ["A", "B"],
  "annotation": ["0011", "1100"],
  "sample": ["0011", "0101", "1001", "0110", "1010", "1100"]
})";
    GradientDescent searcher(sst, 4);
    searcher.run(std::cerr);
}

} // namespace likeligrid
