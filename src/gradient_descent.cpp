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

#include <boost/filesystem.hpp>

namespace likeligrid {

namespace fs = boost::filesystem;

// std::unique_ptr needs to know GenotypeModel implementation
GradientDescent::~GradientDescent() = default;

GradientDescent::GradientDescent(
    std::istream& ist,
    const size_t max_sites,
    const std::pair<size_t, size_t>& epistasis_pair,
    const bool pleiotropy,
    const unsigned int concurrency)
    : model_(std::make_unique<GenotypeModel>(ist, max_sites)),
      concurrency_(concurrency)
    {HERE;
    model_->set_epistasis(epistasis_pair, pleiotropy);
}

GradientDescent::GradientDescent(
    const std::string& infile,
    const size_t max_sites,
    const std::pair<size_t, size_t>& epistasis_pair,
    const bool pleiotropy,
    const unsigned int concurrency)
    : concurrency_(concurrency)
    {HERE;

    std::string genotype_file = infile;
    if (wtl::endswith(infile, ".tsv.gz")) {// previous result
        std::string previous = fs::path(infile).filename().string();
        wtl::izfstream ist(infile);
        size_t prev_max_sites;
        std::tie(genotype_file, prev_max_sites, std::ignore, std::ignore) = read_metadata(ist);
        std::tie(std::ignore, std::ignore, starting_point_) = read_body(ist);
        std::ostringstream oss;
        oss << "grad-from-s" << prev_max_sites << "-" << previous;
        outfile_ = oss.str();
    } else {
        outfile_ = "grad-from-center.tsv.gz";
    }
    model_ = std::make_unique<GenotypeModel>(genotype_file, max_sites);
    model_->set_epistasis(epistasis_pair, pleiotropy);
}

void GradientDescent::run(std::ostream& ost) {HERE;
    auto at_exit = wtl::scope_exit([&ost,this](){
        std::cerr << "\n" << *const_max_iterator() << std::endl;
        write(ost);
    });

    std::valarray<double> new_start(1.0, model_->names().size());
    std::copy(std::begin(starting_point_), std::end(starting_point_), std::begin(new_start));
    history_.emplace(new_start, model_->calc_loglik(new_start));
    std::cerr << "start: " << *history_.begin() << std::endl;

    for (auto it = max_iterator();
         it != history_.end();
         it = find_better(it)) {
    }
}

MapGrid::iterator GradientDescent::find_better(const MapGrid::iterator& prev_it) {
    auto better_it = prev_it;
    auto task = [](GenotypeModel model, const std::valarray<double> theta) {
        // arguments are copied for each thread
        return std::make_pair(theta, model.calc_loglik(theta));
    };
    std::vector<std::future<std::pair<std::valarray<double>, double>>> futures;
    futures.reserve(concurrency_);
    wtl::ThreadPool pool(concurrency_);
    const auto candidates = empty_neighbors_of(prev_it->first);
    for (const auto& theta: candidates) {
        futures.push_back(pool.submit(task, *model_, theta));
        if ((futures.size() == concurrency_) || (&theta == &candidates.back())) {
            for (auto& ftr: futures) {
                auto result_it = history_.insert(ftr.get()).first;
                std::cerr << "." << std::flush;
                if (result_it->second > better_it->second) {
                    better_it = result_it;
                } else if (wtl::approx(result_it->second, better_it->second)) {
                    if (d2_from_neutral(result_it->first) < d2_from_neutral(better_it->first)) {
                        better_it = result_it;
                    }
                }
            }
            if (better_it != prev_it) {
                std::cerr << "*" << std::flush;
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
    std::shuffle(std::begin(empty_neighbors), std::end(empty_neighbors), wtl::sfmt64());
    return empty_neighbors;
}

void GradientDescent::write(std::ostream& ost) {HERE;
    ost << "##genotype_file=" << model_->filename() << "\n";
    ost << "##max_sites=" << model_->max_sites() << "\n";
    ost << "##max_count=" << 0u << "\n";
    ost << "##step=" << 0.01 << "\n";
    ost << "loglik\t";
    wtl::join(model_->names(), ost, "\t") << "\n";
    for (const auto& p: history_) {
        ost << p.second << "\t";
        wtl::join(p.first, ost, "\t") << "\n";
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
    GradientDescent searcher(sst, 4, {0, 1}, false);
    searcher.run(std::cout);
    std::cout << *searcher.max_iterator() << std::endl;
}

} // namespace likeligrid
