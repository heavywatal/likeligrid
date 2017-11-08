/*! @file genotype.hpp
    @brief Interface of GenotypeModel class
*/
#pragma once
#ifndef LIKELIGRID_GENOTYPE_HPP_
#define LIKELIGRID_GENOTYPE_HPP_

#include "typedef.hpp"
#include "util.hpp"

#include <wtl/exception.hpp>

#include <string>
#include <vector>
#include <valarray>

namespace likeligrid {

class GenotypeModel {
  public:
    GenotypeModel(std::istream& ist, const size_t max_sites) {
        init(ist, max_sites);
    }
    GenotypeModel(std::istream&& ist, const size_t max_sites)
    : GenotypeModel(ist, max_sites) {}
    GenotypeModel(const std::string&, const size_t max_sites);

    void set_epistasis(const std::pair<size_t, size_t>& pair);

    double calc_loglik(const std::valarray<double>& theta);
    void benchmark(const size_t);

    // getter
    const std::string& filename() const {return filename_;}
    const std::vector<std::string>& names() const {return names_;}
    const std::pair<size_t, size_t>& epistasis_pair() const {return epistasis_pair_;}
    size_t max_sites() const {return max_sites_;}

    static void test();

  private:
    void init(std::istream&, const size_t max_sites);

    double lnp_sample(const bits_t& genotype) const;

    void mutate(const bits_t& genotype, const bits_t& pathtype, const double anc_p) {
        const auto s = genotype.count() + 1u;
        for (size_t j=0u; j<num_genes_; ++j) {
            if (genotype[j]) continue;
            const bits_t& mut_path = effects_[j];
            double p = anc_p;
            p *= w_gene_[j];
            p *= discount_if_subset(pathtype, mut_path);
            p *= epistasis(pathtype, mut_path);
            denoms_[s] += p;
            if (s < max_sites_) {
                if (wtl::SIGINT_RAISED()) {throw wtl::KeyboardInterrupt();}
                mutate(bits_t(genotype).set(j), pathtype | mut_path, p);
            }
        }
    }

    double discount_if_subset(const bits_t& pathtype, const bits_t& mut_path) const {
        double p = 1.0;
        for (size_t i=0u; i<num_pathways_; ++i) {
            if (mut_path[i]) {
                if (pathtype[i]) {
                    p *= theta_[i];
                } else {
                    return 1.0;
                }
            }
        }
        return p;
    }

    double epistasis(const bits_t& pathtype, const bits_t& mut_path) const {
        if (!with_epistasis_) return 1.0;
        if (pathtype[epistasis_pair_.first]) {
            if (pathtype[epistasis_pair_.second]) return 1.0;
            if (mut_path[epistasis_pair_.second]) return theta_[num_pathways_];
        }
        if (pathtype[epistasis_pair_.second]) {
            if (mut_path[epistasis_pair_.first]) return theta_[num_pathways_];
        }
        if (mut_path[epistasis_pair_.first]) {
            if (mut_path[epistasis_pair_.second]) return theta_[num_pathways_];
        }
        return 1.0;
    }

    double discount(const std::valarray<size_t>& mut_route) const {
        double p = 1.0;
        bits_t pathtype;
        for (const auto j: mut_route) {
            const auto& mut_path = effects_[j];
            p *= discount_if_subset(pathtype, mut_path);
            p *= epistasis(pathtype, mut_path);
            pathtype |= mut_path;
        }
        return p;
    }

    bits_t translate(const size_t& mut_idx) const {
        bits_t mut_path;
        for (size_t j=0u; j<num_pathways_; ++j) {
            mut_path.set(j, annot_[j][mut_idx]);
        }
        return mut_path;
    }

    // initialized in constructor
    std::string filename_ = "-";
    std::vector<std::string> names_;
    size_t num_pathways_;
    std::vector<bits_t> annot_;
    std::vector<bits_t> genot_;
    std::valarray<double> w_gene_;
    size_t num_genes_;
    std::vector<size_t> nsam_with_s_;
    size_t max_sites_;
    std::vector<bits_t> effects_;

    // updated in calc_loglik()
    std::valarray<double> theta_;
    std::valarray<double> denoms_;
    std::pair<size_t, size_t> epistasis_pair_;
    bool with_epistasis_ = false;
};

} // namespace likeligrid

#endif // LIKELIGRID_GENOTYPE_HPP_
