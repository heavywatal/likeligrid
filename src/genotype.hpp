// -*- mode: c++; coding: utf-8 -*-
/*! @file genotype.hpp
    @brief Interface of GenotypeModel class
*/
#pragma once
#ifndef LIKELIGRID_GENOTYPE_HPP_
#define LIKELIGRID_GENOTYPE_HPP_

#include "typedef.hpp"
#include "util.hpp"

#include <string>
#include <vector>
#include <valarray>

#include <wtl/exception.hpp>

namespace likeligrid {

class GenotypeModel {
  public:
    GenotypeModel(std::istream&, const size_t max_sites);

    double calc_loglik(const std::valarray<double>& th_path);
    void benchmark(const size_t);

    // getter
    const std::vector<std::string>& names() const {return names_;}
    size_t max_sites() const {return max_sites_;}

    static void unit_test();

  private:
    double lnp_sample(const bits_t& genotype) const;

    void mutate(const bits_t& genotype, const bits_t& pathtype, const double anc_p) {
        const auto s = genotype.count() + 1;
        for (size_t j=0; j<num_genes_; ++j) {
            if (genotype[j]) continue;
            const bits_t& mut_path = effects_[j];
            double p = anc_p;
            p *= w_gene_[j];
            p *= discount_if_subset(pathtype, mut_path);
            denoms_[s] += p;
            if (s < max_sites_) {
                if (SIGINT_RAISED) {throw wtl::KeyboardInterrupt();}
                mutate(bits_t(genotype).set(j), pathtype | mut_path, p);
            }
        }
    }

    double discount_if_subset(const bits_t& pathtype, const bits_t& mut_path) const {
        double p = 1.0;
        for (size_t i=0; i<num_pathways_; ++i) {
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
        for (size_t j=0; j<num_pathways_; ++j) {
            mut_path.set(j, annot_[j][mut_idx]);
        }
        return mut_path;
    }

    // initialized in constructor
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
    std::valarray<double> th_path_;
    std::valarray<double> denoms_;
};

} // namespace likeligrid

#endif // LIKELIGRID_GENOTYPE_HPP_
