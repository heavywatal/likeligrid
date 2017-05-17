// -*- mode: c++; coding: utf-8 -*-
/*! @file genotype.hpp
    @brief Interface of GenotypeModel class
*/
#pragma once
#ifndef LIKELIGRID_GENOTYPE_HPP_
#define LIKELIGRID_GENOTYPE_HPP_

#include "typedef.hpp"

#include <string>
#include <vector>
#include <valarray>

namespace likeligrid {

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

class GenotypeModel {
  public:
    GenotypeModel(std::istream&, const size_t max_sites);

    const std::vector<std::string>& names() const {return names_;}
    size_t max_sites() const {return max_sites_;}
    size_t num_genes() const {return w_gene_.size();}

    double calc_loglik(const std::valarray<double>& th_path) {
        denoms_ = 0.0;
        th_path_ = th_path;
        mutate(bits_t(), bits_t(), 1.0);
        // std::cerr << "denoms_: " << denoms_ << std::endl;
        double loglik = 0.0;
        for (const auto& genotype: genot_) {
            loglik += lnp_sample(genotype);
        }
        const auto lnD = std::log(denoms_);
        // std::cerr << "lnD: " << lnD << std::endl;
        // -inf, 0, D2, D3, ...
        for (size_t s=2; s<=max_sites_; ++s) {
            loglik -= nsam_with_s_[s] * lnD[s];
        }
        return loglik;
    }
    static void unit_test();

  private:
    double lnp_sample(const bits_t& genotype) const {
        double p = 0.0;
        const double p_basic = slice_prod(w_gene_, genotype);
        auto mut_route = to_indices(genotype);
        do {
            p += p_basic * discount(mut_route);
        } while (std::next_permutation(std::begin(mut_route), std::end(mut_route)));
        return std::log(p);
    }

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
    std::vector<std::string> names_;
    std::vector<bits_t> annot_;
    std::vector<bits_t> genot_;
    std::valarray<double> w_gene_;
    std::vector<size_t> nsam_with_s_;
    size_t max_sites_;
    std::valarray<double> th_path_;
    std::valarray<double> denoms_;
    std::vector<bits_t> effects_;
};

} // namespace likeligrid

#endif // LIKELIGRID_GENOTYPE_HPP_
