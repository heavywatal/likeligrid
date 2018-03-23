/*! @file genotype.hpp
    @brief Interface of GenotypeModel class
*/
#pragma once
#ifndef LIKELIGRID_GENOTYPE_HPP_
#define LIKELIGRID_GENOTYPE_HPP_

#include <string>
#include <vector>
#include <valarray>
#include <bitset>

namespace likeligrid {

using bits_t = std::bitset<128>;

class GenotypeModel {
  public:
    GenotypeModel(std::istream& ist, size_t max_sites) {
        init(ist, max_sites);
    }
    GenotypeModel(std::istream&& ist, size_t max_sites)
    : GenotypeModel(ist, max_sites) {}
    GenotypeModel(const std::string&, size_t max_sites);

    bool set_epistasis(const std::pair<size_t, size_t>& pair, bool pleiotropy=false);

    double calc_loglik(const std::valarray<double>& theta);
    void benchmark(size_t);

    // getter
    const std::string& filename() const {return filename_;}
    const std::vector<std::string>& names() const {return names_;}
    const std::pair<size_t, size_t>& epistasis_pair() const {return epistasis_pair_;}
    size_t max_sites() const {return max_sites_;}

  private:
    void init(std::istream&, size_t max_sites);

    double lnp_sample(const bits_t& genotype) const;

    void mutate(const bits_t& genotype=bits_t(), const bits_t& pathtype=bits_t(),
                double anc_lnp=0.0, double open_lnp=0.0);

    double ln_theta_if_subset(const bits_t& pathtype, const bits_t& mut_path) const {
        double lnp = 0.0;
        for (size_t i=0u; i<num_pathways_; ++i) {
            if (mut_path[i]) {
                if (pathtype[i]) {
                    lnp += ln_theta_[i];
                } else {
                    return 0.0;
                }
            }
        }
        return lnp;
    }

    double ln_theta_if_paired(const bits_t& pathtype, const bits_t& mut_path) const {
        if (pathtype[epistasis_pair_.first]) {
            if (pathtype[epistasis_pair_.second]) return 0.0;
            if (mut_path[epistasis_pair_.second]) return ln_theta_[epistasis_idx_];
        }
        if (pathtype[epistasis_pair_.second]) {
            if (mut_path[epistasis_pair_.first]) return ln_theta_[epistasis_idx_];
        }
        if (mut_path[epistasis_pair_.first]) {
            if (mut_path[epistasis_pair_.second]) return ln_theta_[pleiotropy_idx_];
        }
        return 0.0;
    }

    double sum_ln_theta(const std::valarray<size_t>& mut_route) const {
        double lnp = 0.0;
        bits_t pathtype;
        for (const auto j: mut_route) {
            const auto& mut_path = effects_[j];
            lnp += ln_theta_if_subset(pathtype, mut_path);
            if (epistasis_) {lnp += ln_theta_if_paired(pathtype, mut_path);}
            pathtype |= mut_path;
        }
        return lnp;
    }

    bits_t translate(size_t mut_idx) const {
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
    std::valarray<double> ln_w_gene_;
    size_t num_genes_;
    std::vector<size_t> nsam_with_s_;
    size_t max_sites_;
    std::vector<bits_t> effects_;

    // updated in calc_loglik()
    std::valarray<double> ln_theta_;
    std::valarray<double> ln_denoms_;
    std::pair<size_t, size_t> epistasis_pair_;
    bool epistasis_ = false;
    size_t epistasis_idx_ = -1u;
    size_t pleiotropy_idx_ = epistasis_idx_;
};

} // namespace likeligrid

#endif // LIKELIGRID_GENOTYPE_HPP_
