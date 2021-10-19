//
// Created by alice on 08.09.21.
//

#ifndef LDPC4QKD_CODE_SIMULATION_HELPERS_HPP
#define LDPC4QKD_CODE_SIMULATION_HELPERS_HPP

#include "read_scsmat_format.hpp"

namespace LDPC4QKD::CodeSimulationHelpers {
    template<typename T>
    void noise_bitstring_inplace(std::mt19937_64 &rng, std::vector<T> &src, double err_prob) {
        std::bernoulli_distribution distribution(err_prob);

        for (std::size_t i = 0; i < src.size(); i++) {
            if (distribution(rng)) {
                src[i] = !src[i];
            } else {
                src[i] = src[i];
            }
        }
    }


    // Shannon binary entropy
    template<typename T>
    double h2(T p) {
        return -p * ::log2(p) - (1 - p) * ::log2(1 - p);
    }


    template<typename T>
    double avg(const std::vector<T> &in) {
        double tmp{};
        for (auto i : in) {
            tmp += static_cast<double>(i);
        }
        return tmp / in.size();
    }


    /*!
     * Loads LDPC code (and optionally also rate adaption) from files.
     * WARNING: if the templated types are too small, the numbers in the files are static_cast down!
     *
     * @tparam Bit type of matrix entires (only values zero and one are used, default to bool)
     * @tparam colptr_t unsigned integer type that fits ("number of non-zero matrix entries" + 1)
     * @tparam idx_t unsigned integer type fitting number of columns N (thus also number of rows M)
     * @param cscmat_file_path path to CSCMAT file, where the LDPC code is loaded from.
     *      (no QC-exponents allowed, just binary compressed sparse column (CSC) representation!)
     * @param rate_adaption_file_path Path to load rate adaption from.
     * (csv file of line index pairs, which are combined at each rate adaption step)
     *      If unspecified, no rate adaption will be available.
     * @return Rate adaptive code
     */
    template<typename Bit=bool,
            typename colptr_t=std::uint32_t,
            typename idx_t=std::uint32_t>
    LDPC4QKD::RateAdaptiveCode<Bit, colptr_t, idx_t> load_ldpc(
            const std::string &cscmat_file_path, const std::string &rate_adaption_file_path=""
    ) {
        auto pair = LDPC4QKD::read_matrix_from_cscmat<colptr_t, idx_t>(cscmat_file_path);
        auto colptr = pair.first;
        auto row_idx = pair.second;

        if(rate_adaption_file_path.empty()) {
            return LDPC4QKD::RateAdaptiveCode<Bit, colptr_t, idx_t>(colptr, row_idx);
        } else {
            std::vector<idx_t> rows_to_combine = read_rate_adaption_from_csv<idx_t>(rate_adaption_file_path);
            return LDPC4QKD::RateAdaptiveCode<Bit, colptr_t, idx_t>(colptr, row_idx, rows_to_combine);
        }
    }
}

#endif //LDPC4QKD_CODE_SIMULATION_HELPERS_HPP
