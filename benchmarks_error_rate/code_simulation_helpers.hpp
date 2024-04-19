//
// Created by alice on 08.09.21.
//

#ifndef LDPC4QKD_CODE_SIMULATION_HELPERS_HPP
#define LDPC4QKD_CODE_SIMULATION_HELPERS_HPP

#include <filesystem> // C++17

#include "external/json-6af826d/json.hpp"  // external json parser library

#include "read_ldpc_file_formats.hpp"
#include "src/rate_adaptive_code.hpp"


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
        return tmp / static_cast<double>(in.size());
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
    LDPC4QKD::RateAdaptiveCode<colptr_t, idx_t> load_ldpc_from_cscmat(
            const std::string &cscmat_file_path, const std::string &rate_adaption_file_path=""
    ) {
        auto pair = LDPC4QKD::read_matrix_from_cscmat<colptr_t, idx_t>(cscmat_file_path);
        auto colptr = pair.first;
        auto row_idx = pair.second;

        if(rate_adaption_file_path.empty()) {
            return LDPC4QKD::RateAdaptiveCode<colptr_t, idx_t>(colptr, row_idx);
        } else {
            std::vector<idx_t> rows_to_combine = read_rate_adaption_from_csv<idx_t>(rate_adaption_file_path);
            return LDPC4QKD::RateAdaptiveCode<colptr_t, idx_t>(colptr, row_idx, rows_to_combine);
        }
    }

    /*!
     * Loads LDPC code (and optionally also rate adaption) from files.
     * WARNING: if the templated types are too small, the numbers in the files are static_cast down!
     *
     * @tparam Bit type of matrix entires (only values zero and one are used, default to bool)
     * @tparam colptr_t unsigned integer type that fits ("number of non-zero matrix entries" + 1)
     * @tparam idx_t unsigned integer type fitting number of columns N (thus also number of rows M)
     * @param json_file_path path to json file, where the LDPC code is loaded from.
     *      (no QC-exponents allowed, just binary compressed sparse column (CSC) representation!)
     * @param rate_adaption_file_path Path to load rate adaption from.
     * (csv file of line index pairs, which are combined at each rate adaption step)
     *      If unspecified, no rate adaption will be available.
     * @return Rate adaptive code
     */
    template<typename Bit=bool,
            typename colptr_t=std::uint32_t,
            typename idx_t=std::uint32_t>
    LDPC4QKD::RateAdaptiveCode<colptr_t, idx_t> load_ldpc_from_json(
            const std::string &json_file_path, const std::string &rate_adaption_file_path = ""
    ) {
        try {
            std::ifstream fs(json_file_path);
            if (!fs) {
                throw std::runtime_error("Invalid file.");
            }

            // parse JSON file
            using json = nlohmann::json;
            std::ifstream f(json_file_path);
            json data = json::parse(f);

            if (data["format"] == "BINCSCJSON") { // binary matrix stored
                std::vector<colptr_t> colptr = data["colptr"];
                std::vector<idx_t> rowval = data["rowval"];

                if (rate_adaption_file_path.empty()) {
                    return LDPC4QKD::RateAdaptiveCode<colptr_t, idx_t>(colptr, rowval);
                } else {
                    std::vector<idx_t> rows_to_combine = read_rate_adaption_from_csv<idx_t>(rate_adaption_file_path);
                    return LDPC4QKD::RateAdaptiveCode<colptr_t, idx_t>(colptr, rowval, rows_to_combine);
                }
            } else if (data["format"] == "COMPRESSED_SPARSE_COLUMN") { // quasi-cyclic exponents stored
//                std::vector<std::uint64_t> colptr = data["colptr"];
//                std::vector<std::uint64_t> rowval = data["rowval"];
//                std::vector<std::uint64_t> nzval = data["nzval"];
                    // TODO reconstructing the matrix from qc-exponents and storing as sparse needs some work...
                throw std::runtime_error("Reading QC-exponents not supported yet! Use LDPCStorage.jl to convert to bincsc.json format!");
            } else {
                throw std::runtime_error("Unexpected format within json file.");
            }
        } catch (const std::exception &e) {
            std::stringstream s;
            s << "Failed to read LDPC code from file '" << json_file_path << "'. Reason:\n" << e.what() << "\n";
            throw std::runtime_error(s.str());
        } catch (...) {
            std::stringstream s;
            s << "Failed to read LDPC code from file '" << json_file_path << "' due to unknown error.";
            throw std::runtime_error(s.str());
        }
    }

    /*!
     * Loads LDPC code (and optionally also rate adaption) from files. Parser chosen depends on File ending!
     * WARNING: if the templated types are too small, the numbers in the files are static_cast down!
     *
     * @tparam Bit type of matrix entires (only values zero and one are used, default to bool)
     * @tparam colptr_t unsigned integer type that fits ("number of non-zero matrix entries" + 1)
     * @tparam idx_t unsigned integer type fitting number of columns N (thus also number of rows M)
     * @param file_path path to file, where the LDPC code is loaded from.
     *      (no QC-exponents allowed, just binary compressed sparse column (CSC) representation!)
     * @param rate_adaption_file_path Path to load rate adaption from.
     * (csv file of line index pairs, which are combined at each rate adaption step)
     *      If unspecified, no rate adaption will be available.
     * @return Rate adaptive code
     */
    template<typename Bit=bool,
            typename colptr_t=std::uint32_t,
            typename idx_t=std::uint32_t>
    LDPC4QKD::RateAdaptiveCode<colptr_t, idx_t> load_ldpc(
            const std::string &file_path, const std::string &rate_adaption_file_path=""
    ) {
        std::filesystem::path filePath = file_path;
        if (filePath.extension() == ".cscmat") {
            return load_ldpc_from_cscmat(file_path, rate_adaption_file_path);
        } else if (filePath.extension() == ".json") {
            return load_ldpc_from_json(file_path, rate_adaption_file_path);
        } else {
            throw std::runtime_error("Expected file with extension .cscmat or .json");
        }
    }

}

#endif //LDPC4QKD_CODE_SIMULATION_HELPERS_HPP
