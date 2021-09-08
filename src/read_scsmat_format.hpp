//
// Created by alice on 02.08.21.
//

#ifndef LDPC4QKD_READ_SCSMAT_FORMAT_HPP
#define LDPC4QKD_READ_SCSMAT_FORMAT_HPP

#include <fstream>
#include <sstream>

namespace LDPC4QKD {

    namespace HelpersReadSCSMAT {

        /// Helper method to split a single-line string consisting of space-separated integers into a vector of integers.
        /// CHANGES INPUT STRING!!
        template<typename idx_t=std::size_t>
        std::vector <idx_t> helper_parse_space_sep_ints(std::string &in) {
            std::vector <idx_t> ret;

            std::size_t start = 0;
            std::size_t end = 0;
            constexpr auto trimmed_whitespace = " \t\n\r\f\v";
            constexpr auto delimiter = ' ';

            // trim leading and trailing whitespace
            in.erase(in.find_last_not_of(trimmed_whitespace) + 1);
            in.erase(0, in.find_first_not_of(trimmed_whitespace));

            while (std::string::npos != (end = in.find(delimiter, start))) {
                ret.push_back(static_cast<idx_t>(std::stoll(in.substr(start, end - start))));
                start = end + 1;
                while (::isspace(in[start]))
                    ++start;
            }
            ret.push_back(static_cast<idx_t>(std::stoll(in.substr(start))));
            return ret;
        }

    }

    /// reads two arrays of integers from a .cscmat file.
    /// These are called `colptr` and `rowval`.
    /// They specify a binary LDPC matrx stored in compressed sparse column (CSC) format
    template<typename colptr_t=std::uint32_t, // integer type that fits ("number of non-zero matrix entries" + 1)
            typename idx_t=std::uint16_t>
    std::pair <std::vector<colptr_t>, std::vector<idx_t>> read_matrix_from_cscmat(const std::string &file_path) {
        try {
            std::ifstream fs(file_path);
            if (!fs) {
                throw std::runtime_error("Stream object invalid.");
            }

            std::string current_line;

            while (getline(fs, current_line) && current_line[0] == '#') {
                // Do nothing. Ignores comments at the beginning of the file.
            }

            getline(fs, current_line); // ignores another line containing meta-information
            getline(fs, current_line); // ignores empty line
            std::vector<colptr_t> colptr = HelpersReadSCSMAT::helper_parse_space_sep_ints<colptr_t>(current_line);

            getline(fs, current_line); // ignores empty line
            getline(fs, current_line);
            std::vector<idx_t> rowval = HelpersReadSCSMAT::helper_parse_space_sep_ints<idx_t>(current_line);

            return std::make_pair(colptr, rowval);
        }
        catch (const std::exception &e) {
            std::stringstream s;
            s << "Failed to read LDPC code from file '" << file_path << "'. Reason:\n" << e.what() << "\n";
            throw std::runtime_error(s.str());
        }
        catch (...) {
            std::stringstream s;
            s << "Failed to read LDPC code from file '" << file_path << "' due to unknown error.";
            throw std::runtime_error(s.str());
        }
    }
}

#endif //LDPC4QKD_READ_SCSMAT_FORMAT_HPP
