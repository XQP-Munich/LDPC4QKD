//
// Created by alice on 23.04.21.
//

#ifndef LDPC4QKD_HELPERS_FOR_TESTING_HPP
#define LDPC4QKD_HELPERS_FOR_TESTING_HPP

#include <fstream>
#include <iostream>
#include <random>
#include <array>

namespace HelpersForTests {

    using Bit = bool;

/*!
 * This is only used for the tests to verify agreement between vectors loaded by
 * this code and Luhn/Freiwang Python code.
 * Note: due to the bitsize conversions this hash function has no guarantees about any properties.
 * Adapted from https://stackoverflow.com/a/27216842
 *
 * Corresponding Python code:
    ```python3
    # numpy v1.19.0
    def hash_vector(vec):
        assert len(vec.shape) == 1, "only accepts 1d vectors"
        seed = np.uint32(vec.shape[0])
        for i in vec:
            seed ^= np.uint32(i) + np.uint32(0x9e3779b9) + (seed << np.uint32(6)) + (seed >> np.uint32(2))
        return seed
    ```
* Corresponding Julia code:
    ```julia
    # Julia v1.6
    function hash_vector(vec::AbstractArray{T} where T <: Integer)
        seed = UInt32(length(vec))
        for i in vec
            seed = xor(seed, UInt32(i) + UInt32(0x9e3779b9) + (seed << UInt32(6)) + (seed >> UInt32(2)))
        end
        return seed
    end
    ```
 * @tparam T should be primitive type std::uint[x]_t
 * @param vec Vector
 * @return a hash of all vector entries
 */
    template<typename T>
    std::uint32_t hash_vector(const std::vector<T> &vec) {
        auto seed = static_cast<std::uint32_t>(vec.size());
        for (auto i : vec) {
            seed ^= static_cast<std::uint32_t>(i) + 0x9e3779b9 + (seed << 6u) + (seed >> 2u);
        }
        return seed;
    }


    template<typename T>
    void write_vector_to_csv(const std::string &filepath, const std::vector<T> &vec, bool cast_to_long) {
        std::ofstream myfile(filepath);
        for (std::size_t n = 0; n < vec.size(); n++) {
            if (cast_to_long) {
                myfile << static_cast<long>(vec[n]) << '\n';
            } else {
                myfile << vec[n] << '\n';
            }
        }
    }


    template<typename T>
    void print_arr(T *p, std::size_t len) {
        std::cout << std::endl;
        for (int i = 0; i < len; ++i) {
            std::cout << p[i] << ' ';
        }
        std::cout << std::endl << std::endl;
    }


    template<typename T>
    std::ostream &operator<<(std::ostream &s, const std::vector<T> &v) {
        for (T i : v)
            s << i;
        return s;
    }

    template<typename Bit=Bit>
    std::vector<Bit> get_bitstring(std::size_t n) {
        std::vector<Bit> initial{0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1};

        std::vector<Bit> toencode(n);
        for (std::size_t i = 0; i < std::min(initial.size(), n); ++i)
            toencode[i] = static_cast<bool>(initial.at(i));

        return toencode;
    }


    template<typename T, std::size_t N>
    void noise_bitstring_inplace(std::array<T, N> &src, double err_prob, unsigned int seed = 0) {
        std::mt19937_64 rng{seed}; // hard-coded seed for testing purposes.

        std::bernoulli_distribution distribution(err_prob);

        for (std::size_t i = 0; i < src.size(); i++) {
            if (distribution(rng)) {
                src[i] = !src[i];
            } else {
                src[i] = src[i];
            }
        }
    }


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


    template<typename T>
    void noise_bitstring_inplace(std::vector<T> &src, double err_prob, unsigned int seed = 0) {
        std::mt19937_64 rng{seed}; // hard-coded seed for testing purposes.

        std::bernoulli_distribution distribution(err_prob);

        for (std::size_t i = 0; i < src.size(); i++) {
            if (distribution(rng)) {
                src[i] = !src[i];
            } else {
                src[i] = src[i];
            }
        }
    }

    template<typename T>
    std::uint32_t hash_vector(const T &vec) {
        auto seed = static_cast<std::uint32_t>(vec.size());
        for (auto i : vec) {
            seed ^= static_cast<std::uint32_t>(i) + 0x9e3779b9 + (seed << 6u) + (seed >> 2u);
        }
        return seed;
    }

    template<typename T, std::size_t N>
    std::vector<T> arr_to_vec(std::array<T, N> const &in) {
        std::vector<T> out(N);

        for (int i = 0; i < in.size(); ++i) {
            out[i] = in[i];
        }

        return out;
    }

    template<typename T, std::size_t N>
    void vec_to_arr(std::vector<T> const &in, std::array<T, N> &out) {
        if (out.size() == in.size()) {
            for (std::size_t i = 0; i < in.size(); ++i) {
                out[i] = in[i];
            }
        } else {
            std::cerr << "Warning!!! Size mismatch. Do nothing." << std::endl;
        }
    }


    template<class T>
    void print_nz_inds(const T &vec) {
        std::size_t i{};
        for (auto val : vec) {
            if (val != 0) {
                std::cout << i << ' ';
            }
            i++;
        }
        std::cout << std::endl;
    }

}

#endif //LDPC4QKD_HELPERS_FOR_TESTING_HPP
