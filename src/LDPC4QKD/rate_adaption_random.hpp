//
// Created by Adomas Baliuka on 31.07.26.
//
// Defines Linear Congruential Generator (LCG) used to determine row indices to combine for rate adaption.
//

#ifndef LDPC4QKD_RATE_ADAPTION_RANDOM_HPP
#define LDPC4QKD_RATE_ADAPTION_RANDOM_HPP

#include <cstdint>
#include <stdexcept>
#include <numeric>

namespace LDPC4QKD {

    namespace RateAdaptLCG {

        /*!
         * Trial-division primality test. Only ever called with arguments up to a few million (mother matrix row
         * counts), so trial division (rather than a more sophisticated algorithm) is fast enough.
         */
        constexpr bool is_prime(std::uint64_t n) {
            if (n < 2) {
                return false;
            }
            if (n % 2 == 0) {
                return n == 2;
            }
            for (std::uint64_t d = 3; d * d <= n; d += 2) {
                if (n % d == 0) {
                    return false;
                }
            }
            return true;
        }

        //! Smallest prime >= n (matches Julia `Primes.nextprime`).
        constexpr std::uint64_t nextprime(std::uint64_t n) {
            if (n <= 2) {
                return 2;
            }
            std::uint64_t candidate = (n % 2 == 0) ? (n + 1) : n;
            while (!is_prime(candidate)) {
                candidate += 2;
            }
            return candidate;
        }

        //! Product of distinct prime factors of n (the "radical" / squarefree kernel), matches `Primes.radical`.
        constexpr std::uint64_t radical(std::uint64_t n) {
            std::uint64_t result = 1;
            for (std::uint64_t p = 2; p * p <= n; ++p) {
                if (n % p == 0) {
                    result *= p;
                    while (n % p == 0) {
                        n /= p;
                    }
                }
            }
            if (n > 1) {
                result *= n;
            }
            return result;
        }

        //! Smallest value >= start that is coprime to m. Always terminates: gcd(m+1, m) == 1 for any m, so the
        //! search never scans more than `m` candidates (in practice far fewer, for the highly composite mother
        //! matrix sizes used here).
        constexpr std::uint64_t find_coprime_at_least(std::uint64_t start, std::uint64_t m) {
            std::uint64_t c = start;
            while (std::gcd(c, m) != 1) {
                ++c;
            }
            return c;
        }

        /*!
         * Linear congruential generator (LCG) with period exactly `m` (assuming Hull-Dobell conditions hold,
         * checked in `get_LCG_with_period`).
         * See https://en.wikipedia.org/wiki/Linear_congruential_generator#c_%E2%89%A0_0
         */
        struct LCG {
            std::uint64_t seed;
            std::uint64_t A;
            std::uint64_t C;
            std::uint64_t m;

            //! Advances the generator and returns the new value (matches Julia's `next!`, which updates then returns).
            //! `A, seed < m` and mother matrices have nowhere near 2^32 rows in practice, so `A * seed` (< ~1e16
            //! for `m` up to 1e8) fits comfortably in `std::uint64_t` (max ~1.8e19) without needing a wider type.
            constexpr std::uint64_t next() {
                seed = (A * seed + C) % m;
                return seed;
            }
        };

        /*!
         * Get a linear congruential generator with period exactly `m`, for any `m >= 1`.
         * Chooses `C` coprime to `m` and `A` such that `(A - 1)` is divisible by all prime factors of `m`
         * (and by 4 if `m` is divisible by 4), satisfying the Hull-Dobell theorem for full period.
         *
         * `C` is picked by searching forward from `nextprime(m/4)` for the first value coprime to `m` (see
         * `find_coprime_at_least`). This is a strict generalization of just using `nextprime(m/4)` directly:
         * whenever that candidate already happens to be coprime to `m` (true for every mother matrix size used
         * in this codebase) the search terminates immediately and returns exactly that value; it only searches
         * further for the rare `m` where the plain `nextprime(m/4)` choice would not have given full period
         * (e.g. `m == 2` or `m == 4`). Unlike an earlier version of this function, this construction always
         * succeeds -- there is no longer a degenerate `m` for which it throws.
         *
         * @param m period of the generator. Must be >= 1.
         * @param seed initial state of the generator.
         */
        constexpr LCG get_LCG_with_period(std::uint64_t m, std::uint64_t seed = 0) {
            if (m < 1) {
                throw std::domain_error("get_LCG_with_period: `m` must be positive.");
            }
            if (m == 1) {
                // degenerate but well-defined: the only residue mod 1 is 0.
                return LCG{0, 1, 0, 1};
            }

            const std::uint64_t C = find_coprime_at_least(nextprime(m / 4), m);

            std::uint64_t A_minus_1 = radical(m);
            if (m % 4 == 0 && (A_minus_1 % 4) != 0) {
                A_minus_1 *= 4;
            }
            const std::uint64_t A = A_minus_1 + 1;

            // gcd(C, m) == 1 holds unconditionally by construction of `C` above.
            return LCG{seed % m, A, C, m};
        }

    } // namespace detail

} // namespace LDPC4QKD

#endif //LDPC4QKD_RATE_ADAPTION_RANDOM_HPP
