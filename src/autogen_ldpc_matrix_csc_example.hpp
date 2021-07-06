// This is an automatically generated file.
// A sparse LDPC matrix (containing only zeros and ones) is saved in compressed sparse column (CSC) format.
// Since the matrix (and LDPC code) is known at compile time, there is no need to save it separately in a file.
// This significantly blows up the executable size (the memory would still have to be used when saving the matrix).
// The method seems to be reasonably fast (on a standard laptop).

#include <cstdint>
#include <array>

namespace AutogenLDPC {

    constexpr std::size_t M = 5;
    constexpr std::size_t N = 10;
    constexpr std::size_t num_nz = 5;

    constexpr std::array<std::uint16_t, N + 1> colptr = {
            0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x5, 0x5, 0x5, 0x5, 0x5
    };

// ------------------------------------------------------- 

    constexpr std::array<std::uint16_t, num_nz> row_idx = {
            0x0, 0x1, 0x2, 0x3, 0x4
    };


} // namespace RALDPC
