#pragma once

#include <array>
#include <cstdint>

inline constexpr std::size_t batch_size = 8;

template <typename Value>
using batch = std::array<Value, batch_size>;

using particle_id = std::size_t;
using particle_batch = batch<particle_id>;
