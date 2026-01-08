#pragma once
#include "grid/particle_container/fwd.hpp"

namespace detail {

    class periodic_iterator {
    private:
        using difference_type_priv = particle_container::difference_type;
        using signed_index_t = vec_3d<difference_type_priv>;


        particle_container* container = nullptr;
        signed_index_t current_virtual_idx{};
        signed_index_t target_cell_idx{};
        std::uint8_t displacement_idx = 0;

        static constexpr std::uint8_t displacement_count = 13;
        vec_3d<difference_type_priv> signed_grid;

    public:
        // Standard iterator boilerplate for C++ ranges
            using iterator_category = std::forward_iterator_tag;
        using value_type = particle_container::index;
        using difference_type = particle_container::difference_type;
        using pointer = value_type*;
        using reference = value_type;

        // Start Constructor
        constexpr periodic_iterator(particle_container& c, vec_3d<difference_type> v_idx)
                : container(&c)
                , current_virtual_idx(v_idx)
                , target_cell_idx(v_idx) // Initialize this!
                , signed_grid{static_cast<difference_type_priv>(c.grid_size().x),
                              static_cast<difference_type_priv>(c.grid_size().y),
                              static_cast<difference_type_priv>(c.grid_size().z)}
        {
            ++*this;
        }

        // End Constructor
        constexpr periodic_iterator(particle_container& c, interactions_iterator::get_end_tag)
                : container(&c)
                , current_virtual_idx{} // Required for constexpr
                , target_cell_idx{}     // Required for constexpr
                , displacement_idx(displacement_count)
                , signed_grid{}         // Required for constexpr
        {}
        template <axis Axis>
        constexpr bool do_displacement(const vec_3d<difference_type>& displacement) noexcept {

            if (displacement[Axis] == 0) return true;
            TRACE_INTERACTION_ITER("Doing displacement {} on {}", displacement, current_virtual_idx);

            // We do not need to use __builtin_add_overflow() here, because:
            assert(displacement[Axis] == -1 || displacement[Axis] == 1);
            target_cell_idx[Axis] = current_virtual_idx[Axis] + displacement[Axis];

            // If we are moving backwards, don't go beyond the first ghost layer
            if (displacement[Axis] < 0 && target_cell_idx[Axis] < 0) {
                return false;
            }

            const bool outside_domain_bc_max =
                    out_of_bounds<axis_to_boundary_type(Axis, extremum::max)>(
                            target_cell_idx, signed_grid
                    );

            return !outside_domain_bc_max;
        }

        constexpr periodic_iterator& operator++() noexcept {
            // Try all displacements for the current cell
            while (displacement_idx < displacement_count) {
                const auto& displacement = displacements[displacement_idx];
                target_cell_idx = current_virtual_idx;

                const bool success = do_displacement<axis::x>(displacement) &&
                                     do_displacement<axis::y>(displacement) &&
                                     do_displacement<axis::z>(displacement);

                ++displacement_idx;

                if (success) {
                    TRACE_INTERACTION_ITER(
                            "Found valid interaction: {} -> {}", current_cell_idx, target_cell_idx
                    );
                    return *this;
                }
            }

            return *this;
        }

        constexpr periodic_iterator operator++(int) noexcept {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        constexpr bool operator==(const periodic_iterator& other) const noexcept {
            assert(container == other.container);
            return displacement_idx == other.displacement_idx;
        }

        constexpr bool operator!=(const periodic_iterator& other) const noexcept {
            return displacement_idx != other.displacement_idx;
        }

        constexpr particle_container::index operator*() const noexcept {
            using uint_t = particle_container::size_type;
            return {
                    static_cast<uint_t>(target_cell_idx.x),
                    static_cast<uint_t>(target_cell_idx.y),
                    static_cast<uint_t>(target_cell_idx.z)
            };
        }

        static constexpr std::array<vec_3d<difference_type>, 13> displacements = {
                {{0, 0, +1},    // i,     j,     k + 1
                 {0, +1, -1},   // i,     j + 1, k - 1
                 {0, +1, 0},    // i,     j + 1, k
                 {0, +1, +1},   // i,     j + 1, k + 1
                 {1, -1, -1},   // i + 1, j - 1, k - 1
                 {+1, -1, 0},   // i + 1, j - 1, k
                 {+1, -1, +1},  // i + 1, j - 1, k + 1
                 {+1, 0, -1},   // i + 1, j,     k - 1
                 {+1, 0, 0},    // i + 1, j,     k
                 {+1, 0, +1},   // i + 1, j,     k + 1
                 {+1, +1, -1},  // i + 1, j + 1, k - 1
                 {+1, +1, 0},   // i + 1, j + 1, k
                 {+1, +1, +1}}  // i + 1, j + 1, k + 1
        };
    };
} // namespace detail

    /**
 * @brief A lightweight range for periodic neighbor interactions.
 * Input: A signed virtual index.
 * Output: Yields unsigned target indices.
 */
class periodic_range {
    particle_container& container;
    // We store the "Lie" here so begin() can use it
    vec_3d<particle_container::difference_type> v_idx;

public:
    constexpr periodic_range(particle_container& c,
                             vec_3d<particle_container::difference_type> start_v_idx) noexcept
            : container(c), v_idx(start_v_idx) {}

    constexpr detail::periodic_iterator begin() noexcept {
        return {container, v_idx};
    }

    constexpr detail::periodic_iterator end() noexcept {
        return {container, detail::interactions_iterator::get_end_tag{}};
    }
};

// This will tell you exactly which line of the iterator requirement is missing
static_assert(std::forward_iterator<detail::periodic_iterator>);