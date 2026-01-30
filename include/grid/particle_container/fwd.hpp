#pragma once

#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "grid/particle_container/system.hpp"
#include "physics/vec_3d.hpp"
#include "simulation/entities.hpp"
#include "utility/concepts.hpp"

// TODO(tuna): add particle_system related docs
/**
 * @brief The container type used to hold particles.
 *
 * This type is used throughout the simulation to store all particle objects.
 */
class particle_container {
private:
	static constexpr auto max_particle_type_count = 8;

public:
	/// The fundamental element of the grid structure of the simulation.
	using cell = std::vector<particle_id>;
	using value_type = cell;
	/// The type representing the index of a cell on an axis.
	// Not std::size_t because this type needs to operate nicely with positions of doubles, and
	// std::size_t indices cannot be fully represented in a double.
	using size_type = std::uint32_t;
	/// The type used to represent differences between indices on an axis.
	using difference_type = std::int32_t;
	/// The type used to index into this container.
	using index = vec_3d<size_type>;
	/// The type used to represent differences between indices.
	using signed_index = vec_3d<difference_type>;

private:
	/**
	 * @brief Determines the linear index of the cell of containing a position.
	 * @param pos Position in domain.
	 * @return Linear index of the cell containing the position.
	 */
	constexpr size_type pos_to_linear_index(const vec& pos) const noexcept;

	/**
	 * @brief Check if the given position runs afoul of domain maximum boundaries.
	 *
	 * @throws std::domain_error If the given position is outside the domain.
	 * @param pos Position to check against.
	 */
	constexpr void check_if_out_of_domain_max(const vec& pos) const noexcept(false);

	/// Storage for particles.
	std::vector<cell> grid;
	/// Domain size in terms of simulation coordinates. Used for bounds checking.
	vec domain_;
	/// Number of cells on each axis. Used for indexing.
	index grid_size_;
	/// Dimension of each cell (a cube). Forms the correspondence between the position and index
	/// world.
	double cutoff_radius_;
	particle_system system_;
	std::vector<particle_id> membrane_particles_;
	std::vector<particle_id> upwards_moving_membrane_members_;
	index membrane_scale_;

public:
	/**
	 * @brief Get linear index, associated with a cell's coordinates on the grid.
	 * @param x X coordinate.
	 * @param y Y coordinate.
	 * @param z Z coordinate.
	 * @return Index into the 1D representation of the grid.
	 */
	constexpr size_type linear_index(size_type x, size_type y, size_type z) const;

	constexpr particle_container::index index_from_linear(size_type linear) const;

	// TODO(tuna): make all constructors specify the object they are constructing in their @brief
	// docs
	/**
	 * @brief Parameterized constructor for a particle container.
	 *
	 * Initializes the domain, grid size and cutoff radius.
	 *
	 * @tparam Vec Deduced forwarding reference type for `domain_arg`.
	 * @param domain_arg Width, height and depth of the domain.
	 * @param cutoff_radius_arg Cutoff radius, which becomes the dimensions for each cell.
	 */
	template <fwd_reference_to<vec> Vec>
	constexpr particle_container(Vec&& domain_arg, double cutoff_radius_arg);

	/**
	 * @brief Constructs a particle inside the container.
	 *
	 * Takes a pack of arguments that can construct a particle, and constructs it inside its
	 * corresponding cell.
	 *
	 * @tparam Vec Deduced forwarding reference type for `position`.
	 * @tparam Args Pack of types of @ref Particle constructor arguments.
	 *
	 * @param position Position to emplace particle at.
	 * @param args Arguments to construct particle from.
	 */
	constexpr void
	add_particle(const vec& position, const vec& velocity, const material_description& material, const sim_configuration& config);

	constexpr void reload_particle_state(
		const particle_state_parameters& parameters, const vec& velocity,
		const material_description& material, const sim_configuration& config
	);

	/**
	 * @brief Add a mesh of particles, representing one body, into the container.
	 *
	 * @param origin Position of the lower left front-side corner.
	 * @param scale Number of particles in each direction.
	 * @param meshwidth Relative distance between two particles.
	 * @param velocity Initial additional velocity of each particle.
	 * @param mass Mass of one particle.
	 * @param browninan_mean Average velocity of the Brownian Motion.
	 * @param seq_no Mutable reference to a sequence number, used to achieve stateless random number
	 * generation.
	 */
	template <std::size_t N>
	constexpr void add_cuboid(
		const cuboid_parameters<3>& cuboid, const body_common_parameters& body_parameters,
		const vec& velocity, const material_description& material, const sim_configuration& config, std::size_t& seq_no
	);

	constexpr void add_membrane(
		const cuboid_parameters<3>& membrane, const body_common_parameters& body_parameters,
		const vec& velocity, const material_description& material, const sim_configuration& config, std::size_t& seq_no
	);

	/**
	 * @brief Add a 2D disc of particles into the container.
	 *
	 * @param center Position of the center of the disc.
	 * @param radius Radius of the created disc.
	 * @param meshwidth Relative distance between two particles.
	 * @param velocity Initial additional velocity of each particle.
	 * @param mass Mass of one particle.
	 * @param browninan_mean Average velocity of the Brownian Motion.
	 * @param seq_no Mutable reference to a sequence number, used to achieve stateless random number
	 * generation.
	 */

	constexpr void add_disc(
		const disc_parameters<3>& disc, const body_common_parameters& body_parameters,
		const vec& velocity, const material_description& material, const sim_configuration& config, std::size_t& seq_no
	);

	// TODO(tuna): see if the return type specified conflicts when the container is const
	/**
	 * @brief Get a range over pair of cell indices that are to interact according to Newton's third
	 * law.
	 * @return Range of pair of indices of interacting cells.
	 */
	constexpr range_of<std::pair<const index&, const index&>> auto
	directional_interactions() noexcept;

	// TODO(tuna): change this once the return type of the new border cell is fixed
	/**
	 * @brief Get a range over indices of border cells and their types.
	 * @return Range of pair of pair of border cell index and corresponding border type.
	 */
	constexpr /* range_of<border_cell_range::iterator::value_type> */ auto border_cells() noexcept;

	// The return type of this contains references to cells, which are conditionally const, which
	// cannot be expressed in here: so we leave it as auto.
	/**
	 * @brief Get a range over indices and the cells they reference.
	 * @return Range of indices and corresponding cells.
	 */
	constexpr auto enumerate_cells() noexcept;

	/**
	 * @brief Get a range over all cells in the container.
	 * @return Range to const over all cells in container.
	 */
	constexpr const range_of<const particle_container::cell> auto& cells() const noexcept;

	/**
	 * @brief Get a range over all cells in the container.
	 * @return Range to non-const over all cells in container.
	 */
	constexpr range_of<particle_container::cell> auto& cells() noexcept;

	/**
	 * @brief Indexes into the grid using the given index.
	 * @param three_d_index Index triple of requested cell.
	 * @return Mutable reference to cell with given index.
	 */
	constexpr cell& operator[](const index& three_d_index) noexcept;

	/**
	 * @brief Indexes into the grid using the given index.
	 * @param three_d_index Index triple of requested cell.
	 * @return Const reference to cell with given index.
	 */
	constexpr const cell& operator[](const index& three_d_index) const noexcept;

	/**
	 * @brief Indexes into the grid using the given index.
	 * @param x Index of the cell on the X axis.
	 * @param y Index of the cell on the Y axis.
	 * @param z Index of the cell on the Z axis.
	 * @return Mutable reference to cell with given index.
	 */
	constexpr cell& operator[](size_type x, size_type y, size_type z) noexcept;

	/**
	 * @brief Indexes into the grid using the given index.
	 * @param x Index of the cell on the X axis.
	 * @param y Index of the cell on the Y axis.
	 * @param z Index of the cell on the Z axis.
	 * @return Const reference to cell with given index.
	 */
	constexpr const cell& operator[](size_type x, size_type y, size_type z) const noexcept;

	/**
	 * @brief Calculates the cell containing the given position.
	 * @param position Position in domain.
	 * @return Mutable reference to cell containing the given position.
	 */
	constexpr cell& cell_containing(const vec& position) noexcept;

	/**
	 * @brief Calculates the cell containing the given position.
	 * @param position Position in domain.
	 * @return Const reference to cell containing the given position.
	 */
	constexpr const cell& cell_containing(const vec& pos) const noexcept;

	/**
	 * @brief Get the domain adhered to by this container.
	 * @return Vector describing the domain size.
	 */
	constexpr const vec& domain() const noexcept;

	/**
	 * @brief Get the size of the grid used by this container.
	 * @return Vector describing the grid size.
	 */
	constexpr const index& grid_size() const noexcept;

	/**
	 * @brief Get the size of the cells of the grid used by this container.
	 * @return Cutoff radius, i.e. cell dimensions.
	 */
	constexpr double cutoff_radius() const noexcept;

	constexpr const particle_system& system() const noexcept;
	constexpr particle_system& system() noexcept;

	constexpr std::span<const particle_id> membrane_particles() const noexcept;
	constexpr const index& membrane_scale() const noexcept;
	constexpr std::span<const particle_id> upwards_moving_membrane_members() const noexcept;
};
