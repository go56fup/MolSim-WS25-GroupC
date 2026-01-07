#pragma once

#include "simulation/config/parse.hpp"

// TODO(tuna): move to cmakepresets
#define ENABLE_CHECKPOINT

#ifdef ENABLE_VTK_OUTPUT
#include "output_writers/vtk.hpp"
#define WRITE_VTK_OUTPUT(plot_f, ...) plot_f(output_writer::vtk::plot_particles, __VA_ARGS__)
#else
#define WRITE_VTK_OUTPUT(...) (void)0
#endif

#ifdef ENABLE_CHECKPOINT
#define WRITE_CHECKPOINT(container, output_path)                                                   \
	config::write_state_to_file(config::dump_state(container), output_path)
#else
#define WRITE_CHECKPOINT(...) (void)0
#endif

// TODO(tuna): fix docs
/**
 * @brief Concept that defines a particle I/O provider.
 *
 * A callable satisfies this concept if it can be invoked as-if it shares this signature:
 * @verbatim
 void operator()(std::span<const Particle> particles, std::string_view base_name, int iteration)
 @endverbatim
 * where @p particles is the particles to output, @p base_name is the base name for the output files
 being
 * generated and @p iteration is the current simulation iteration.
 *
 * @tparam Candidate Type being tested for the particle I/O provider interface.
 *
 * @return `true` @a iff the callable can be invoked as specified.
 */
template <typename Candidate>
concept particle_io_provider = requires(
	Candidate f, particle_container& particles, std::string_view out_name, unsigned iteration
) {
	{ std::invoke(f, particles, out_name, iteration) } -> std::same_as<void>;
};
