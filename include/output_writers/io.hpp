#pragma once

#ifdef ENABLE_VTK_OUTPUT
#include "output_writers/vtk.hpp"
#define WRITE_VTK_OUTPUT(...) output_writer::vtk::plot_particles(__VA_ARGS__)
#else
#define WRITE_VTK_OUTPUT(...) (void)0
#endif

#ifdef ENABLE_VTK_OUTPUT
#define ENABLE_STATISTICS_OUTPUT
#endif

#ifdef ENABLE_STATISTICS_OUTPUT
#include "output_writers/statistics_output.hpp"
#define WRITE_STATISTICS_OUTPUT(...) write_statistics_output(__VA_ARGS__)
#else
#define WRITE_STATISTICS_OUTPUT(...) (void)0
#endif
