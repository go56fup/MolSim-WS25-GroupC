#pragma once

#ifdef ENABLE_VTK_OUTPUT
#include "output_writers/vtk.hpp"
#define WRITE_VTK_OUTPUT(...) output_writer::vtk::plot_particles(__VA_ARGS__)
#else
#define WRITE_VTK_OUTPUT(...) (void)0
#endif
