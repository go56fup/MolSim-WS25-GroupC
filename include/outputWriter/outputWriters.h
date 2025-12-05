#pragma once

#ifdef ENABLE_VTK_OUTPUT
#include "outputWriter/VTKWriter.h"
#define OUTPUT_WRITER outputWriter::VTKWriter
#else
#include "outputWriter/XYZWriter.h"
#define OUTPUT_WRITER outputWriter::XYZWriter
#endif
