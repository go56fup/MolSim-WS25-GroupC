#pragma once

#ifdef ENABLE_VTK_OUTPUT
#include "outputWriter/VTKWriter.h"
#define OUTPUT_WRITER outputWriter::VTKWriter
#else
#include "outputWriter/TXTWriter.h"
// #define OUTPUT_WRITER outputWriter::XYZWriter
#define OUTPUT_WRITER outputWriter::TXTWriter
#endif
