#pragma once

#include <spdlog/spdlog.h>

#include "utility/tracing/config.hpp"

#if LOG_FORCES
#define TRACE_FORCES(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_FORCES(...) (void)0
#endif

#if LOG_INTERACTION_ITER
#define TRACE_INTERACTION_ITER(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_INTERACTION_ITER(...) (void)0
#endif

#if LOG_PARTICLE_CONTAINER
#define TRACE_PARTICLE_CONTAINER(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_PARTICLE_CONTAINER(...) (void)0
#endif

#if LOG_BORDER_CELL_ITER
#define TRACE_BORDER_CELL_ITER(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_BORDER_CELL_ITER(...) (void)0
#endif

#if LOG_SPECIAL_MEMFUNS
#define TRACE_SPECIAL_MEMFUNS(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_SPECIAL_MEMFUNS(...) (void)0
#endif

#if LOG_SIM
#define TRACE_SIM(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_SIM(...) (void)0
#endif

#if LOG_SIM_STATE
#define TRACE_SIM_STATE(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_SIM_STATE(...) (void)0
#endif

#if LOG_GRID
#define TRACE_GRID(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_GRID(...) (void)0
#endif

#if LOG_INPUT_PARSING
#define TRACE_INPUT_PARSING(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_INPUT_PARSING(...) (void)0
#endif

#if LOG_PERIODIC
#define TRACE_PERIODIC(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_PERIODIC(...) (void)0
#endif

#if LOG_THERMOSTAT
#define TRACE_THERMOSTAT(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_THERMOSTAT(...) (void)0
#endif

#if LOG_CHECKPOINT
#define TRACE_CHECKPOINT(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_CHECKPOINT(...) (void)0
#endif

#if LOG_RANDOM
#define TRACE_RANDOM(...) SPDLOG_TRACE(__VA_ARGS__)
#else
#define TRACE_RANDOM(...) (void)0
#endif
