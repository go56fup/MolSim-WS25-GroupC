#pragma once


#include <omp.h>


inline int get_thread_num() {
#if SINGLETHREADED
	return -1;
#else
	return omp_get_thread_num();
#endif
}

inline int get_max_threads() {
#if SINGLETHREADED
	return 1;
#else
	return omp_get_max_threads();
#endif
}

inline int get_num_threads() {
#if SINGLETHREADED
	return 1;
#else
	return omp_get_num_threads();
#endif
}

inline void set_num_threads() {
	return omp_set_num_threads(2);;
}


inline void debug_omp() {
	omp_set_num_threads(4);
	spdlog::info("--- OpenMP Debug ---");
	spdlog::info("Num threads (omp_get_num_threads): {}", omp_get_num_threads());
	spdlog::info("Max threads (omp_get_max_threads): {}", get_max_threads());
	spdlog::info("Num procs available: {}", omp_get_num_procs());
	spdlog::info("Thread limit (omp_get_thread_limit): {}", omp_get_thread_limit());
	spdlog::info("Nested level (omp_get_level): {}", omp_get_level());
	spdlog::info("Is dynamic adjustment enabled: {}", omp_get_dynamic());
}