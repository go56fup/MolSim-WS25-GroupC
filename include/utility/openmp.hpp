#pragma once

#if !SINGLETHREADED
#include <omp.h>
#endif

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
