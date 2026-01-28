#pragma once

#ifndef SINGLETHREADED
#include <omp.h>
#endif

inline int get_thread_num() {
#ifdef SINGLETHREADED
	return -1;
#else
	return omp_get_thread_num();
#endif
}

