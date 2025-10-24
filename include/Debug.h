#pragma once

template <typename... Ts>
void debug_print_impl(const Ts&... msgs) {
	std::cout << "[DEBUG]\t";
	(..., (std::cout << msgs)) << '\n';
}

#ifndef NDEBUG
#define DEBUG_LOG(...) debug_print_impl(__VA_ARGS__)
#else
#define DEBUG_LOG(...) (void)0
#endif
