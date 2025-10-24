#pragma once

namespace detail {
/**
 * @brief Prints debug messages composed of arbitrary arguments.
 *
 * This helper function prints all provided arguments to `std::cout`,
 * prefixed with a `[DEBUG]` tag and followed by a newline.
 *
 * Do not call directly, use the @ref DEBUG_LOG macro,
 * which enables or disables debug output depending on the build mode.
 *
 * Example usage:
 * @code
 * debug_print_impl("Iteration: ", 10, ", Velocity: ", v);
 * // Output:
 * // [DEBUG]    Iteration: 10, Velocity: (1.0, 0.5, 0.0)
 * @endcode
 *
 * @tparam Ts Variadic template parameter pack representing the types of the message components.
 * @param msgs The message components to stream to `std::cout` (via `operator<<`).
 */
template <typename... Ts>
void debug_print_impl(const Ts&... msgs) {
	std::cout << "[DEBUG]\t";
	(..., (std::cout << msgs)) << '\n';
}
}  // namespace detail

/**
 * @def DEBUG_LOG(...)
 * @brief Conditional debug logging macro.
 *
 * Prints messages using detail::debug_print_impl only when debugging is enabled.
 * If the program is compiled in release mode (i.e., `NDEBUG` is defined),
 * this macro expands to a no-op.
 *
 * @param ... Arguments to be forwarded to detail::debug_print_impl.
 */
#ifndef NDEBUG
#define DEBUG_LOG(...) detail::debug_print_impl(__VA_ARGS__)
#else
#define DEBUG_LOG(...) (void)0
#endif
