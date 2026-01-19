#pragma once

#if !defined(__clang__) && defined(__GNUC__)
#define IS_GCC 1
#else
#define IS_GCC 0
#endif

/**
 * @brief Macro expanding to `constexpr` only when compiling with GCC.
 *
 * Constexpr behavior is different between GCC and clang, specifically with `<cmath>` and the
 * accompanying `__builtin_sqrt` etc. This means some mathematical calculations can be done
 * in compile-time with GCC but not clang, where these macros come in handy.
 *
 * @note <a href="https://github.com/llvm/llvm-project/issues/167874">This issue</a> tracks Clang
 * support for the above functions.
 */
#if IS_GCC
#define CONSTEXPR_IF_GCC constexpr
#else
#define CONSTEXPR_IF_GCC
#endif

#if __cpp_pp_embed && !HAS_EMBED
#error "#embed support is claimed by the compiler but was not detected by CMake."
#endif

#if defined(__clang__)
#define IS_CLANG 1
#else
#define IS_CLANG 0
#endif
