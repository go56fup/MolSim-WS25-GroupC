#pragma once

#include "gtest_constexpr/macros.hpp"

#include "utility/compiler_traits.hpp"
#include "utility/debug.hpp"

#include "physics/vec_3d.hpp"

// TODO(tuna): change the name of this file to something more descriptive along the lines of
// common_something etc., move the GCC_* asserts and GTEST_CXP_GCC here
template class trace_special_memfuns<"particle">;
template class trace_special_memfuns<"vec_3d">;
template struct vec_3d<double>;

#define INTERNAL_COMPONENT_CMP(i, comparison, a, b, ...)                                           \
	comparison((a).i, (b).i __VA_OPT__(, ) __VA_ARGS__)
#define INTERNAL_VEC_COMPONENTWISE_CMP(...)                                                        \
	INTERNAL_COMPONENT_CMP(x, __VA_ARGS__);                                                        \
	INTERNAL_COMPONENT_CMP(y, __VA_ARGS__);                                                        \
	INTERNAL_COMPONENT_CMP(z, __VA_ARGS__)

#define GCC_STATIC_EXPECT_VEC_NEAR(a, b, abs_error)                                                \
	INTERNAL_VEC_COMPONENTWISE_CMP(GCC_STATIC_EXPECT_NEAR, a, b, abs_error)
#define STATIC_EXPECT_VEC_DOUBLE_EQ(a, b)                                                          \
	INTERNAL_VEC_COMPONENTWISE_CMP(STATIC_EXPECT_DOUBLE_EQ, a, b)
#define GCC_STATIC_EXPECT_VEC_DOUBLE_EQ(a, b)                                                      \
	INTERNAL_VEC_COMPONENTWISE_CMP(GCC_STATIC_EXPECT_DOUBLE_EQ, a, b)
#define EXPECT_VEC_DOUBLE_EQ(a, b) INTERNAL_VEC_COMPONENTWISE_CMP(EXPECT_DOUBLE_EQ, a, b)

// TODO(tuna): maybe not only bool but also implicitly convertible to bool, also maybe constrain the
// auto in Results
template <auto Results>
void expect_all() {
	[&]<std::size_t... I>(std::index_sequence<I...>) {
		(..., ([&] { STATIC_EXPECT_TRUE(Results[I]); }()));
	}(std::make_index_sequence<std::tuple_size_v<decltype(Results)>>{});
}

template <std::size_t N>
void expect_all(const std::array<bool, N>& results) {
	for (std::size_t i = 0; i < results.size(); ++i) {
		EXPECT_TRUE(results[i]) << "testing at index: " << i;
	}
}

#ifndef GTEST_CONFIG_RUNTIME_STATIC_EXPECT
#define STATIC_EXPECT_ALL(arr) expect_all<arr>()
#else
#define STATIC_EXPECT_ALL(arr) expect_all(arr)
#endif

#if IS_GCC
#define GCC_STATIC_EXPECT_ALL(arr) STATIC_EXPECT_ALL(arr)
#else
#define GCC_STATIC_EXPECT_ALL(arr) expect_all(arr)
#endif

#if IS_GCC
#define GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(gtest_macro) STATIC_##gtest_macro
#else
#define GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(gtest_macro) gtest_macro
#endif

/**
 * @name GCC_STATIC_{ASSERT,EXPECT}_*
 * @{
 * @brief Assertion macros that are only STATIC when compiling with GCC.
 *
 * For rationale, see @ref CONSTEXPR_IF_GCC.
 */

/** @gtest_equ_gcc_val{ASSERT_EQ} */
#define GCC_STATIC_ASSERT_EQ GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(ASSERT_EQ)

/** @gtest_equ_gcc_val{EXPECT_EQ} */
#define GCC_STATIC_EXPECT_EQ GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(EXPECT_EQ)

/** @gtest_equ_gcc{ASSERT_TRUE,condition} */
#define GCC_STATIC_ASSERT_TRUE GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(ASSERT_TRUE)

/** @gtest_equ_gcc{EXPECT_TRUE,condition} */
#define GCC_STATIC_EXPECT_TRUE GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(EXPECT_TRUE)

/** @gtest_equ_gcc_val{ASSERT_DOUBLE_EQ} */
#define GCC_STATIC_ASSERT_DOUBLE_EQ GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(ASSERT_DOUBLE_EQ)

/** @gtest_equ_gcc_val{EXPECT_DOUBLE_EQ} */
#define GCC_STATIC_EXPECT_DOUBLE_EQ GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(EXPECT_DOUBLE_EQ)

/** @gtest_equ_gcc_val{ASSERT_NEAR} */
#define GCC_STATIC_ASSERT_NEAR GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(ASSERT_NEAR)

/** @gtest_equ_gcc_val{EXPECT_NEAR} */
#define GCC_STATIC_EXPECT_NEAR GTEST_CONSTEXPR_INTERNAL_STAMP_GCC(EXPECT_NEAR)
/// @}

/**
 * @brief `static constexpr` @a iff `GTEST_CONFIG_RUNTIME_STATIC_EXPECT` is not set.
 *
 * Used to push calculations to runtime, as `static constexpr` has to be opted out of.
 */
#ifndef GTEST_CONFIG_RUNTIME_STATIC_EXPECT
#define GTEST_CXP static constexpr
#else
#define GTEST_CXP
#endif

/**
 * @brief Equivalent to `GTEST_CXP` when compiling with GCC.
 *
 * For rationale, see @ref CONSTEXPR_IF_GCC.
 */
#if IS_GCC
#define GTEST_CXP_GCC GTEST_CXP
#else
#define GTEST_CXP_GCC
#endif

#if COMPILE_SIM_TESTS
#define GTEST_CXP_SIM GTEST_CXP
#else
#define GTEST_CXP_SIM
#endif

#if IS_GCC
#define GTEST_CXP_GCC_SIM GTEST_CXP_SIM
#else
#define GTEST_CXP_GCC_SIM
#endif

namespace detail {
template <std::size_t MaxSize, typename ValueT, typename Builder>
constexpr std::span<const ValueT> put_into_static_storage(Builder&& builder) {
	static auto data = std::invoke(std::forward<Builder>(builder));
	return data;
}
}

#ifndef GTEST_CONFIG_RUNTIME_STATIC_EXPECT
#define CONSTEXPR_TWO_STEP constexpr_two_step
#else
#define CONSTEXPR_TWO_STEP detail::put_into_static_storage
#endif
