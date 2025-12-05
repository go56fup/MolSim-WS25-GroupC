#pragma once

#include "Debug.h"

// TODO(tuna): change the name of this file to something more descriptive along the lines of common_something
// etc., move the GCC_* asserts and GTEST_CXP_GCC here
template class flag_special_member_funcs<"Particle">;
template class flag_special_member_funcs<"vec_3d">;

#define INTERNAL_COMPONENT_CMP(i, comparison, a, b, ...) comparison((a).i, (b).i __VA_OPT__(, ) __VA_ARGS__)
#define INTERNAL_VEC_COMPONENTWISE_CMP(...)                                                                            \
	INTERNAL_COMPONENT_CMP(x, __VA_ARGS__);                                                                            \
	INTERNAL_COMPONENT_CMP(y, __VA_ARGS__);                                                                            \
	INTERNAL_COMPONENT_CMP(z, __VA_ARGS__)

#define GCC_STATIC_EXPECT_VEC_NEAR(a, b, abs_error)                                                                    \
	INTERNAL_VEC_COMPONENTWISE_CMP(GCC_STATIC_EXPECT_NEAR, a, b, abs_error)
#define GCC_STATIC_EXPECT_VEC_DOUBLE_EQ(a, b) INTERNAL_VEC_COMPONENTWISE_CMP(GCC_STATIC_EXPECT_DOUBLE_EQ, a, b)
#define EXPECT_VEC_DOUBLE_EQ(a, b) INTERNAL_VEC_COMPONENTWISE_CMP(EXPECT_DOUBLE_EQ, a, b)

// TODO(tuna): maybe not only bool but also implicitly convertible to bool, also maybe constrain the auto in
// Results
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

#ifdef IS_GCC
#define GCC_STATIC_EXPECT_ALL(arr) STATIC_EXPECT_ALL(arr)
#else
#define GCC_STATIC_EXPECT_ALL(arr) expect_all(arr)
#endif
