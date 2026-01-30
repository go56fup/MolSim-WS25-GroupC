#pragma once

#include <gtest/gtest.h>

#include "constexpr_comparators.hpp"

// Prefix all gtest implementation info with
// https://github.com/google/googletest/blob/b2b9072/googletest/
// to get direct links to GitHub.

#ifndef GTEST_CONFIG_RUNTIME_STATIC_EXPECT
#define GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_OP2(op, gtest_macro, a, b)                                         \
	static_assert(a op b, #gtest_macro "(" #a ", " #b ")");                                                            \
	SUCCEED() << #gtest_macro "(" #a ", " #b ")"

#define GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_FUNC(f, gtest_macro, ...)                                          \
	static_assert(f(__VA_ARGS__), #gtest_macro "(" #__VA_ARGS__ ")");                                                  \
	SUCCEED() << #gtest_macro "(" #__VA_ARGS__ ")"
#define GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_PREDN(f, gtest_macro, ...)                                         \
	static_assert(f(__VA_ARGS__), #gtest_macro "(" #f ", " #__VA_ARGS__ ")");                                          \
	SUCCEED() << #gtest_macro "(" #f ", " #__VA_ARGS__ ")"
#else
#define GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_OP2(op, gtest_macro, a, b) gtest_macro(a, b)
#define GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_FUNC(f, gtest_macro, ...) gtest_macro(__VA_ARGS__)
#define GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_PREDN(f, gtest_macro, ...) gtest_macro(f, __VA_ARGS__)
#endif

#define GTEST_CONSTEXPR_INTERNAL_STAMP_ASSERTION_FUNC(type, f, gtest_suffix, ...)                                      \
	GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_FUNC(f, type##_##gtest_suffix, __VA_ARGS__)

// These are generated via a macro:
// include/gtest/gtest.h#L1460
// The macro just evaluates val1 op val2, so skip the middleman and just
// static_assert(val1 op val2).
#define GTEST_CONSTEXPR_INTERNAL_STAMP_COMPARISON(type, op, name, ...)                                                 \
	GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_OP2(op, type##_##name, __VA_ARGS__)

#define GTEST_CONSTEXPR_INTERNAL_STAMP_TRUE(type, ...)                                                                 \
	GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_FUNC(, type##_TRUE, __VA_ARGS__)
#define GTEST_CONSTEXPR_INTERNAL_STAMP_FALSE(type, ...)                                                                \
	GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_FUNC(!, type##_FALSE, __VA_ARGS__)

// Google's implementation and this are not strictly equivalent, but avoiding
// the overload resolution results in better output in the form of "expression
// evaluates to '1 == 2'" which makes issues with assertions much clearer.
#ifndef GTEST_CONFIG_EQ_NO_SHORTCUT
// EQ:
// ASSERT: include/gtest/gtest.h#L1900
// EXPECT: include/gtest/gtest.h#L1887

#define GTEST_CONSTEXPR_INTERNAL_STAMP_EQ(type, ...)                                                                   \
	GTEST_CONSTEXPR_INTERNAL_STAMP_COMPARISON(type, ==, EQ, __VA_ARGS__)
#else
#define GTEST_CONSTEXPR_INTERNAL_STAMP_EQ(type, ...)                                                                   \
	GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_FUNC(gtest_constexpr::EqHelper::Compare, type##_EQ, __VA_ARGS__)
#endif

// NE:
// ASSERT: include/gtest/gtest.h#L1889
// EXPECT: include/gtest/gtest.h#L1921
#define GTEST_CONSTEXPR_INTERNAL_STAMP_NE(type, ...)                                                                   \
	GTEST_CONSTEXPR_INTERNAL_STAMP_COMPARISON(type, !=, NE, __VA_ARGS__)

// LE:
// ASSERT: include/gtest/gtest.h#L1925
// EXPECT: include/gtest/gtest.h#L1891
#define GTEST_CONSTEXPR_INTERNAL_STAMP_LE(type, ...)                                                                   \
	GTEST_CONSTEXPR_INTERNAL_STAMP_COMPARISON(type, <=, LE, __VA_ARGS__)

// LT:
// ASSERT: include/gtest/gtest.h#L1929
// EXPECT: include/gtest/gtest.h#L1893
#define GTEST_CONSTEXPR_INTERNAL_STAMP_LT(type, ...) GTEST_CONSTEXPR_INTERNAL_STAMP_COMPARISON(type, <, LT, __VA_ARGS__)

// GE:
// ASSERT: include/gtest/gtest.h#L1933
// EXPECT: include/gtest/gtest.h#L1895
#define GTEST_CONSTEXPR_INTERNAL_STAMP_GE(type, ...)                                                                   \
	GTEST_CONSTEXPR_INTERNAL_STAMP_COMPARISON(type, >=, GE, __VA_ARGS__)

// GT:
// ASSERT: include/gtest/gtest.h#L1937
// EXPECT: include/gtest/gtest.h#L1897
#define GTEST_CONSTEXPR_INTERNAL_STAMP_GT(type, ...) GTEST_CONSTEXPR_INTERNAL_STAMP_COMPARISON(type, >, GT, __VA_ARGS__)

#define GTEST_CONSTEXPR_INTERNAL_STAMP_VIA_CMPHELPER(type, cmphelper, ...)                                             \
	GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_FUNC(                                                                  \
		gtest_constexpr::CmpHelper##cmphelper, ASSERT_##cmphelper, __VA_ARGS__                                         \
	)

// STREQ:
// ASSERT: include/gtest/gtest.h#L1965
// EXPECT: include/gtest/gtest.h#L1956
#define GTEST_CONSTEXPR_INTERNAL_STAMP_STREQ(type, ...)                                                                \
	GTEST_CONSTEXPR_INTERNAL_STAMP_VIA_CMPHELPER(type, STREQ, __VA_ARGS__)

// STRNE:
// ASSERT: include/gtest/gtest.h#L1967
// EXPECT: include/gtest/gtest.h#L1958
#define GTEST_CONSTEXPR_INTERNAL_STAMP_STRNE(type, ...)                                                                \
	GTEST_CONSTEXPR_INTERNAL_STAMP_VIA_CMPHELPER(type, STRNE, __VA_ARGS__)

// STRCASEEQ:
// ASSERT: include/gtest/gtest.h#L1969
// EXPECT: include/gtest/gtest.h#L1960
#define GTEST_CONSTEXPR_INTERNAL_STAMP_STRCASEEQ(type, ...)                                                            \
	GTEST_CONSTEXPR_INTERNAL_STAMP_VIA_CMPHELPER(type, STRCASEEQ, __VA_ARGS__)

// STRCASENE:
// ASSERT: include/gtest/gtest.h#L1971
// EXPECT: include/gtest/gtest.h#L1962
#define GTEST_CONSTEXPR_INTERNAL_STAMP_STRCASENE(type, ...)                                                            \
	GTEST_CONSTEXPR_INTERNAL_STAMP_VIA_CMPHELPER(type, STRCASENE, __VA_ARGS__)

#define GTEST_CONSTEXPR_INTERNAL_STAMP_FLOATS(type, float_type, gtest_macro, ...)                                      \
	GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_FUNC(                                                                  \
		gtest_constexpr::CmpHelperFloatingPointEQ<float_type>, type##_##gtest_macro, __VA_ARGS__                       \
	)

// FLOAT_EQ:
// ASSERT: include/gtest/gtest.h#L1996
// EXPECT: include/gtest/gtest.h#L1988
#define GTEST_CONSTEXPR_INTERNAL_STAMP_FLOAT_EQ(type, ...)                                                             \
	GTEST_CONSTEXPR_INTERNAL_STAMP_FLOATS(type, float, FLOAT_EQ, __VA_ARGS__)

// DOUBLE_EQ:
// ASSERT: include/gtest/gtest.h#L2000
// EXPECT: include/gtest/gtest.h#L1992
#define GTEST_CONSTEXPR_INTERNAL_STAMP_DOUBLE_EQ(type, ...)                                                            \
	GTEST_CONSTEXPR_INTERNAL_STAMP_FLOATS(type, double, DOUBLE_EQ, __VA_ARGS__)

// NEAR:
// ASSERT: include/gtest/gtest.h#L2008
// EXPECT: include/gtest/gtest.h#L2004
#define GTEST_CONSTEXPR_INTERNAL_STAMP_NEAR(type, ...)                                                                 \
	GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_FUNC(gtest_constexpr::DoubleNearPredFormat, type##_NEAR, __VA_ARGS__)

#define GTEST_CONSTEXPR_INTERNAL_STAMP_PREDN(type, n, f, ...)                                                          \
	GTEST_CONSTEXPR_INTERNAL_GENERATE_ASSERTION_PREDN(f, type##_PRED##n, __VA_ARGS__)
// PRED1:
// ASSERT: include/gtest/gtest_pred_impl.h#L113
// EXPECT: include/gtest/gtest_pred_impl.h#L110
#define GTEST_CONSTEXPR_INTERNAL_STAMP_PRED1(type, ...) GTEST_CONSTEXPR_INTERNAL_STAMP_PREDN(type, 1, __VA_ARGS__)

// PRED2:
// ASSERT: include/gtest/gtest_pred_impl.h#L149
// EXPECT: include/gtest/gtest_pred_impl.h#L145
#define GTEST_CONSTEXPR_INTERNAL_STAMP_PRED2(type, ...) GTEST_CONSTEXPR_INTERNAL_STAMP_PREDN(type, 2, __VA_ARGS__)

// PRED3:
// ASSERT: include/gtest/gtest_pred_impl.h#L188
// EXPECT: include/gtest/gtest_pred_impl.h#L184
#define GTEST_CONSTEXPR_INTERNAL_STAMP_PRED3(type, ...) GTEST_CONSTEXPR_INTERNAL_STAMP_PREDN(type, 3, __VA_ARGS__)

// PRED4:
// ASSERT: include/gtest/gtest_pred_impl.h#L229
// EXPECT: include/gtest/gtest_pred_impl.h#L225
#define GTEST_CONSTEXPR_INTERNAL_STAMP_PRED4(type, ...) GTEST_CONSTEXPR_INTERNAL_STAMP_PREDN(type, 4, __VA_ARGS__)

// PRED5:
// ASSERT: include/gtest/gtest_pred_impl.h#L274
// EXPECT: include/gtest/gtest_pred_impl.h#L270
#define GTEST_CONSTEXPR_INTERNAL_STAMP_PRED5(type, ...) GTEST_CONSTEXPR_INTERNAL_STAMP_PREDN(type, 5, __VA_ARGS__)
