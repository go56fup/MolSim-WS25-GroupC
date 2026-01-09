#include <functional>

#include <gtest/gtest.h>
#include "gtest_constexpr/macros.hpp"
#include "testing_utilities.hpp"

#include "physics/vec_3d.hpp"

// NOLINTBEGIN(*magic-numbers)

// Check that vector addition works
TEST(VectorTests, Addition) {
	GTEST_CXP vec a{1, 2, 3};
	GTEST_CXP vec b{4, 5, 6};
	GTEST_CXP vec result{5, 7, 9};
	STATIC_EXPECT_EQ(a + b, result);
}

// Check that vector subtraction works
TEST(VectorTests, Subtraction) {
	GTEST_CXP vec a{5, 7, 9};
	GTEST_CXP vec b{4, 5, 6};
	GTEST_CXP vec result{1, 2, 3};
	STATIC_EXPECT_EQ(a - b, result);
}

// Check that vector plus-equals works
TEST(VectorTests, AdditionAssignment) {
	GTEST_CXP vec result = std::invoke([] {
		vec a{1, 2, 3};
		vec b{4, 5, 6};
		a += b;
		return a;
	});
	GTEST_CXP vec expected{5, 7, 9};
	STATIC_EXPECT_EQ(result, expected);
}

// Check that vector minus-equals works
TEST(VectorTests, SubtractionAssignment) {
	GTEST_CXP vec result = std::invoke([] {
		vec a{4, 5, 6};
		vec b{1, 2, 3};
		a -= b;
		return a;
	});
	GTEST_CXP vec expected{3, 3, 3};
	STATIC_EXPECT_EQ(result, expected);
}

// Check that vector scaling works
TEST(VectorTests, Scalar) {
	GTEST_CXP vec a{5, 7, 9};
	GTEST_CXP vec result{10, 14, 18};
	STATIC_EXPECT_EQ(a * 2, result);
	STATIC_EXPECT_EQ(2 * a, result);
}

// Check that vector times-equals works
TEST(VectorTests, ScalarAssignment) {
	GTEST_CXP vec result = std::invoke([] {
		vec a{5, 7, 9};
		a *= 2;
		return a;
	});

	GTEST_CXP vec expected{10, 14, 18};
	STATIC_EXPECT_EQ(result, expected);
}

// Check that vector normalization works
TEST(VectorTests, Norm) {
	GTEST_CXP vec a{2, 3, 6};
	GTEST_CXP double result = 7;
	GCC_STATIC_EXPECT_EQ(a.euclidian_norm(), result);
}

// Check that vector equality works
TEST(VectorTests, Equal) {
	GTEST_CXP vec a{1, 2, 3};
	GTEST_CXP vec b{1, 2, 3};
	STATIC_EXPECT_EQ(a, b);
}

// NOLINTEND(*magic-numbers)
