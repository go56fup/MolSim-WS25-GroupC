#pragma once

#include <cmath>
#include <concepts>
#include <limits>

#include <gtest/gtest.h>

// Prefix all gtest implementation info with
// https://github.com/google/googletest/blob/b2b9072/googletest/
// to get direct links to GitHub.

namespace gtest_constexpr {
// CmpHelperEQ is just lhs == rhs, so elide that:
// include/gtest/gtest.h#L1394

// This is just here to allow for future extension: in case I ever get around to
// implementing support for C++26's custom static_assert messages, one can
// actually transform the gtest code 1-to-1 to get the exact runtime messages at
// compile-time failures as well.
#define CmpHelperEQ(val1, val2) val1 == val2

class EqHelper {
public:
	// This templatized version is for the general case.
	template <
		typename T1, typename T2,
		// Disable this overload for cases where one argument is a pointer
	    // and the other is the null pointer constant.
		typename std::enable_if<!std::is_integral<T1>::value || !std::is_pointer<T2>::value>::type* = nullptr>
	static constexpr bool Compare(const T1& lhs, const T2& rhs) {
		return CmpHelperEQ(lhs, rhs);
	}

	// With this overloaded version, we allow anonymous enums to be used
	// in {ASSERT|EXPECT}_EQ when compiled with gcc 4, as anonymous
	// enums can be implicitly cast to BiggestInt.
	//
	// Even though its body looks the same as the above version, we
	// cannot merge the two, as it will make anonymous enums unhappy.
	static constexpr bool Compare(testing::internal::BiggestInt lhs, testing::internal::BiggestInt rhs) {
		return CmpHelperEQ(lhs, rhs);
	}

	template <typename T>
	static constexpr bool Compare(
		// Handle cases where '0' is used as a null pointer literal.
		std::nullptr_t /* lhs */, T* rhs
	) {
		// We already know that 'lhs' is a null pointer.
		return CmpHelperEQ(static_cast<T*>(nullptr), rhs);
	}

#undef CmpHelperEQ
};

namespace String {
// src/gtest.cc#L1268
constexpr bool CStringEquals(const char* lhs, const char* rhs) {
	if (lhs == nullptr) return rhs == nullptr;

	if (rhs == nullptr) return false;

	// strcmp is not constexpr, so use std::string_view.
	return std::string_view{lhs} == rhs;
}

// Locales are not available at compile-time, and strcasecmp does not respect Unicode
// or locales to begin with, so just assume ASCII to ensure parity.
constexpr unsigned char constexpr_tolower(unsigned char c) {
	return (c >= 'A' && c <= 'Z') ? (c + ('a' - 'A')) : c;
}

// From: https://github.com/nmoinvaz/strcasecmp
constexpr int constexpr_strcasecmp(const char* s1, const char* s2) {
	unsigned char c1;
	unsigned char c2;
	
	do {
		c1 = static_cast<unsigned char>(*s1++);
        c2 = static_cast<unsigned char>(*s2++);
		if (c1 != c2) {
			c1 = constexpr_tolower(c1);
			c2 = constexpr_tolower(c2);
			if (c1 != c2) return static_cast<int>(c1) - static_cast<int>(c2);
		}
	} while (c1 != 0);
	return 0;
}

// src/gtest.cc#L2190
constexpr bool CaseInsensitiveCStringEquals(const char* lhs, const char* rhs) {
	if (lhs == nullptr) return rhs == nullptr;
	if (rhs == nullptr) return false;
	return constexpr_strcasecmp(lhs, rhs) == 0;
}

}  // namespace String

// src/gtest.cc#L1788
constexpr bool CmpHelperSTREQ(const char* lhs, const char* rhs) {
	return String::CStringEquals(lhs, rhs);
}

// src/gtest.cc#L1800
constexpr bool CmpHelperSTRCASEEQ(const char* lhs, const char* rhs) {
	return String::CaseInsensitiveCStringEquals(lhs, rhs);
}

// src/gtest.cc#L1812
constexpr bool CmpHelperSTRNE(const char* s1, const char* s2) {
	return !String::CStringEquals(s1, s2);
}

// src/gtest.cc#L1825
constexpr bool CmpHelperSTRCASENE(const char* s1, const char* s2) {
	return !String::CaseInsensitiveCStringEquals(s1, s2);
}

// Edited from:
// include/internal/gtest-internal.h#L245
namespace FloatingPoint {
// include/gtest/internal/gtest-port.h#L2198
template <std::size_t size>
class TypeWithSize {
public:
	// This prevents the user from using TypeWithSize<N> with incorrect
	// values of N.
	using UInt = void;
};

// The specialization for size 4.
template <>
class TypeWithSize<4> {
public:
	using Int = std::int32_t;
	using UInt = std::uint32_t;
};

// The specialization for size 8.
template <>
class TypeWithSize<8> {
public:
	using Int = std::int64_t;
	using UInt = std::uint64_t;
};

template <typename T>
using bits_t = typename TypeWithSize<sizeof(T)>::UInt;

// Converts an integer from the sign-and-magnitude representation to
// the biased representation.  More precisely, let N be 2 to the
// power of (kBitCount - 1), an integer x is represented by the
// unsigned number x + N.
//
// For instance,
//
//   -N + 1 (the most negative number representable using
//          sign-and-magnitude) is represented by 1;
//   0      is represented by N; and
//   N - 1  (the biggest number representable using
//          sign-and-magnitude) is represented by 2N - 1.
//
// Read https://en.wikipedia.org/wiki/Signed_number_representations
// for more details on signed number representations.
template <std::floating_point T>
static constexpr bits_t<T> SignAndMagnitudeToBiased(bits_t<T> sam) {
	static constexpr std::size_t kBitCount = 8 * sizeof(T);
	// The mask for the sign bit.
	static constexpr bits_t<T> kSignBitMask = static_cast<bits_t<T>>(1) << (kBitCount - 1);
	if (kSignBitMask & sam) {
		// sam represents a negative number.
		return ~sam + 1;
	} else {
		// sam represents a positive number.
		return kSignBitMask | sam;
	}
}

// Given two numbers in the sign-and-magnitude representation,
// returns the distance between them as an unsigned number.
template <std::floating_point T>
static constexpr bits_t<T> DistanceBetweenSignAndMagnitudeNumbers(bits_t<T> sam1, bits_t<T> sam2) {
	const bits_t<T> biased1 = SignAndMagnitudeToBiased<T>(sam1);
	const bits_t<T> biased2 = SignAndMagnitudeToBiased<T>(sam2);
	return (biased1 >= biased2) ? (biased1 - biased2) : (biased2 - biased1);
}

// Returns true if and only if this number is at most kMaxUlps ULP's away
// from rhs.  In particular, this function:
//
//   - returns false if either number is (or both are) NAN.
//   - treats really large numbers as almost equal to infinity.
//   - thinks +0.0 and -0.0 are 0 ULP's apart.
template <std::floating_point T>
constexpr bool AlmostEquals(T lhs, T rhs) {
	// How many ULP's (Units in the Last Place) we want to tolerate when
	// comparing two numbers.  The larger the value, the more error we
	// allow.  A 0 value means that two numbers must be exactly the same
	// to be considered equal.
	//
	// The maximum error of a single floating-point operation is 0.5
	// units in the last place.  On Intel CPU's, all floating-point
	// calculations are done with 80-bit precision, while double has 64
	// bits.  Therefore, 4 should be enough for ordinary use.
	//
	// See the following article for more details on ULP:
	// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
	static constexpr std::uint32_t kMaxUlps = 4;

	// The IEEE standard says that any comparison operation involving
	// a NAN must return false.
	if (std::isnan(lhs) || std::isnan(rhs)) return false;

	return DistanceBetweenSignAndMagnitudeNumbers<T>(std::bit_cast<bits_t<T>>(lhs), std::bit_cast<bits_t<T>>(rhs)) <=
	       kMaxUlps;
}
};  // namespace FloatingPoint

// include/gtest/gtest.h#L1582
template <typename RawType>
constexpr bool CmpHelperFloatingPointEQ(RawType lhs_value, RawType rhs_value) {
	return FloatingPoint::AlmostEquals(lhs_value, rhs_value);
}

// src/gtest.cc#L1685
constexpr bool DoubleNearPredFormat(double val1, double val2, double abs_error) {
	// std::{fabs, isinf, signbit} are constexpr since C++23 [P0533].
	// We cannot check against __cpp_lib_constexpr_cmath >= 202202L because
	// support is not complete as of today (2025-10-26) and no one sets the
	// FTM. see: https://en.cppreference.com/w/cpp/compiler_support/23.html
	// for P0533R9.

	// We want to return success when the two values are infinity and at
	// least one of the following is true:
	//  * The values are the same-signed infinity.
	//  * The error limit itself is infinity.
	// This is done here so that we don't end up with a NaN when calculating
	// the difference in values.
	if (std::isinf(val1) && std::isinf(val2) &&
	    (std::signbit(val1) == std::signbit(val2) || (abs_error > 0.0 && std::isinf(abs_error)))) {
		return true;
	}

	const double diff = std::fabs(val1 - val2);
	// All branches from here on evaluate to an assertion failure, shortcut:
	return diff <= abs_error;
}
}  // namespace gtest_constexpr
