#include <array>
#include <bit>
#include <iostream>
#include <span>

#include "definitions.h"

inline constexpr auto numbers = [] consteval {
	const char data[element_count * sizeof(double)] {
#embed "random.bin"
	};
	return std::bit_cast<std::array<double, element_count>>(data);
}();

static_assert(numbers.size() == element_count);

int main() {
	for (int i = 0; i < element_count; ++i) {
		std::cout << i << ' ' << numbers[i] << '\n';
	}
}
