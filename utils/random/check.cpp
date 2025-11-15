#include <array>
#include <bit>
#include <iostream>

#include "Constants.h"

inline constexpr auto numbers = [] consteval {
	const char data[random_table_size * sizeof(double)] {
#embed "random.bin"
	};
	return std::bit_cast<std::array<double, random_table_size>>(data);
}();

static_assert(numbers.size() == random_table_size);

int main() {
	for (std::size_t i = 0; i < random_table_size; ++i) {
		std::cout << i << ' ' << numbers[i] << '\n';
	}
}
