#include <fstream>
#include <iostream>
#include <random>
#include <span>

#include "Constants.h"

int main() {
	std::ofstream file("random.bin", std::ios::binary);
	static constexpr auto seed = 42;
	static std::default_random_engine randomEngine(seed);
	std::normal_distribution<double> normalDistribution{0, 1};
	for (std::size_t i = 0; i < random_table_size; ++i) {
		const double num = normalDistribution(randomEngine);
		std::cout << i << ' ' << num << '\n';
		const std::span num_bytes(reinterpret_cast<const char*>(&num), sizeof(num));
		file.write(num_bytes.data(), num_bytes.size());
	}
}
