#include <gtest/gtest.h>
#include <spdlog/spdlog.h>

int main(int argc, char** argv) {
	spdlog::set_level(spdlog::level::trace);
	SPDLOG_TRACE("Initializing tests");
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
