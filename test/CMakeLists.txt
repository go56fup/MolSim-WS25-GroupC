if(NOT BUILD_TESTING)
	return()
endif()

include(FetchContent)
FetchContent_Declare(
	googletest
	GIT_REPOSITORY https://github.com/google/googletest
	BRANCH main
	GIT_TAG v1.17.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(
	gtest_force_shared_crt
	ON
	CACHE BOOL "" FORCE
)
FetchContent_MakeAvailable(googletest)

add_executable(MolSim_Test Test.cpp)
target_link_libraries(MolSim_Test PRIVATE GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(MolSim_Test)
