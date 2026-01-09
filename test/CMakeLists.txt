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
include(GoogleTest)

option(BUILD_STATIC_TESTING "Enable compile-time tests" OFF)

function(create_test_target)
	set(oneValueArgs TARGET PREFIX)
	set(multiValueArgs SOURCES)
	cmake_parse_arguments(
		PARSE_ARGV
		0
		"tests"
		""
		"${oneValueArgs}"
		"${multiValueArgs}"
	)
	add_executable("${tests_TARGET}" ${tests_SOURCES})
	target_link_libraries("${tests_TARGET}" PRIVATE project_options)
	target_link_libraries("${tests_TARGET}" PRIVATE GTest::gtest)

	# Disable warnings from inside GTest.
	target_include_directories("${tests_TARGET}" SYSTEM PRIVATE "${googletest_SOURCE_DIR}/include")

	target_include_directories(
		"${tests_TARGET}" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include"
	)
	
	if(CXX_COMPILER_ID STREQUAL Clang)
		target_compile_options("${tests_TARGET}" PRIVATE -fmacro-backtrace-limit=0)
	endif()
	gtest_discover_tests("${tests_TARGET}" TEST_PREFIX "${tests_PREFIX}")
endfunction()

file(GLOB_RECURSE TEST_SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

if(ENABLE_STATIC_TESTING)
	create_test_target(
		TARGET
		MolSim_Test_static
		PREFIX
		"compile."
		SOURCES
		${TEST_SRC}
	)
	target_compile_definitions(MolSim_Test_static PRIVATE "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_OFF")
	if(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
		target_compile_options(
			MolSim_Test_static PRIVATE "-fconstexpr-loop-limit=2147483647"
										"-fconstexpr-ops-limit=2147483647"
		)
	elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
		target_compile_options(MolSim_Test_static PRIVATE "-fconstexpr-steps=4294967295")
	else()
		message(WARNING "Cannot detect compiler of id ${CMAKE_CXX_COMPILER_ID} and \
set constexpr evaluation steps to maximum. Compilation may fail."
		)
	endif()
endif()
create_test_target(
	TARGET
	MolSim_Test_runtime
	PREFIX
	"runtime."
	SOURCES
	${TEST_SRC}
)
target_compile_definitions(MolSim_Test_runtime PRIVATE "GTEST_CONFIG_RUNTIME_STATIC_EXPECT")
target_compile_definitions(MolSim_Test_runtime PRIVATE "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE")
if(ENABLE_COVERAGE)
	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		target_compile_options(
			MolSim_Test_runtime PRIVATE "--coverage" "-fno-inline" "-fno-inline-small-functions"
										"-fno-default-inline"
		)
	elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
		target_compile_options(MolSim_Test_runtime PRIVATE "--coverage" "-fno-inline-functions")
	endif()
	target_link_options(MolSim_Test_runtime PRIVATE "--coverage")
endif()

include(doxygen)

set(_1 "\\1")
set(_2 "\\2")

# Used to annotate assertion macros directly, e.g. Compile-time equivalent of `ASSERT_EQ`.
doxygen_generate_variadic(
	NAME
	gtest_equ
	NORMAL_ARGC
	1
	VARIADIC_UP_TO
	6
	PATTERN
	"@def STATIC_${_1}(**)^^@brief Compile-time equivalent of `${_1}`."
)
define_param_shorthands(gtest_equ 1)

doxygen_generate_variadic(
	NAME
	gtest_equ_gcc
	NORMAL_ARGC
	1
	VARIADIC_UP_TO
	6
	PATTERN
	"@def GCC_STATIC_${_1}(**)^^@brief Compile-time equivalent of `${_1}` (only on GCC)."
)
define_param_shorthands(gtest_equ_gcc 1)

# Begins a section for given macro suffix, e.g. STATIC_{ASSERT,EXPECT}_EQ(val1, val2)
doxygen_generate_variadic(
	NAME
	gtest_header
	NORMAL_ARGC
	1
	VARIADIC_UP_TO
	3
	PATTERN
	"@name STATIC_\\{ASSERT\\,EXPECT\\}_${_1}(**)^^@\\{"
	SEPARATOR
	"\\, "
)
define_param_shorthands(gtest_header 1)

# Annotation for @param's, e.g. @param val1 Left-hand (...) @param val2 Right-hand (...)
doxygen_alias(
	gtest_params
	2
	"@param ${_1} Left-hand side of the comparison.^^@param ${_2} Right-hand side of the comparison."
)
define_param_shorthands(gtest_params 0)

# Section for macros of comparison operators, e.g. STATIC_{ASSERT,EXPECT}_GT(val1, val2) Verifies
# that `val1 > val2`.
doxygen_alias(
	gtest_section_cmp 2 "@gtest_header_val{${_1}}^^@brief Verifies that `val1` ${_2} `val2`."
)

define_whole_section_shorthand(gtest_section_cmp val)

# Section for macros for comparing strings, e.g. STATIC_{ASSERT,EXPECT}_STREQ(str1, str2) Verifies
# that two C strings `str1` and `str2` have the same contents.
doxygen_alias(
	gtest_section_str 2
	"@gtest_header_str{${_1}}^^@brief Verifies that two C strings `str1` and `str2` have ${_2}."
)
define_whole_section_shorthand(gtest_section_str str)

doxygen_alias(same 0 "the same contents")
doxygen_alias(different 0 "different contents")
doxygen_alias(case_insens 0 "\\, ignoring case")
# @gtest_section_str{STREQ,@same} @gtest_section_str{STRCASNE,@different@case_insens}

# Section for macros comparing floating point values, e.g. STATIC_{ASSERT,EXPECT}_FLOAT_EQ(val1,
# val2) Verifies that two `float` values `val1` and `val2` are approximately equal.
doxygen_alias(
	gtest_section_fp
	2
	"@gtest_header_val{${_1}}^^@brief Verifies that two `${_2}` values `val1` and `val2` are approximately equal."
)
define_whole_section_shorthand(gtest_section_fp val)

# Used to annotate PARAMN @param's, e.g. @param val3 Third argument to pass into `pred` (for
# `*PRED3` or higher).
doxygen_alias(
	gtest_nth_param 2
	"@param val${_1} ${_2} argument to pass into `pred` (for `*PRED${_1}` or higher)."
)
