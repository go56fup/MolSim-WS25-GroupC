function(get_third_party)
	set(oneValueArgs TARGET)
	cmake_parse_arguments(
		PARSE_ARGV
		0
		"THIRD"
		""
		"${oneValueArgs}"
		""
	)
	add_library("${THIRD_TARGET}" INTERFACE)
	include(FetchContent)

	FetchContent_Declare(
		spdlog
		GIT_REPOSITORY https://github.com/gabime/spdlog.git
		GIT_TAG v1.16.0
	)
	FetchContent_MakeAvailable(spdlog)

	FetchContent_Declare(
		fmt
		GIT_REPOSITORY https://github.com/fmtlib/fmt
		GIT_TAG 12.1.0
	)
	FetchContent_MakeAvailable(fmt)
	target_link_libraries("${THIRD_TARGET}" INTERFACE fmt::fmt)

	if(SPDLOG_BUILD_EXAMPLE_HO)
		target_link_libraries("${THIRD_TARGET}" INTERFACE spdlog::spdlog_header_only)
	else()
		target_link_libraries("${THIRD_TARGET}" INTERFACE spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)
	endif()
	target_include_directories("${THIRD_TARGET}" SYSTEM INTERFACE "${spdlog_SOURCE_DIR}/include")
endfunction()
