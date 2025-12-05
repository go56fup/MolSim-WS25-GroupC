cmake_minimum_required(VERSION 3.31)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
	MolSim_GroupC
	VERSION 0.0.1
	LANGUAGES CXX
)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
	set(
		CMAKE_BUILD_TYPE
		"Release"
		CACHE STRING
				"Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
				FORCE
	)
endif(NOT CMAKE_BUILD_TYPE)

# collect all cpp files
file(
	GLOB_RECURSE
	MY_SRC
	"${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
	# header don't need to be included but this might be necessary for some IDEs
	"${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)

# create make target
add_executable(MolSim ${MY_SRC})
add_library(project_options INTERFACE)

target_include_directories(
	project_options INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/input
)

target_compile_options(
	project_options
	INTERFACE $<$<CXX_COMPILER_ID:GNU>:-Wall
				-Wextra
				-Wpedantic
				-Wconversion
				# We use #embed, which is technically C++26, but we allow C23.
				-Wno-c++26-extensions
				>
				# disable some spam warnings for icpc...
				$<$<CXX_COMPILER_ID:Intel>:-w3
				-wd383,981,1418,1572,2259>
				$<$<CXX_COMPILER_ID:Clang>:-Wall
				-Wextra
				-Wpedantic
				-Wconversion
				-Wno-c23-extensions
				>
)

if(ENABLE_SANITIZERS)
	target_compile_options(project_options INTERFACE "-fsanitize=undefined,address")
	target_link_options(project_options INTERFACE "-fsanitize=undefined,address")
endif()

if(ENABLE_VTK_OUTPUT)
	target_compile_definitions(project_options INTERFACE "ENABLE_VTK_OUTPUT")
endif()

if(BENCHMARK_PRESET)
	target_compile_definitions(MolSim PRIVATE "BENCHMARK_PRESET=${BENCHMARK_PRESET}")
endif()

target_compile_definitions(MolSim PRIVATE "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_${LOG_LEVEL}")

# add additional modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(vtk)
include(clang-format)
include(doxygen)
include(clang-tidy)
include(third-party)

get_third_party(TARGET third_party)
target_link_libraries(project_options INTERFACE third_party)

include(CTest)
add_subdirectory(test)
generate_documentation_target()

target_link_libraries(MolSim PRIVATE project_options)
