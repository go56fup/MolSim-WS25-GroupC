set(EMBED_TEST_DIR "${CMAKE_BINARY_DIR}/embed_test")
file(MAKE_DIRECTORY "${EMBED_TEST_DIR}")

# Write payload file to be embedded
file(WRITE "${EMBED_TEST_DIR}/payload.bin" "ABC")

# Write test source
file(WRITE "${EMBED_TEST_DIR}/test_embed.cpp" "
static constexpr unsigned char data[] = {
#embed \"payload.bin\"
};

int main() {
    return data[0];
}
")

# Required so try_compile doesn't try to link an executable
# set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

try_compile(
    HAS_EMBED
    "${EMBED_TEST_DIR}/build"
    "${EMBED_TEST_DIR}/test_embed.cpp"
    CMAKE_FLAGS "-DCMAKE_CXX_STANDARD=23"
)

if(HAS_EMBED)
	message(STATUS "#embed support detected")
else()
	message(WARNING "#embed not supported - some tests will not be run")
endif()

target_compile_definitions(project_options INTERFACE HAS_EMBED=$<IF:$<BOOL:${HAS_EMBED}>,1,0>)
