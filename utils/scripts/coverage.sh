#!/bin/sh

BUILD="./build/debug-gcc"

ctest --test-dir "$BUILD" -T Coverage
rm -rf "$BUILD/coverage"
mkdir "$BUILD/coverage"
lcov --output-file "$BUILD/coverage/coverage.info" --capture --ignore-errors inconsistent \
	 --directory "$BUILD" --demangle-cpp --base-directory "." --no-external --exclude "*_deps*" \
     --exclude "FixedString.h"
genhtml -o "$BUILD/coverage" "$BUILD/coverage/coverage.info"
xdg-open "$BUILD/coverage/index.html"
