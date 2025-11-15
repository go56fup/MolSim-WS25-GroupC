#!/bin/sh

set -e

compile() {
	clang++ -std=c++23 -Wall -Wextra -Wpedantic -Werror -Wno-c23-extensions $1 -o $2
}

dir="$(mktemp -d)"
project="$(pwd)"
cd "$dir"
compile "$project/utils/random/generate.cpp" "$dir/gen"
"$dir/gen" > "$dir/gen_human"
compile "$project/utils/random/check.cpp" "$dir/check"
"$dir/check" > "$dir/check_human"
diff "$dir/gen_human" "$dir/check_human"
mv "$dir/random.bin" "$project/include"
cd "$project"
rm -rf "$dir"
