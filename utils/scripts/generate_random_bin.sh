#!/bin/sh

set -e
dir="$(mktemp -d)"
project="$(pwd)"

compile() {
	clang++ -std=c++23 -Wall -Wextra -Wpedantic -Werror -Wno-c23-extensions -I$project/include $1 -o $2
}

cd "$dir"
compile "$project/utils/random/generate.cpp" "$dir/gen"
"$dir/gen" > "$dir/gen_human"
compile "$project/utils/random/check.cpp" "$dir/check"
"$dir/check" > "$dir/check_human"
diff "$dir/gen_human" "$dir/check_human"
mv "$dir/random.bin" "$project/include/utils"
cd "$project"
rm -rf "$dir"
