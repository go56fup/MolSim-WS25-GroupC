#!/bin/sh

do_benchmark() {
	implementation="$1"
	type="$2"
	rm -rf /tmp/assignment3_benchmarking
	benchmark_type="$implementation/$type"
	for file in input/assignment3/task2/benchmark/$benchmark_type/*; do
		echo "Performing $benchmark_type/$(basename $file)"
		hyperfine -w 3 "build/release/MolSim $file -o /tmp/assignment3_benchmarking"
	done
}

cmake --workflow release
do_benchmark linked_cell constant_domain
do_benchmark linked_cell constant_squaresize
do_benchmark linked_cell fitted_domain
