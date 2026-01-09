set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR x86_64)

set(
	CMAKE_CXX_FLAGS_INIT
	# knl is deprecated in GCC 14.
	"-Wno-error=deprecated -march=knl -mmmx -mno-3dnow -msse -msse2 -msse3 -mssse3 -mno-sse4a -mcx16 -msahf -mmovbe -maes -msha -mpclmul -mpopcnt -mabm -mno-lwp -mfma -mno-fma4 -mno-xop -mbmi -msgx -mbmi2 -mno-tbm -mavx -mavx2 -msse4.2 -msse4.1 -mlzcnt -mrtm -mhle -mrdrnd -mf16c -mfsgsbase -mrdseed -mprfchw -madx -mfxsr -mxsave -mxsaveopt -mavx512f -mno-avx512er -mavx512cd -mno-avx512pf -mno-prefetchwt1 -mclflushopt -mxsavec -mxsaves -mavx512dq -mavx512bw -mavx512vl -mavx512ifma -mavx512vbmi -mno-avx5124fmaps -mno-avx5124vnniw -mclwb -mno-mwaitx -mno-clzero -mpku -mrdpid --param l1-cache-size=48 --param l1-cache-line-size=64 --param l2-cache-size=107520 -mtune=generic"
)  
