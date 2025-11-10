#pragma once

#if ! defined(__clang__) && defined(__GNUC__)
#define IS_GCC
#endif

#ifdef IS_GCC
#define CONSTEXPR_IF_GCC constexpr
#else
#define CONSTEXPR_IF_GCC
#endif
