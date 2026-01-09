#pragma once

#ifdef IS_TESTING
#define STATIC_IF_NOT_TESTING
#else
#define STATIC_IF_NOT_TESTING static
#endif
