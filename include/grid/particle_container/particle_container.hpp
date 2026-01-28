#pragma once

#include "grid/particle_container/fwd.hpp"
#include "grid/particle_container/impl.hpp"

static_assert(std::ranges::range<
			  decltype((std::declval<particle_container>()).directional_interactions())>);
