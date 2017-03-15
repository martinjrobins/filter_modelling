#pragma once

#include "filter.h"
#include "particular_MAPS.h"

void solve_stokes_LMAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0);
