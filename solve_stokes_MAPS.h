#pragma once
#include "filter.h"
#include "particular_MAPS.h"

double solve_stokes_MAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0);
