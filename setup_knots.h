#pragma once
#include "filter.h"

void calculate_c(KnotsType &knots, double c0, const double nx, vdouble2 domain_min, vdouble2 domain_max);
void setup_knots(KnotsType &knots, ParticlesType &fibres, const double fibre_radius, const double ibre_resolution_factor, const double nx, vdouble2 domain_min, vdouble2 domain_max, const double c0, const double k, const bool periodic, const int nbucket);
