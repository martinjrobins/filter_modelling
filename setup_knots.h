#pragma once
#include "filter.h"

void calculate_c(KnotsType &knots, double c0, const double nx, double2 domain_min, double2 domain_max);
void setup_knots(KnotsType &knots, ParticlesType &fibres, const double fibre_radius, const double fibre_resolution_factor, const double nx, double2 domain_min, double2 domain_max, const double c0, const double k, const double epsilon_strength, const double epsilon_falloff, const bool periodic);
