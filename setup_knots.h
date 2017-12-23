#pragma once
#include "filter.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#include <iostream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;


void calculate_c(KnotsType &knots, double c0, const double nx, vdouble2 domain_min, vdouble2 domain_max);
CDT setup_knots(KnotsType &knots, ParticlesType &fibres, const double fibre_radius, const double ibre_resolution_factor, const double nx, vdouble2 domain_min, vdouble2 domain_max, const double c0, const double k, const bool periodic, const int nbucket);
