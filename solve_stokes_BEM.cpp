#include "solve_stokes_BEM.h"

void setup_elements(ElementsType &elements, ElementsType &boundarye,
                    ParticlesType &fibres, vdouble2 domain_min,
                    vdouble2 domain_max, const unsigned int nx,
                    const double fibre_radius) {
  const double L = domain_max[0] - domain_min[0];
  const double Ly = domain_max[1] - domain_min[1];
  const double dtheta = 2 * PI / nx;

  elements.resize(fibres.size() * nx);
  vdouble2 ns_buffer(0, L / 10.0);

  std::cout << "setup elements: domain_min = " << domain_min
            << " domain_max = " << domain_max
            << " fibre_radius = " << fibre_radius << " nx = " << nx
            << std::endl;

  for (int ii = 0; ii < fibres.size(); ++ii) {
    const vdouble2 origin = get<position>(fibres)[ii];
    const double charge = get<charge>(fibres)[ii];
    bool outside, started;
    started = false;
    for (int kk = 0; kk < nx; ++kk) {
      ElementsType::reference p = elements[ii * nx + kk];
      get<position>(p) =
          origin +
          fibre_radius * vdouble2(std::cos(kk * dtheta), std::sin(kk * dtheta));
      get<point_a>(p) =
          origin + fibre_radius * vdouble2(std::cos((kk - 0.5) * dtheta),
                                           std::sin((kk - 0.5) * dtheta));
      get<point_b>(p) =
          origin + fibre_radius * vdouble2(std::cos((kk + 0.5) * dtheta),
                                           std::sin((kk + 0.5) * dtheta));
      get<boundary>(p) = false;
      get<charge>(p) = charge;
      get<normal>(p) =
          eigen_vector(std::cos(kk * dtheta), std::sin(kk * dtheta));

      for (int i = 0; i < 2; ++i) {
        if (get<position>(p)[i] < domain_min[i]) {
          get<position>(p)[i] += domain_max[i] - domain_min[i];
          get<point_a>(p)[i] += domain_max[i] - domain_min[i];
          get<point_b>(p)[i] += domain_max[i] - domain_min[i];
        } else if ((get<position>(p)[i] >= domain_max[i])) {
          get<position>(p)[i] -= domain_max[i] - domain_min[i];
          get<point_a>(p)[i] -= domain_max[i] - domain_min[i];
          get<point_b>(p)[i] -= domain_max[i] - domain_min[i];
        }
      }
      // std::cout << "element "<<ii*nx+kk<<" has p = "<<get<position>(p)<<" p1
      // = "<<get<point_a>(p)<<" and p2 = "<<get<point_b>(p) << std::endl;
    }
  }

  /*
  const double dx_aim =
  (get<position>(elements)[0]-get<position>(elements)[1]).norm(); const int
  n_inlet = std::ceil(L/dx_aim); const double dx = L/n_inlet;
  //elements.resize(fibres.size()*nx + 2*n_inlet);
  boundarye.resize(2*n_inlet);

  for (int i = 0; i < n_inlet; ++i) {
      //ElementsType::reference p = elements[fibres.size()*nx + 2*i];
      ElementsType::reference p = boundarye[2*i];
      get<position>(p) =
  vdouble2((i+0.5)*dx,domain_min[1]-2*(domain_max[1]-domain_min[1]));
      get<point_a>(p) =
  vdouble2((i)*dx,domain_min[1]-2*(domain_max[1]-domain_min[1]));
      get<point_b>(p) =
  vdouble2((i+1)*dx,domain_min[1]-2*(domain_max[1]-domain_min[1]));
      get<boundary>(p) = true;

      //ElementsType::reference p2 = elements[fibres.size()*nx + 2*i+1];
      ElementsType::reference p2 = boundarye[2*i+1];
      get<position>(p2) =
  vdouble2((i+0.5)*dx,domain_max[1]+2*(domain_max[1]-domain_min[1]));
      get<point_a>(p2) =
  vdouble2((i)*dx,domain_max[1]+2*(domain_max[1]-domain_min[1]));
      get<point_b>(p2) =
  vdouble2((i+1)*dx,domain_max[1]+2*(domain_max[1]-domain_min[1]));
      get<boundary>(p2) = true;

  }
  */

  elements.init_neighbour_search(domain_min, domain_max, vbool2(true, false));
  // boundarye.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,vbool2(true,false));
}

double solve_stokes_BEM(KnotsType &knots, ElementsType &elements,
                        ElementsType &boundarye, const double alpha,
                        const int nlambda, const int nmu) {
  std::cout << "solving stokes with BEM..." << elements.size() << std::endl;

  const double flow_rate = 1.0;
  const double mu = 1.0;

  typedef typename KnotsType::position position;
  typedef typename position::value_type const &const_position_reference;
  typedef typename KnotsType::const_reference const_particle_reference;

  const vdouble2 min = elements.get_min();
  const vdouble2 max = elements.get_max();
  const double h =
      (get<point_b>(elements)[0] - get<point_a>(elements)[0]).norm();

  auto kernel = make_greens_kernel_2d1p(alpha, nlambda, nmu, h, min, max, true);
  auto A = create_dense_operator(elements, elements, kernel);
  /*
  auto Abl = create_dense_operator(elements,boundarye,
          make_boundary_layer_kernel_2d1p(mu,min,max,false));
          */

  std::cout << "setup equations..." << std::endl;

  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Map<vector_type> map_type;

  const size_t N = elements.size();
  // const size_t Nb = boundarye.size();

  vector_type source(2 * N);
  vector_type alphas(2 * N);
  // vector_type vel(Nb);
  matrix_type A_eigen(2 * N, 2 * N);

  std::cout << "assemble matrix..." << std::endl;
  // c[i] = c0;
  A.assemble(A_eigen);

  for (int ii = 0; ii < N; ++ii) {
    source(2 * ii) = 0;
    source(2 * ii + 1) = flow_rate * 4 * PI * mu;
  }

  std::cout << "solve w BEM ..." << std::endl;
  alphas = A_eigen.householderQr().solve(source);
  // alphas = A_eigen.colPivHouseholderQr().solve(source);
  double relative_error = (A_eigen * alphas - source).norm() / source.norm();
  std::cout << "The relative error is:\n" << relative_error << std::endl;

  std::cout << "done solve..." << std::endl;

  for (int ii = 0; ii < N; ++ii) {
    get<traction>(elements)[ii][0] = alphas[2 * ii];
    get<traction>(elements)[ii][1] = alphas[2 * ii + 1];
  }

  vtkWriteGrid("BEMelements", 0, elements.get_grid(true));

  std::cout << "assemble knot matrix..." << std::endl;
  matrix_type A_knot(2 * knots.size(), 2 * N);
  auto kernel2 =
      make_greens_kernel_2d1p(alpha, nlambda, nmu, h, min, max, false);
  auto Aknots = create_dense_operator(knots, elements, kernel2);

  vector_type svel = Aknots * alphas / (4 * PI * mu);

  for (int ii = 0; ii < knots.size(); ++ii) {
    get<velocity>(knots)[ii][0] = svel[2 * ii];
    get<velocity>(knots)[ii][1] = svel[2 * ii + 1] - flow_rate;
  }

  vtkWriteGrid("BEMknots", 0, knots.get_grid(true));

  std::cout << "done solving stokes" << std::endl;
  return relative_error;
}
