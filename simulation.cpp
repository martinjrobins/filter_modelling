#include "filter.h"
#include "solvers.h"
//#include "setup_knots.h"
#include <Eigen/Dense>
#include <boost/archive/xml_oarchive.hpp>
#include <fstream>

int main(int argc, char **argv) {

  unsigned int nout, max_iter_linear, restart_linear, nx;
  int fibre_number, seed, fibre_arrangement;

  double fibre_resolution, fibre_radius, particle_rate, particle_radius,
      react_rate, D, L;
  double dt_aim, k, gamma, rf, c0, particles_charge, fibres_charge,
      fibres_charge_stddev;
  unsigned int solver_in;
  bool periodic;
  int electrostatics_fibre;
  std::string base_dir;

  po::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")(
      "max_iter_linear",
      po::value<unsigned int>(&max_iter_linear)->default_value(2000),
      "maximum iterations for linear solve")(
      "restart_linear",
      po::value<unsigned int>(&restart_linear)->default_value(2001),
      "iterations until restart for linear solve")(
      "linear_solver", po::value<unsigned int>(&solver_in)->default_value(2),
      "linear solver")("nout",
                       po::value<unsigned int>(&nout)->default_value(100),
                       "number of output points")(
      "k", po::value<double>(&k)->default_value(1.00), "spring constant")(
      "D", po::value<double>(&D)->default_value(0.01), "diffusion constant")(
      "particle_rate", po::value<double>(&particle_rate)->default_value(500.0),
      "particle rate")("periodic",
                       po::value<bool>(&periodic)->default_value(true),
                       "periodic in x")(
      "react_rate", po::value<double>(&react_rate)->default_value(1.0),
      "particle reaction rate")

      ("c0", po::value<double>(&c0)->default_value(0.0835),
       "kernel constant")("nx", po::value<unsigned int>(&nx)->default_value(10),
                          "number of knots around each fibre")(
          "fibre_arrangement",
          po::value<int>(&fibre_arrangement)->default_value(0),
          "(0=regular, 1=hexigonal 2=random)")(
          "fibre_radius", po::value<double>(&fibre_radius)->default_value(0.3),
          "radius of fibres")(
          "particle_radius",
          po::value<double>(&particle_radius)->default_value(0.3 / 100.0),
          "radius of fibres")(
          "electrostatics_fibre",
          po::value<int>(&electrostatics_fibre)->default_value(0),
          "particles are electrostaticly attracted to the fibres")(
          "particles_charge",
          po::value<double>(&particles_charge)->default_value(1.0),
          "charge on particles")(
          "fibres_charge",
          po::value<double>(&fibres_charge)->default_value(1.0),
          "fibres_charge_stddev",
          po::value<double>(&fibres_charge_stddev)->default_value(1.0),

          "charge on fibres")("seed", po::value<int>(&seed)->default_value(10),
                              "seed")(
          "base_dir",
          po::value<std::string>(&base_dir)->default_value("default"),
          "base filename for output")(
          "fibre_number", po::value<int>(&fibre_number)->default_value(5),
          "number of fibres")("domain_size",
                              po::value<double>(&L)->default_value(5),
                              "domain size")(
          "dt", po::value<double>(&dt_aim)->default_value(0.001), "timestep");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  KnotsType knots;
  ElementsType elements;
  ElementsType boundary;
  ParticlesType particles;
  ParticlesType dead_particles;
  ParticlesType fibres;

  const int electrostatics_sum = 5;
  const double c2 = std::pow(c0, 2);
  const double c3 = std::pow(c0, 3);
  const double c4 = std::pow(c0, 4);
  const bool reflective = false;
  const double mu = 1.0;
  const double flow_rate = 1.0;
  const double Tf = 2.0 * L;
  const int timesteps = Tf / dt_aim;
  const double fibre_radius2 = std::pow(fibre_radius, 2);
  const double dt = Tf / timesteps;
  const double dt_adapt = (1.0 / 100.0) * PI / sqrt(2 * k);
  const vdouble2 domain_min(0, -1);
  const vdouble2 domain_max(L, L + 1);
  const vdouble2 ns_buffer_fibres(L / 3, 2 * particle_radius);
  const vdouble2 ns_buffer_particles(static_cast<int>(reflective) * L / 3,
                                     L / 2);

  std::default_random_engine rnd_generator(seed);

  // fibres
  {

    std::uniform_real_distribution<double> xrange(domain_min[0], domain_max[0]);
    std::uniform_real_distribution<double> yrange(
        domain_min[1] + 1 + fibre_radius, domain_max[1] - 1 - fibre_radius);
    std::uniform_real_distribution<double> U(0, 1);
    std::normal_distribution<double> normal(0, 1);
    typename ParticlesType::value_type p;
    if (fibre_arrangement == 2) {
      for (int jj = 0; jj < fibre_number; ++jj) {
        bool free_position = false;
        while (free_position == false) {
          get<position>(p) =
              vdouble2(xrange(rnd_generator), yrange(rnd_generator));
          free_position = true;
          /*
             for (auto tpl: euclidean_search(fibres.get_query(),
             get<position>(p),
             2*fibre_radius+2*particle_radius)) {
             */
          for (auto f : fibres) {
            if ((get<position>(p) - get<position>(f)).norm() <
                2 * fibre_radius + 2 * particle_radius) {
              free_position = false;
              break;
            }
          }
        }
        fibres.push_back(p);
      }
    } else if (fibre_arrangement == 1) {
      for (int jj = 0; jj < fibre_number; ++jj) {
        for (int ii = 0; ii < ((jj % 2 == 1) ? fibre_number + 1 : fibre_number);
             ++ii) {
          const double dx = L / fibre_number;
          vdouble2 origin =
              vdouble2(((jj % 2 == 1) ? ii : ii + 0.5) * dx, (jj + 0.5) * dx);
          get<position>(p) = origin;
          fibres.push_back(p);
        }
      }
    } else {
      for (int ii = 0; ii < fibre_number; ++ii) {
        for (int jj = 0; jj < fibre_number; ++jj) {
          const double dx = L / fibre_number;
          const vdouble2 origin = vdouble2((ii + 0.5) * dx, (jj + 0.5) * dx);
          get<position>(p) = origin;
          fibres.push_back(p);
        }
      }
    }
    for (size_t i = 0; i < fibres.size(); ++i) {
      if (electrostatics_fibre == 1) {
        get<charge>(fibres)[i] =
            U(rnd_generator) < 0.5 ? -fibres_charge : fibres_charge;
      } else if (electrostatics_fibre == 2) {
        get<charge>(fibres)[i] = fibres_charge * normal(rnd_generator);
      }
    }
    if (electrostatics_fibre == 2) {
      const double variance = std::pow(fibres_charge_stddev, 2);
      Eigen::MatrixXd C(fibres.size(), fibres.size());
      for (size_t i = 0; i < fibres.size(); ++i) {
        for (size_t j = 0; j < fibres.size(); ++j) {
          C(i, j) =
              std::exp(-(get<position>(fibres)[i] - get<position>(fibres)[j])
                            .squaredNorm() /
                       (2 * variance));
        }
      }
      Eigen::LLT<Eigen::MatrixXd> lltOfC(C);
      Eigen::Map<Eigen::VectorXd> ch(get<charge>(fibres).data(), fibres.size());
      ch = lltOfC.matrixL() * ch;
      for (size_t i = 0; i < fibres.size(); ++i) {
        const double y = get<position>(fibres)[i][1];
        get<charge>(fibres)[i] *=
            0.5 * (std::tanh(y - 1.5) + std::tanh(L - 1.5 - y));
      }
    }

    particles.init_neighbour_search(domain_min,
                                    domain_max + ns_buffer_particles,
                                    vbool2(!reflective, false));

    fibres.init_neighbour_search(domain_min - ns_buffer_fibres,
                                 domain_max + ns_buffer_fibres, vbool2(false));

    std::cout << "added " << fibres.size() << " fibres" << std::endl;

    // vtkWriteGrid("init_knots_fibres",0,fibres.get_grid(true));
    {
      std::ofstream ofs("init_knots_fibres");
      boost::archive::xml_oarchive oa(ofs);
      oa << BOOST_SERIALIZATION_NVP(fibres);
    }
  }

  //
  // SETUP KNOTS
  //
  // setup_knots(knots, fibres, fibre_radius, fibre_resolution, nx, domain_min,
  // domain_max, k,false,10);

  //
  // SETUP ELEMENTS
  //
  setup_elements(elements, boundary, fibres, domain_min, domain_max, nx,
                 fibre_radius);

  max_iter_linear = knots.size() * 4;
  restart_linear = max_iter_linear + 1;

  //
  // SOLVE STOKES
  //
  const double alpha =
      (elements.get_max() - elements.get_min()).prod() / (4.0 * PI);
  const double h =
      (get<point_b>(elements)[0] - get<point_a>(elements)[0]).norm();
  const int nlambda = 2;
  const int nmu = 2;
  const double relative_error =
      solve_stokes_BEM(knots, elements, boundary, alpha, nlambda, nmu);
  solve_laplace_BEM(knots, elements, fibres_charge);

  auto stokesSLPkernel = make_greens_kernel_2d1p(alpha, nlambda, nmu, h,
                                                 domain_min, domain_max, false);
  auto stokesSLP = create_dense_operator(particles, elements, stokesSLPkernel);
  auto laplaceGradSLPkernel =
      make_laplace_gradSLP_2d1p(domain_min, domain_max, false);
  auto laplaceGradSLP =
      create_dense_operator(particles, elements, laplaceGradSLPkernel);

  std::cout << "starting timesteps!" << std::endl;
  const int timesteps_per_out = timesteps / nout;
  auto t0 = Clock::now();
  auto t1 = Clock::now();
  double time_vel_eval = 0;
  double time_vel_rest = 0;
  for (int ii = 0; ii < timesteps; ++ii) {
    // add new particles
    std::poisson_distribution<int> poisson(particle_rate * dt);
    std::uniform_real_distribution<double> uniform(2 * particle_radius,
                                                   L - 2 * particle_radius);
    const int new_n = poisson(rnd_generator);
    for (int jj = 0; jj < new_n; ++jj) {
      ParticlesType::value_type p;
      get<position>(p) =
          vdouble2(uniform(rnd_generator), domain_max[1] - particle_radius);
      particles.push_back(p);
    }

    t0 = Clock::now();

#pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i) {
      get<stokes_velocity>(particles)[i] = eigen_vector::Zero();
      get<electro_velocity>(particles)[i] = eigen_vector::Zero();
    }

    // evaluate velocity field
    stokesSLP.get_first_kernel().evaluate(get<stokes_velocity>(particles),
                                          get<traction>(elements));
    if (electrostatics_fibre) {
      laplaceGradSLP.get_first_kernel().evaluate(
          get<electro_velocity>(particles), get<gradP>(elements));
    }

#pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i) {
      get<stokes_velocity>(particles)[i] *= 1.0 / (4.0 * PI * mu);
      get<stokes_velocity>(particles)[i][1] -= flow_rate;
      get<electro_velocity>(particles)[i] *= -1;
      get<velocity>(particles)[i] = get<stokes_velocity>(particles)[i] +
                                    get<electro_velocity>(particles)[i];
    }

    t1 = Clock::now();
    time_vel_eval += (t1 - t0).count();
    t0 = Clock::now();

/*
   if (electrostatics_particles) {
   make fmm here, sum into velocity
   }
   */
#pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i) {
      ParticlesType::reference p = particles[i];
      /*
      if (electrostatics_fibre) {
          for (int i = 0; i < fibres.size(); ++i) {
              ParticlesType::reference f = fibres[i];
              for (int i = -electrostatics_sum+1; i < electrostatics_sum ; ++i)
      { const eigen_vector dx( get<position>(f)[0]-get<position>(p)[0]+i*L,
                          get<position>(f)[1]-get<position>(p)[1]);
                  const double scalar =
      fibres_charge/std::pow(dx.squaredNorm()+fibre_radius2,1.5);
                  get<velocity>(p) += dx*scalar;
              }
          }
      }
      */
      std::normal_distribution<double> normal;
      const vdouble2 N(normal(get<generator>(particles)[i]),
                       normal(get<generator>(particles)[i]));

      get<position>(p) += std::sqrt(2.0 * D * dt) * N;
      for (int i = 0; i < 2; ++i) {
        get<position>(p)[i] += dt * get<velocity>(p)[i];
      }
      /*
      if (get<position>(p)[1] > domain_max[1]+ns_buffer_particles[1] ||
      get<position>(p)[1] < domain_min[1]-ns_buffer_particles[1]) { std::cout <<
      "error: particle at "<<get<position>(p)<<" with velocity =
      "<<get<velocity>(p) << std::endl;
      }
      */
    }

    if (ii % timesteps_per_out == 0) {
      std::cout << "timestep " << ii << " of " << timesteps
                << " (time_vel_eval = " << time_vel_eval
                << " time_vel_rest = " << time_vel_rest << ") "
                << particles.size() << " particles alive" << std::endl;
      /*
      {
          std::ostringstream ostr;
          ostr << std::setfill('0') << std::setw(5) << ii;
          std::ofstream ofs(
                  base_dir + "simulation" + ostr.str() + ".xml");
          boost::archive::xml_oarchive oa(ofs);
          oa << BOOST_SERIALIZATION_NVP(particles);
          oa << BOOST_SERIALIZATION_NVP(fibres);
          oa << BOOST_SERIALIZATION_NVP(dead_particles);
      }
      */
      vtkWriteGrid((base_dir + "particles").c_str(), ii,
                   particles.get_grid(true));
      vtkWriteGrid((base_dir + "fibres").c_str(), ii, fibres.get_grid(true));
      vtkWriteGrid((base_dir + "dead_particles").c_str(), ii,
                   dead_particles.get_grid(true));
    }

    // react with fibres
    // alive_[a] = !any(bf,U[a]<react_rate);
    std::uniform_real_distribution<double> uni(0, 1);
    //#pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i) {
      ParticlesType::reference p = particles[i];
      for (const auto &i :
           euclidean_search(fibres.get_query(), get<position>(p),
                            fibre_radius + particle_radius)) {
        ParticlesType::reference f = std::get<0>(i);
        const vdouble2 &dx = std::get<1>(i);
        if (uni(get<Aboria::generator>(p)) < react_rate) {
          get<angle>(p) = std::atan2(-dx[1], -dx[0]);
          dead_particles.push_back(p);
          get<alive>(p) = false;
          ++get<count>(f);
        } else {
          get<position>(p) +=
              ((fibre_radius + particle_radius) / dx.norm() - 1) * dx;
        }
      }
      if (reflective && get<alive>(p) == true) {
        if (get<position>(p)[0] > L - particle_radius) {
          get<position>(p) += 2 * (L - particle_radius) - get<position>(p)[0];
        } else if (get<position>(p)[0] < particle_radius) {
          get<position>(p) += 2 * particle_radius - get<position>(p)[0];
        }
      }
    }
    particles.update_positions();

    t1 = Clock::now();
    time_vel_rest += (t1 - t0).count();
  }
}
