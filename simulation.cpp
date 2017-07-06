#include "filter.h"
#include "solve_stokes_MAPS.h"
#include "setup_knots.h"

#define MAPS


int main(int argc, char **argv) {

    unsigned int nout,max_iter_linear,restart_linear,nx;
    int fibre_number,seed;
    double fibre_resolution,fibre_radius,particle_rate,react_rate,D;
    double dt_aim,k,gamma,rf,c0;
    unsigned int solver_in;
    bool periodic;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("max_iter_linear", po::value<unsigned int>(&max_iter_linear)->default_value(2000), "maximum iterations for linear solve")
        ("restart_linear", po::value<unsigned int>(&restart_linear)->default_value(2001), "iterations until restart for linear solve")
        ("linear_solver", po::value<unsigned int>(&solver_in)->default_value(2), "linear solver")
        ("nout", po::value<unsigned int>(&nout)->default_value(100), "number of output points")
        ("k", po::value<double>(&k)->default_value(1.00), "spring constant")
        ("D", po::value<double>(&D)->default_value(0.01), "diffusion constant")
        ("particle_rate", po::value<double>(&particle_rate)->default_value(500.0), "particle rate")
        ("periodic", po::value<bool>(&periodic)->default_value(false), "periodic in x")
        ("react_rate", po::value<double>(&react_rate)->default_value(0.5), "particle reaction rate")

        ("c0", po::value<double>(&c0)->default_value(0.0835), "kernel constant")
        ("nx", po::value<unsigned int>(&nx)->default_value(10), "nx")
        ("fibre_resolution", po::value<double>(&fibre_resolution)->default_value(0.2), "number of knots around each fibre")
        ("fibre_radius", po::value<double>(&fibre_radius)->default_value(0.3), "radius of fibres")
        ("seed", po::value<int>(&seed)->default_value(10), "seed")
        ("fibre_number", po::value<int>(&fibre_number)->default_value(5), "number of fibres")
        ("dt", po::value<double>(&dt_aim)->default_value(0.001), "timestep")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }


    KnotsType knots;
    ParticlesType particles;
    ParticlesType fibres;

      const double c2 = std::pow(c0,2);
      const double c3 = std::pow(c0,3);
      const double c4 = std::pow(c0,4);
      const double mu = 1.0;
      const double flow_rate = 1.0;
      const double L = fibre_number*1.0;
      const double Tf = 2.0*L;
      const double delta = L/nx;
      const double boundary_layer = delta/5;
      const double s = 1.1*delta;
      const int timesteps = Tf/dt_aim;
      const double dt = Tf/timesteps;
      const double dt_adapt = (1.0/100.0)*PI/sqrt(2*k);
      const double2 domain_min(0,-1);
      const double2 domain_max(L,L+1);
      const double2 ns_buffer(L/3,0.5);

      std::default_random_engine generator(seed);


      // fibres
      {
        typename ParticlesType::value_type p;
        for (int ii=0; ii<fibre_number; ++ii) {
          for (int jj=0; jj<fibre_number; ++jj) {
            const double dx = L/fibre_number;
            const double2 origin = double2(
              (ii+0.5)*dx,
              (jj+0.5)*dx
            );
            get<position>(p) = origin;
            fibres.push_back(p);
          }
        }
        fibres.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,bool2(false));
        particles.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,bool2(false));

        std::cout << "added "<<fibres.size()<<" fibres"<<std::endl;
      }

      //
      // SETUP KNOTS
      //
      setup_knots(knots, fibres, fibre_radius, fibre_resolution, nx, domain_min, domain_max, c0, k,periodic,10);

      //
      // CALCULATE C
      //
      calculate_c(knots,c0,nx,domain_min,domain_max);

      max_iter_linear = knots.size()*4;
      restart_linear = max_iter_linear+1;

      //
      // SOLVE STOKES
      //
#ifdef MAPS
      const double relative_error = solve_stokes_MAPS(knots,max_iter_linear,restart_linear,solver_in,c0);
#endif
#ifdef COMPACT
      const double relative_error = solve_stokes_Compact(knots,max_iter_linear,restart_linear,solver_in,c0);
#endif
      //solve_stokes_fMAPS(knots,max_iter_linear,restart_linear,solver_in,c0,ncheb);
      //solve_stokes_LMAPS(knots,max_iter_linear,restart_linear,solver_in,c0);
      //
      
      const size_t Nk = knots.size();

      typedef typename position::value_type const & const_position_reference;
      typedef typename KnotsType::const_reference const_knot_reference;
      typedef typename KnotsType::reference knot_reference;

#ifdef MAPS
      auto psol_u1_kernel = 
            [&](const double2& dx, const double2&, const double2&) {
                return psol_u1(dx,c0);
            };
      auto psol_u2_kernel = 
            [&](const double2& dx, const double2&, const double2&) {
                return psol_u2(dx,c0);
            };
      auto psol_v1_kernel = 
            [&](const double2& dx, const double2&, const double2&) {
                return psol_v1(dx,c0);
            };
      auto psol_v2_kernel = 
            [&](const double2& dx, const double2&, const double2&) {
                return psol_v2(dx,c0);
            };
      auto psol_p1_kernel = 
            [&](const double2& dx, const double2&, const double2&) {
                return psol_p1(dx,c0);
            };
      auto psol_p2_kernel = 
            [&](const double2& dx, const double2&, const double2&) {
                return psol_p2(dx,c0);
            };

#endif

      std::cout << "making fmm_queries" <<std::endl;
      auto psol_u1_fmm = make_fmm_query(knots.get_query(),make_black_box_expansion<2,8>(psol_u1_kernel));
      auto psol_u2_fmm = make_fmm_query(knots.get_query(),make_black_box_expansion<2,8>(psol_u2_kernel));
      auto psol_v1_fmm = make_fmm_query(knots.get_query(),make_black_box_expansion<2,8>(psol_v1_kernel));
      auto psol_v2_fmm = make_fmm_query(knots.get_query(),make_black_box_expansion<2,8>(psol_v2_kernel));

      std::cout << "calculating expansions" <<std::endl;
      psol_u1_fmm.calculate_expansions(get<alpha1>(knots));
      psol_u2_fmm.calculate_expansions(get<alpha2>(knots));
      psol_v1_fmm.calculate_expansions(get<alpha1>(knots));
      psol_v2_fmm.calculate_expansions(get<alpha2>(knots));
      std::cout << "done calculating expansions" <<std::endl;

      Symbol<boundary> is_b;
      Symbol<inlet> is_in;
      Symbol<outlet> is_out;
      Symbol<interior> is_i;
      Symbol<position> r;
      Symbol<alive> alive_;
      Symbol<kernel_constant> c;

      Symbol<velocity_u> vu;
      Symbol<velocity_v> vv;
      Symbol<pressure> pr;
      Symbol<alpha1> al1;
      Symbol<alpha2> al2;

      Label<0,KnotsType> i(knots);
      Label<1,KnotsType> j(knots);
      Label<0,ParticlesType> a(particles);
      Label<1,ParticlesType> b(particles);
      Label<0,ParticlesType> af(fibres);
      Label<1,ParticlesType> bf(fibres);
      auto dx = create_dx(i,j);
      auto dkf = create_dx(i,bf);
      auto dpf = create_dx(a,bf);
      auto dpk = create_dx(a,j);
      AccumulateWithinDistance<std::plus<double2> > sumv(fibre_radius);
      AccumulateWithinDistance<std::bit_or<bool> > any(fibre_radius);
      any.set_init(false);
      VectorSymbolic<double,2> vector;
      Normal N;
      Uniform U;

#undef NDEBUG
      std::cout << "starting timesteps!"<<std::endl;
      const int timesteps_per_out = timesteps/nout;
      auto t0 = Clock::now();
      auto t1 = Clock::now();
      double time_vel_eval = 0;
      double time_vel_rest = 0;
      for (int ii=0; ii<timesteps; ++ii) {
          // add new particles
          std::poisson_distribution<int> poisson(particle_rate*dt);
          std::uniform_real_distribution<double> uniform(L/10.0,L-L/10.0);
          const int new_n = poisson(generator);
          for (int jj=0; jj<new_n; ++jj) {
              ParticlesType::value_type p;
              get<position>(p) = double2(uniform(generator),domain_max[1]);
              get<kernel_constant>(p) = c0;
              particles.push_back(p);
          }

          // evaluate velocity field
          
          t0 = Clock::now();
          for (ParticlesType::reference p: particles) {
              get<velocity_u>(p) = psol_u1_fmm.evaluate_expansion(get<position>(p),
                      get<alpha1>(knots)) 
                  + psol_u2_fmm.evaluate_expansion(get<position>(p),
                          get<alpha2>(knots));
              get<velocity_v>(p) = psol_v1_fmm.evaluate_expansion(get<position>(p),
                      get<alpha1>(knots)) 
                  + psol_v2_fmm.evaluate_expansion(get<position>(p),
                          get<alpha2>(knots));
          }
          t1 = Clock::now();
          time_vel_eval += (t1 - t0).count();
          t0 = Clock::now();

          r[a] += std::sqrt(2.0*D*dt)*vector(N[a],N[a]) + dt*vector(vu[a],vv[a]);

          if (ii % timesteps_per_out == 0) {
              std::cout << "timestep "<<ii<<" of "<<timesteps<<" (time_vel_eval = "<<time_vel_eval<<" time_vel_rest = "<<time_vel_rest<<std::endl;
              vtkWriteGrid("particles",ii,particles.get_grid(true));
          }

          // react with fibres
          alive_[a] = !any(bf,U[a]<react_rate);

          // react with side walls
          alive_[a] = !if_else(r[a][0] > L
                  ,U[a]<react_rate
                  ,if_else(r[a][0] < 0
                      ,U[a]<react_rate
                      ,false
                      )
                  );

          // reflect off fibres (if still alive)
          r[a] += sumv(bf,(fibre_radius/norm(dpf)-1)*dpf);

          // reflect off side walls (if still alive)
          r[a] = vector(
                  if_else(r[a][0] > L
                      ,2*L-r[a][0]
                      ,if_else(r[a][0] < 0
                          ,-r[a][0]
                          ,r[a][0]
                          )
                      )
                  ,r[a][1]
                  );

          t1 = Clock::now();
          time_vel_rest += (t1 - t0).count();

      }

}
