#include "filter.h"
#include "setup_knots.h"
#include "solve_stokes_MAPS.h"


int main(int argc, char **argv) {
    unsigned int nout,max_iter_linear,restart_linear,nx,nc0;
    int fibre_number,seed;
    double fibre_radius,particle_rate,react_rate,D,fibre_resolution;
    double dt_aim,h0_factor,k,gamma,rf,c0_min,c0_max,epsilon_strength,epsilon_falloff;
    unsigned int solver_in;
    std::string filename;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("max_iter_linear", po::value<unsigned int>(&max_iter_linear)->default_value(2000), "maximum iterations for linear solve")
        ("restart_linear", po::value<unsigned int>(&restart_linear)->default_value(2001), "iterations until restart for linear solve")
        ("linear_solver", po::value<unsigned int>(&solver_in)->default_value(2), "linear solver")
        ("nout", po::value<unsigned int>(&nout)->default_value(100), "number of output points")
        ("k", po::value<double>(&k)->default_value(1.00), "spring constant")
        ("D", po::value<double>(&D)->default_value(0.01), "diffusion constant")
        ("particle_rate", po::value<double>(&particle_rate)->default_value(1000.0), "particle rate")
        ("react_rate", po::value<double>(&react_rate)->default_value(0.5), "particle reaction rate")
        ("epsilon_strength", po::value<double>(&epsilon_strength)->default_value(1), "boundary clustering strength")
        ("epsilon_falloff", po::value<double>(&epsilon_falloff)->default_value(0.3), "boundary clustering fall-off")
        ("c0_min", po::value<double>(&c0_min)->default_value(0.1), "kernel constant")
        ("c0_max", po::value<double>(&c0_max)->default_value(0.1), "kernel constant")
        ("nx", po::value<unsigned int>(&nx)->default_value(20), "nx")
        ("nc0", po::value<unsigned int>(&nc0)->default_value(20), "number of c0 samples")
        ("fibre_resolution", po::value<double>(&fibre_resolution)->default_value(0.75), "knot resolution around fibre")
        ("seed", po::value<int>(&seed)->default_value(10), "seed")
        ("fibre_number", po::value<int>(&fibre_number)->default_value(5), "number of fibres")
        ("fibre_radius", po::value<double>(&fibre_radius)->default_value(0.3), "radius of fibres")
        ("dt", po::value<double>(&dt_aim)->default_value(0.001), "timestep")
        ("filename", po::value<std::string>(&filename)->default_value("optim.out"), "filename")
    ;


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    ComsolType comsol;
    read_data_files(comsol);

    const double c0_increment = (c0_max-c0_min)/nc0;

    std::ofstream file;
    file.open(filename.c_str());
    file   << std::setw(15) << "c0"
           << std::setw(15) << "rms_error_u"
           << std::setw(15) << "rms_error_v"
           << std::setw(15) << "rms_error_p"
           << std::setw(15) << "max_error_u"
           << std::setw(15) << "max_error_v"
           << std::setw(15) << "max_error_p"
           << std::endl;

    for (double c0 = c0_min; c0 < c0_max; c0+=c0_increment) {

      KnotsType knots;
      ParticlesType particles;
      ParticlesType fibres;

      const double c2 = std::pow(c0,2);
      const double c3 = std::pow(c0,3);
      const double c4 = std::pow(c0,4);
      const double mu = 1.0;
      const double flow_rate = 1.0;
      const double Tf = 2.0;
      const double L = fibre_number*1.0;
      const double delta = L/nx;
      const double boundary_layer = delta/5;
      const double s = 1.1*delta;
      const double h0 = h0_factor*delta;
      const int timesteps = Tf/dt_aim;
      const double dt = Tf/timesteps;
      const double dt_adapt = (1.0/100.0)*PI/sqrt(2*k);
      const double2 domain_min(0,-3);
      const double2 domain_max(L,L+1);
      const double2 ns_buffer(L/3,L/3);

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
        fibres.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,fibre_radius+boundary_layer,bool2(false));

        std::cout << "added "<<fibres.size()<<" fibres"<<std::endl;
      }

      //
      // SETUP KNOTS
      //
      setup_knots(knots, fibres, fibre_radius, fibre_resolution, nx, domain_min, domain_max, c0, k,epsilon_strength,epsilon_falloff);

      //
      // CALCULATE C
      //
      calculate_c(knots,c0,nx,domain_min,domain_max);

      max_iter_linear = knots.size()*4;
      restart_linear = max_iter_linear+1;

      //
      // SOLVE STOKES
      //
      solve_stokes_MAPS(knots,max_iter_linear,restart_linear,solver_in,c0);
      //solve_stokes_LMAPS(knots,max_iter_linear,restart_linear,solver_in,c0);

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

      Symbol<dvelocity_u> dvu;
      Symbol<dvelocity_v> dvv;
      Symbol<dpressure> dpr;
      Symbol<alpha> al;

      Label<0,ComsolType> i(comsol);
      Label<1,KnotsType> j(knots);
      auto dx = create_dx(i,j);
      Accumulate<std::plus<double> > sum;
      Accumulate<std::plus<double2> > sumv;
      VectorSymbolic<double,2> vector;

      Accumulate<Aboria::max<double> > max;
      max.set_init(0);

      auto psol_u1 = gen_psol_u1(i,j,c);
      auto psol_v1 = gen_psol_v1(i,j,c);
      auto psol_p1 = gen_psol_p1(i,j,c);
      auto psol_u2 = gen_psol_u2(i,j,c);
      auto psol_v2 = gen_psol_v2(i,j,c);
      auto psol_p2 = gen_psol_p2(i,j,c);

      //c[i] = c0;


      // calculate solution at comsol pts
      /*
      vu[i] = sum(j,true,gen_psol_u1(i,j,c)*al[j][0] + gen_psol_u2(i,j,c)*al[j][1]);
      vv[i] = sum(j,true,gen_psol_v1(i,j,c)*al[j][0] + gen_psol_v2(i,j,c)*al[j][1]);
      pr[i] = sum(j,true,gen_psol_p1(i,j,c)*al[j][0] + gen_psol_p2(i,j,c)*al[j][1]);
      */
      vu[i] = sum(j,true,psol_u1*al[j][0] + psol_u2*al[j][1]);
      vv[i] = sum(j,true,psol_v1*al[j][0] + psol_v2*al[j][1]);
      pr[i] = sum(j,true,psol_p1*al[j][0] + psol_p2*al[j][1]);

      //compare
      const double rms_error_u = std::sqrt(eval(sum(i,true,pow(vu[i]-dvu[i],2)))/comsol.size());
      const double rms_error_v = std::sqrt(eval(sum(i,true,pow(vv[i]-dvv[i],2)))/comsol.size());
      const double rms_error_p = std::sqrt(eval(sum(i,true,pow(pr[i]-dpr[i],2)))/comsol.size());
      const double max_error_u = eval(max(i,true,abs(vu[i]-dvu[i])));
      const double max_error_v = eval(max(i,true,abs(vv[i]-dvv[i])));
      const double max_error_p = eval(max(i,true,abs(pr[i]-dpr[i])));

      file << std::setw(15) << c0
           << std::setw(15) << rms_error_u
           << std::setw(15) << rms_error_v
           << std::setw(15) << rms_error_p
           << std::setw(15) << max_error_u
           << std::setw(15) << max_error_v
           << std::setw(15) << max_error_p
           << std::endl;

      std::cout << "rms errors = "
      << rms_error_u << ' '
      << rms_error_v << ' '
      << rms_error_p << ' '
      << std::endl;

      std::cout << "max errors = "
      << max_error_u << ' '
      << max_error_v << ' '
      << max_error_p << ' '
      << std::endl;
    }
    file.close();
}
