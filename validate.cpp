#include "filter.h"
#include "solve_stokes_MAPS.h"
#include "solve_stokes_LMAPS.h"
#include "setup_knots.h"

int main(int argc, char **argv) {
    unsigned int nout,max_iter_linear,restart_linear,nx_min,nx_max,nnx;
    int fibre_number,seed;
    double fibre_radius,particle_rate,react_rate,D,fibre_resolution;
    double dt_aim,h0_factor,k,gamma,rf,c0,epsilon_strength,epsilon_falloff;
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
        ("c0", po::value<double>(&c0)->default_value(0.1), "kernel constant")
        ("nx_min", po::value<unsigned int>(&nx_min)->default_value(20), "nx min")
        ("nx_max", po::value<unsigned int>(&nx_max)->default_value(21), "nx max")
        ("nnx", po::value<unsigned int>(&nnx)->default_value(1), "number of nx samples")
        ("fibre_resolution", po::value<double>(&fibre_resolution)->default_value(0.75), "knot resolution around fibre")
        ("seed", po::value<int>(&seed)->default_value(10), "seed")
        ("fibre_number", po::value<int>(&fibre_number)->default_value(5), "number of fibres")
        ("fibre_radius", po::value<double>(&fibre_radius)->default_value(0.3), "radius of fibres")
        ("dt", po::value<double>(&dt_aim)->default_value(0.001), "timestep")
        ("filename", po::value<std::string>(&filename)->default_value("validate.out"), "filename")
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

    std::ofstream file;
    file.open(filename.c_str());
    file   << std::setw(15) << "nx"
           << std::setw(15) << "rms_error_u"
           << std::setw(15) << "rms_error_v"
           << std::setw(15) << "rms_error_p"
           << std::setw(15) << "max_error_u"
           << std::setw(15) << "max_error_v"
           << std::setw(15) << "max_error_p"
           << std::endl;

    const size_t nx_increment = (nx_max-nx_min)/nnx;
    for (size_t nx = nx_min; nx < nx_max; nx+=nx_increment) {

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
      const double2 domain_min(0,-1);
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
        fibres.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,bool2(false));

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
      //solve_stokes_MAPS(knots,max_iter_linear,restart_linear,solver_in,c0);
      solve_stokes_LMAPS(knots,max_iter_linear,restart_linear,solver_in,c0);

      
      const size_t Nk = knots.size();
      const size_t Nc = comsol.size();


      typedef typename position::value_type const & const_position_reference;
      typedef typename KnotsType::const_reference const_knot_reference;
      typedef typename ComsolType::const_reference const_comsol_reference;

      auto psol_u1_op = create_dense_operator(comsol,knots,
            [](const_position_reference dx,
               const_comsol_reference a,
               const_knot_reference b) {
            return psol_u1(dx,get<kernel_constant>(b));
            });

      auto psol_v1_op = create_dense_operator(comsol,knots,
            [](const_position_reference dx,
               const_comsol_reference a,
               const_knot_reference b) {
            return psol_v1(dx,get<kernel_constant>(b));
            });

      auto psol_p1_op = create_dense_operator(comsol,knots,
            [](const_position_reference dx,
               const_comsol_reference a,
               const_knot_reference b) {
            return psol_p1(dx,get<kernel_constant>(b));
            });

      auto psol_u2_op = create_dense_operator(comsol,knots,
            [](const_position_reference dx,
               const_comsol_reference a,
               const_knot_reference b) {
            return psol_u2(dx,get<kernel_constant>(b));
            });

      auto psol_v2_op = create_dense_operator(comsol,knots,
            [](const_position_reference dx,
               const_comsol_reference a,
               const_knot_reference b) {
            return psol_v2(dx,get<kernel_constant>(b));
            });

      auto psol_p2_op = create_dense_operator(comsol,knots,
            [](const_position_reference dx,
               const_comsol_reference a,
               const_knot_reference b) {
            return psol_p2(dx,get<kernel_constant>(b));
            });

      std::cout << "assembling comsol matricies...." << std::endl;

      typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
      typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
      typedef Eigen::Map<vector_type> map_type;

      matrix_type psol_u1_mat(Nc,Nk);
      psol_u1_op.assemble(psol_u1_mat);
      matrix_type psol_v1_mat(Nc,Nk);
      psol_v1_op.assemble(psol_v1_mat);
      matrix_type psol_p1_mat(Nc,Nk);
      psol_p1_op.assemble(psol_p1_mat);

      matrix_type psol_u2_mat(Nc,Nk);
      psol_u2_op.assemble(psol_u2_mat);
      matrix_type psol_v2_mat(Nc,Nk);
      psol_v2_op.assemble(psol_v2_mat);
      matrix_type psol_p2_mat(Nc,Nk);
      psol_p2_op.assemble(psol_p2_mat);
      std::cout << "done assembling comsol matricies." << std::endl;


      // calculate solution at comsol pts
      /*
      vu[i] = sum(j,true,gen_psol_u1(i,j,c)*al[j][0] + gen_psol_u2(i,j,c)*al[j][1]);
      vv[i] = sum(j,true,gen_psol_v1(i,j,c)*al[j][0] + gen_psol_v2(i,j,c)*al[j][1]);
      pr[i] = sum(j,true,gen_psol_p1(i,j,c)*al[j][0] + gen_psol_p2(i,j,c)*al[j][1]);
      */

      std::cout << "calculating solution at comsol points...." << std::endl;
    map_type(get<velocity_u>(comsol).data(),Nc) = 
        psol_u1_mat*map_type(get<alpha1>(knots).data(),Nk)+
        psol_u2_mat*map_type(get<alpha2>(knots).data(),Nk);

    map_type(get<velocity_v>(comsol).data(),Nc) = 
        psol_v1_mat*map_type(get<alpha1>(knots).data(),Nk)+
        psol_v2_mat*map_type(get<alpha2>(knots).data(),Nk);

    map_type(get<pressure>(comsol).data(),Nc) = 
        psol_v1_mat*map_type(get<alpha1>(knots).data(),Nk)+
        psol_v2_mat*map_type(get<alpha2>(knots).data(),Nk);
      std::cout << "done calculating solution at comsol points.<< std::endl;

      //compare
      Symbol<velocity_u> vu;
      Symbol<velocity_v> vv;
      Symbol<pressure> pr;

      Symbol<dvelocity_u> dvu;
      Symbol<dvelocity_v> dvv;
      Symbol<dpressure> dpr;

      Label<0,ComsolType> i(comsol);
      Accumulate<std::plus<double> > sum;
      Accumulate<Aboria::max<double> > max;
      max.set_init(0);

      const double rms_error_u = std::sqrt(eval(sum(i,true,pow(vu[i]-dvu[i],2)))/comsol.size());
      const double rms_error_v = std::sqrt(eval(sum(i,true,pow(vv[i]-dvv[i],2)))/comsol.size());
      const double rms_error_p = std::sqrt(eval(sum(i,true,pow(pr[i]-dpr[i],2)))/comsol.size());
      const double max_error_u = eval(max(i,true,abs(vu[i]-dvu[i])));
      const double max_error_v = eval(max(i,true,abs(vv[i]-dvv[i])));
      const double max_error_p = eval(max(i,true,abs(pr[i]-dpr[i])));

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

      file << std::setw(15) << nx
           << std::setw(15) << rms_error_u
           << std::setw(15) << rms_error_v
           << std::setw(15) << rms_error_p
           << std::setw(15) << max_error_u
           << std::setw(15) << max_error_v
           << std::setw(15) << max_error_p
           << std::endl;


      vtkWriteGrid("error",nx,comsol.get_grid(true));
    }
    file.close();
}
