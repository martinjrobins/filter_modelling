#include "filter.h"
#include "setup_knots.h"
#include <fstream>


//#define FMAPS
//#define MAPS
//#define COMPACT 

#include "solve_stokes_BEM.h"



int main(int argc, char **argv) {
    unsigned int nout,max_iter_linear,restart_linear,nx,nmin,nmax;
    int fibre_number,seed;
    double fibre_radius,particle_rate,react_rate,D,fibre_resolution;
    double dt_aim,h0_factor,k,gamma,rf;
    unsigned int solver_in;
    bool periodic;
    std::string filename = "optim_default.out";

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("max_iter_linear", po::value<unsigned int>(&max_iter_linear)->default_value(2000), "maximum iterations for linear solve")
        ("restart_linear", po::value<unsigned int>(&restart_linear)->default_value(2001), "iterations until restart for linear solve")
        ("linear_solver", po::value<unsigned int>(&solver_in)->default_value(2), "linear solver")
        ("nout", po::value<unsigned int>(&nout)->default_value(100), "number of output points")
        ("k", po::value<double>(&k)->default_value(1.00), "spring constant")
        ("periodic", po::value<bool>(&periodic)->default_value(true), "periodic in x")
        ("D", po::value<double>(&D)->default_value(0.01), "diffusion constant")
        ("particle_rate", po::value<double>(&particle_rate)->default_value(1000.0), "particle rate")
        ("react_rate", po::value<double>(&react_rate)->default_value(0.5), "particle reaction rate")
        ("n_min", po::value<unsigned int>(&nmin)->default_value(2), "min n")
        ("n_max", po::value<unsigned int>(&nmax)->default_value(3), "max n")
        ("nx", po::value<unsigned int>(&nx)->default_value(20), "nx")
        ("fibre_resolution", po::value<double>(&fibre_resolution)->default_value(0.75), "knot resolution around fibre")
        ("seed", po::value<int>(&seed)->default_value(10), "seed")
        ("fibre_number", po::value<int>(&fibre_number)->default_value(5), "number of fibres")
        ("fibre_radius", po::value<double>(&fibre_radius)->default_value(0.3), "radius of fibres")
        ("dt", po::value<double>(&dt_aim)->default_value(0.001), "timestep")
        //("filename", po::value<std::string>(&filename)->default_value("optim.out"), "filename")
    ;


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    ComsolType comsol;
    read_data_files(comsol);


    std::ofstream file;
    file.open(filename.c_str());
    file   << std::setw(15) << "n"
           << std::setw(15) << "rms_error"
           << std::setw(15) << "max_error"
           << std::setw(15) << "solve_error"
           << std::endl;

    for (unsigned int n= nmin; n < nmax ; ++n) {

      KnotsType knots;
      ElementsType elements;
      ParticlesType particles;
      ParticlesType fibres;

      const double mu = 1.0;
      const double flow_rate = 1.0;
      const double Tf = 2.0;
      const double L = fibre_number*1.0;
      const int timesteps = Tf/dt_aim;
      const double dt = Tf/timesteps;
      const double dt_adapt = (1.0/100.0)*PI/sqrt(2*k);
      const vdouble2 domain_min(0,-1);
      const vdouble2 domain_max(L,L+1);
      const vdouble2 ns_buffer(L/3,L/3);

      std::default_random_engine generator(seed);


      // fibres
      {
        typename ParticlesType::value_type p;
        for (int ii=0; ii<fibre_number; ++ii) {
          for (int jj=0; jj<fibre_number; ++jj) {
            const double dx = L/fibre_number;
            const vdouble2 origin = vdouble2(
              (ii+0.5)*dx,
              (jj+0.5)*dx
            );
            get<position>(p) = origin;
            fibres.push_back(p);
          }
        }
        fibres.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,vbool2(false));

        std::cout << "added "<<fibres.size()<<" fibres"<<std::endl;
      }

      //
      // SETUP KNOTS
      //
      //setup_knots(knots, fibres, fibre_radius, fibre_resolution, nx, domain_min, domain_max, k,periodic,8*8);

        //
        // SETUP ELEMENTS 
        //
        setup_elements(elements, fibres, domain_min, domain_max, nx, fibre_radius);


        max_iter_linear = knots.size()*4;
        restart_linear = max_iter_linear+1;

        //
        // SOLVE STOKES
        //
        const double alpha = (elements.get_max()-elements.get_min()).prod()/(4.0*PI);
        const double h = (get<point_b>(elements)[0] - get<point_a>(elements)[0]).norm();
        const int nlambda = n;
        const int nmu = n;
        const double relative_error = solve_stokes_BEM(knots, elements, alpha, nlambda, nmu);
 

      const size_t Nk = knots.size();
      const size_t Nc = comsol.size();

      typedef typename position::value_type const & const_position_reference;
      typedef typename KnotsType::const_reference const_knot_reference;
      typedef typename KnotsType::reference knot_reference;
      typedef typename ComsolType::const_reference const_comsol_reference;
      typedef typename ComsolType::reference comsol_reference;


      std::cout << "alpha = "<<alpha << std::endl;
      auto Akernel = make_greens_kernel_2d1p(alpha,nlambda,nmu,h,
                                   domain_min,domain_max,false);
      auto A = create_dense_operator(comsol,elements,
                Akernel);
     

      std::cout << "calculating solution at comsol points...." << std::endl;
      A.get_first_kernel().evaluate(get<velocity>(comsol),get<traction>(elements));
      for(auto p:comsol) {
          get<velocity>(p) *= -1.0/(4.0*PI*mu);
          get<velocity>(p)[1] -= flow_rate;
      }
      std::cout << "done calculating solution at comsol points."<< std::endl;

      double rms_error_v = 
            std::sqrt(std::accumulate(std::begin(comsol),std::end(comsol),0.0,
                [](double accum, const_comsol_reference p) {
                    return accum + (get<velocity>(p)-get<dvelocity>(p)).squaredNorm();
                }));

      const double max_error_v = 
            std::accumulate(std::begin(comsol),std::end(comsol),0.0,
                [](double accum, const_comsol_reference p) {
                    return std::max(accum,std::max((get<velocity>(p)-get<dvelocity>(p))[0],(get<velocity>(p)-get<dvelocity>(p))[1]));
                });


      file << std::setw(15) << n
           << std::setw(15) << rms_error_v
           << std::setw(15) << max_error_v
           << std::setw(15) << relative_error
           << std::endl;

      std::cout << "rms errors = "
      << rms_error_v << ' '
      << std::endl;

      std::cout << "max errors = "
      << max_error_v << ' '
      << std::endl;

      vtkWriteGrid("error",0,comsol.get_grid(true));
    }
    file.close();
}
