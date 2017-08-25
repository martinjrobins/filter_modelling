#include "filter.h"
#include "setup_knots.h"
#include <fstream>


//#define FMAPS
#define MAPS
//#define COMPACT 

#ifdef FMAPS
#include "solve_stokes_fMAPS.h"
#endif
#ifdef COMPACT 
#include "solve_stokes_Compact.h"
#endif
#ifdef MAPS
#include "solve_stokes_MAPS.h"
#endif
//#include "solve_stokes_LMAPS.h"


template <unsigned int N, typename U1K, typename U2K, 
                          typename V1K, typename V2K, 
                          typename P1K, typename P2K>
void eval_solution(KnotsType& knots, ComsolType& comsol, 
                   U1K& u1k, U2K& u2k,V1K& v1k, V2K& v2k,P1K& p1k, P2K& p2k,
                   double& rms_error_u, double& rms_error_v, double& rms_error_p,
                   double& rms_diff_u, double& rms_diff_v, double& rms_diff_p,
                   double& time_setup, double& time_eval, double& time_direct_eval) {
    const unsigned int D = KnotsType::dimension; 
    auto t0 = Clock::now();
    auto psol_u1_fmm = make_fmm_with_source(knots,make_black_box_expansion<D,N>(u1k),
                                            get<alpha1>(knots));
    auto psol_u2_fmm = make_fmm_with_source(knots,make_black_box_expansion<D,N>(u2k),
                                            get<alpha2>(knots));
    auto psol_v1_fmm = make_fmm_with_source(knots,make_black_box_expansion<D,N>(v1k),
                                            get<alpha1>(knots));
    auto psol_v2_fmm = make_fmm_with_source(knots,make_black_box_expansion<D,N>(v2k),
                                            get<alpha2>(knots));
    auto psol_p1_fmm = make_fmm_with_source(knots,make_black_box_expansion<D,N>(p1k),
                                            get<alpha1>(knots));
    auto psol_p2_fmm = make_fmm_with_source(knots,make_black_box_expansion<D,N>(p2k),
                                            get<alpha2>(knots));
    auto t1 = Clock::now();
    time_setup = (t1 - t0).count();
    t0 = Clock::now();
    for (ComsolType::reference p: comsol) {
        get<velocity_u>(p) = psol_u1_fmm.evaluate_at_point(get<position>(p),
                get<alpha1>(knots)) 
            + psol_u2_fmm.evaluate_at_point(get<position>(p),
                    get<alpha2>(knots));
        get<velocity_v>(p) = psol_v1_fmm.evaluate_at_point(get<position>(p),
                get<alpha1>(knots)) 
            + psol_v2_fmm.evaluate_at_point(get<position>(p),
                    get<alpha2>(knots));
        get<pressure>(p) =   psol_p1_fmm.evaluate_at_point(get<position>(p),
                get<alpha1>(knots)) 
            + psol_p2_fmm.evaluate_at_point(get<position>(p),
                    get<alpha2>(knots));
    }
    t1 = Clock::now();
    time_eval = (t1 - t0).count();

    t0 = Clock::now();
    rms_error_u = rms_error_v = rms_error_p = 0;
    rms_diff_u = rms_diff_v = rms_diff_p = 0;
    double max_error_u, max_error_v, max_error_p;
    max_error_u = max_error_v = max_error_p = 0;
    for (ComsolType::reference i: comsol) {
        double u,v,p;
        u = v = p = 0;
        for (KnotsType::reference j: knots) {
            const vdouble2 dx = get<position>(j)-get<position>(i);
            u += u1k(dx,get<position>(i),get<position>(j))*get<alpha1>(j);
            u += u2k(dx,get<position>(i),get<position>(j))*get<alpha2>(j);
            v += v1k(dx,get<position>(i),get<position>(j))*get<alpha1>(j);
            v += v2k(dx,get<position>(i),get<position>(j))*get<alpha2>(j);
            p += p1k(dx,get<position>(i),get<position>(j))*get<alpha1>(j);
            p += p2k(dx,get<position>(i),get<position>(j))*get<alpha2>(j);
        }
        double diff_u = std::abs(u-get<velocity_u>(i));
        double diff_v = std::abs(v-get<velocity_v>(i));
        double diff_p = std::abs(p-get<pressure>(i));
        double error_u = std::abs(get<velocity_u>(i)-get<dvelocity_u>(i));
        double error_v = std::abs(get<velocity_v>(i)-get<dvelocity_v>(i));
        double error_p = std::abs(get<pressure>(i)-get<dpressure>(i));
        if (error_u > max_error_u) max_error_u = error_u;
        if (error_v > max_error_v) max_error_v = error_v;
        if (error_p > max_error_p) max_error_p = error_p;
        rms_diff_u += std::pow(diff_u,2);
        rms_diff_v += std::pow(diff_v,2);
        rms_diff_p += std::pow(diff_p,2);
        rms_error_u += std::pow(error_u,2);
        rms_error_v += std::pow(error_v,2);
        rms_error_p += std::pow(error_p,2);
    }
    //vtkWriteGrid("fmm_eval_knots",N,knots.get_grid(true));
    //vtkWriteGrid("fmm_eval_comsol",N,comsol.get_grid(true));
    rms_error_u = std::sqrt(rms_error_u/comsol.size());
    rms_error_v = std::sqrt(rms_error_v/comsol.size());
    rms_error_p = std::sqrt(rms_error_p/comsol.size());
    rms_diff_u = std::sqrt(rms_diff_u/comsol.size());
    rms_diff_v = std::sqrt(rms_diff_v/comsol.size());
    rms_diff_p = std::sqrt(rms_diff_p/comsol.size());
    t1 = Clock::now();
    time_direct_eval = (t1 - t0).count();

    /*
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

    rms_error_u = std::sqrt(eval(sum(i,pow(vu[i]-dvu[i],2)))/comsol.size());
    rms_error_v = std::sqrt(eval(sum(i,pow(vv[i]-dvv[i],2)))/comsol.size());
    rms_error_p = std::sqrt(eval(sum(i,pow(pr[i]-dpr[i],2)))/comsol.size());
    const double max_error_u = eval(max(i,abs(vu[i]-dvu[i])));
    const double max_error_v = eval(max(i,abs(vv[i]-dvv[i])));
    const double max_error_p = eval(max(i,abs(pr[i]-dpr[i])));
    */

    std::cout << "rms errors (difference from comsol) = "
      << rms_error_u << ' '
      << rms_error_v << ' '
      << rms_error_p << ' '
      << std::endl;

    std::cout << "max errors (difference from comsol) = "
      << max_error_u << ' '
      << max_error_v << ' '
      << max_error_p << ' '
      << std::endl;

    std::cout << "fmm rms error (difference from direct eval)= "
      << rms_diff_u << ' '
      << rms_diff_v << ' '
      << rms_diff_p << ' '
      << std::endl;

    
    std::cout << "time for setup = " << time_setup << std::endl;
    std::cout << "time for eval = " << time_eval << std::endl;
    std::cout << "time for direct eval = " << time_direct_eval << std::endl;



}


int main(int argc, char **argv) {
    unsigned int nout,max_iter_linear,restart_linear,nx;
    int fibre_number,seed,nbucket_min,nbucket_max,ncheb_min,ncheb_max,fibre_arrangement;
    double fibre_radius,particle_rate,react_rate,D,fibre_resolution,particle_radius;
    double dt_aim,h0_factor,k,gamma,rf,c0;
    unsigned int solver_in;
    bool periodic;
    std::string filename;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("max_iter_linear", po::value<unsigned int>(&max_iter_linear)->default_value(2000), "maximum iterations for linear solve")
        ("restart_linear", po::value<unsigned int>(&restart_linear)->default_value(2001), "iterations until restart for linear solve")
        ("linear_solver", po::value<unsigned int>(&solver_in)->default_value(2), "linear solver")
        ("nout", po::value<unsigned int>(&nout)->default_value(100), "number of output points")
        ("k", po::value<double>(&k)->default_value(1.00), "spring constant")
        ("periodic", po::value<bool>(&periodic)->default_value(false), "periodic in x")
        ("D", po::value<double>(&D)->default_value(0.01), "diffusion constant")
        ("particle_rate", po::value<double>(&particle_rate)->default_value(1000.0), "particle rate")
        ("react_rate", po::value<double>(&react_rate)->default_value(0.5), "particle reaction rate")
        ("nbucket_min", po::value<int>(&nbucket_min)->default_value(10), "number of points in bucket")
        ("nbucket_max", po::value<int>(&nbucket_max)->default_value(100), "number of points in bucket max")
        ("ncheb_min", po::value<int>(&ncheb_min)->default_value(3), "number of cheb points")
        ("ncheb_max", po::value<int>(&ncheb_max)->default_value(10), "number of cheb points max")
        ("nx", po::value<unsigned int>(&nx)->default_value(10), "nx")
        ("fibre_resolution", po::value<double>(&fibre_resolution)->default_value(0.2), "knot resolution around fibre")
        ("seed", po::value<int>(&seed)->default_value(10), "seed")
        ("fibre_number", po::value<int>(&fibre_number)->default_value(5), "number of fibres")
        ("fibre_radius", po::value<double>(&fibre_radius)->default_value(0.3), "radius of fibres")
        ("particle_radius", po::value<double>(&particle_radius)->default_value(0.3/100.0), "radius of fibres")
        ("fibre_arrangement", po::value<int>(&fibre_arrangement)->default_value(0), "(0=regular, 1=hexigonal 2=random)")
        ("c0", po::value<double>(&c0)->default_value(0.0835), "kernel constant")
        ("dt", po::value<double>(&dt_aim)->default_value(0.001), "timestep")
        ("filename", po::value<std::string>(&filename)->default_value("fmm_eval.out"), "filename")
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
    file   << std::setw(15) << "nbucket"
           << std::setw(15) << "rms_diff_u5"
           << std::setw(15) << "rms_diff_u6"
           << std::setw(15) << "rms_diff_u7"
           << std::setw(15) << "rms_diff_u8"
           << std::setw(15) << "rms_diff_u9"
           << std::setw(15) << "rms_diff_u10"
           << std::setw(15) << "rms_diff_v5"
           << std::setw(15) << "rms_diff_v6"
           << std::setw(15) << "rms_diff_v7"
           << std::setw(15) << "rms_diff_v8"
           << std::setw(15) << "rms_diff_v9"
           << std::setw(15) << "rms_diff_v10"
           << std::setw(15) << "rms_diff_p5"
           << std::setw(15) << "rms_diff_p6"
           << std::setw(15) << "rms_diff_p7"
           << std::setw(15) << "rms_diff_p8"
           << std::setw(15) << "rms_diff_p9"
           << std::setw(15) << "rms_diff_p10"
           << std::setw(15) << "rms_error_u5"
           << std::setw(15) << "rms_error_u6"
           << std::setw(15) << "rms_error_u7"
           << std::setw(15) << "rms_error_u8"
           << std::setw(15) << "rms_error_u9"
           << std::setw(15) << "rms_error_u10"
           << std::setw(15) << "rms_error_v5"
           << std::setw(15) << "rms_error_v6"
           << std::setw(15) << "rms_error_v7"
           << std::setw(15) << "rms_error_v8"
           << std::setw(15) << "rms_error_v9"
           << std::setw(15) << "rms_error_v10"
           << std::setw(15) << "rms_error_p5"
           << std::setw(15) << "rms_error_p6"
           << std::setw(15) << "rms_error_p7"
           << std::setw(15) << "rms_error_p8"
           << std::setw(15) << "rms_error_p9"
           << std::setw(15) << "rms_error_p10"
           << std::setw(15) << "timet_setup5"
           << std::setw(15) << "timet_setup6"
           << std::setw(15) << "timet_setup7"
           << std::setw(15) << "timet_setup8"
           << std::setw(15) << "timet_setup9"
           << std::setw(15) << "timet_setup10"
           << std::setw(15) << "timet_eval5"
           << std::setw(15) << "timet_eval6"
           << std::setw(15) << "timet_eval7"
           << std::setw(15) << "timet_eval8"
           << std::setw(15) << "timet_eval9"
           << std::setw(15) << "timet_eval10"
           << std::setw(15) << "timet_direct_eval5"
           << std::setw(15) << "timet_direct_eval6"
           << std::setw(15) << "timet_direct_eval7"
           << std::setw(15) << "timet_direct_eval8"
           << std::setw(15) << "timet_direct_eval9"
           << std::setw(15) << "timet_direct_eval10"
           << std::endl;

    for (int nbucket = nbucket_min; nbucket < nbucket_max; 
            nbucket += static_cast<int>((nbucket_max-nbucket_min)/5.0)) {
      std::cout << "nbucket = "<<nbucket<<std::endl;
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
      const bool reflective = false;
      const double delta = L/nx;
      const double boundary_layer = delta/5;
      const double s = 1.1*delta;
      const double h0 = h0_factor*delta;
      const int timesteps = Tf/dt_aim;
      const double dt = Tf/timesteps;
      const double dt_adapt = (1.0/100.0)*PI/sqrt(2*k);
      const vdouble2 domain_min(0,-1);
      const vdouble2 domain_max(L,L+1);

      const vdouble2 ns_buffer_fibres(L/3,2*particle_radius);
      const vdouble2 ns_buffer_particles(static_cast<int>(reflective)*L/3,2*particle_radius);

      std::default_random_engine generator(seed);


      // fibres
      {

        std::uniform_real_distribution<double> xrange(domain_min[0],domain_max[0]);
        std::uniform_real_distribution<double> yrange(domain_min[1]+1+fibre_radius,domain_max[1]-1-fibre_radius);
        typename ParticlesType::value_type p;
        if (fibre_arrangement == 2) {
            for (int jj=0; jj<fibre_number; ++jj) {
                bool free_position = false;
                while (free_position == false) {
                    get<position>(p) = vdouble2(xrange(generator),yrange(generator));
                    free_position = true;
                    /*
                    for (auto tpl: euclidean_search(fibres.get_query(),
                                get<position>(p),
                                2*fibre_radius+2*particle_radius)) {
                                */
                    for (auto f: fibres) {
                        if ((get<position>(p)-get<position>(f)).norm() < 2*fibre_radius+2*particle_radius) {
                            free_position = false;
                            break;
                        }
                    }
                }
                fibres.push_back(p);

            }
        } else if (fibre_arrangement == 1) {
            for (int jj=0; jj<fibre_number; ++jj) {
              for (int ii=0; ii<((jj%2==1)?fibre_number+1:fibre_number); ++ii) {
                const double dx = L/fibre_number;
                vdouble2 origin = vdouble2(
                  ((jj%2==1)?ii:ii+0.5)*dx,
                  (jj+0.5)*dx
                );
                get<position>(p) = origin;
                fibres.push_back(p);
              }
            }
        } else {
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
        }

        particles.init_neighbour_search(domain_min-ns_buffer_particles,domain_max+ns_buffer_particles,vbool2(!reflective,false));

        fibres.init_neighbour_search(domain_min-ns_buffer_fibres,domain_max+ns_buffer_fibres,vbool2(false));

        std::cout << "added "<<fibres.size()<<" fibres"<<std::endl;

        //vtkWriteGrid("init_knots_fibres",0,fibres.get_grid(true));
      }


      //
      // SETUP KNOTS
      //
      setup_knots(knots, fibres, fibre_radius, fibre_resolution, nx, domain_min, domain_max, c0, k,periodic,nbucket);

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
      const size_t Nc = comsol.size();

      typedef typename position::value_type const & const_position_reference;
      typedef typename KnotsType::const_reference const_knot_reference;
      typedef typename KnotsType::reference knot_reference;
      typedef typename ComsolType::const_reference const_comsol_reference;
      typedef typename ComsolType::reference comsol_reference;

      for (comsol_reference i: comsol) {
          get<kernel_constant>(i) = c0;
      }

      
#ifdef MAPS
      auto psol_u1_kernel = 
            [&](const vdouble2& dx, const vdouble2&, const vdouble2&) {
                return psol_u1(dx,c0);
            };
      auto psol_u2_kernel = 
            [&](const vdouble2& dx, const vdouble2&, const vdouble2&) {
                return psol_u2(dx,c0);
            };
      auto psol_v1_kernel = 
            [&](const vdouble2& dx, const vdouble2&, const vdouble2&) {
                return psol_v1(dx,c0);
            };
      auto psol_v2_kernel = 
            [&](const vdouble2& dx, const vdouble2&, const vdouble2&) {
                return psol_v2(dx,c0);
            };
      auto psol_p1_kernel = 
            [&](const vdouble2& dx, const vdouble2&, const vdouble2&) {
                return psol_p1(dx,c0);
            };
      auto psol_p2_kernel = 
            [&](const vdouble2& dx, const vdouble2&, const vdouble2&) {
                return psol_p2(dx,c0);
            };
#endif


      std::cout << "calculating solution at comsol points...." << std::endl;
      double rms_error_u[6], rms_error_v[6], rms_error_p[6]; 
      double rms_diff_u[6], rms_diff_v[6], rms_diff_p[6]; 
      double time_eval[6], time_setup[6], time_direct_eval[6];
      /*
      eval_solution<5>(knots, comsol, 
                    psol_u1_kernel, psol_u2_kernel,
                    psol_v1_kernel, psol_v2_kernel,
                    psol_p1_kernel, psol_p2_kernel,
                    rms_error_u[0], rms_error_v[0], rms_error_p[0], 
                    rms_diff_u[0], rms_diff_v[0], rms_diff_p[0], 
                    time_setup[0], time_eval[0], time_direct_eval[0]);
      eval_solution<6>(knots, comsol, 
                    psol_u1_kernel, psol_u2_kernel,
                    psol_v1_kernel, psol_v2_kernel,
                    psol_p1_kernel, psol_p2_kernel,
                    rms_error_u[1], rms_error_v[1], rms_error_p[1], 
                    rms_diff_u[1], rms_diff_v[1], rms_diff_p[1], 
                    time_setup[1], time_eval[1], time_direct_eval[1]);
      eval_solution<7>(knots, comsol, 
                    psol_u1_kernel, psol_u2_kernel,
                    psol_v1_kernel, psol_v2_kernel,
                    psol_p1_kernel, psol_p2_kernel,
                    rms_error_u[2], rms_error_v[2], rms_error_p[2], 
                    rms_diff_u[2], rms_diff_v[2], rms_diff_p[2], 
                    time_setup[2], time_eval[2], time_direct_eval[2]);
                    */
      eval_solution<8>(knots, comsol, 
                    psol_u1_kernel, psol_u2_kernel,
                    psol_v1_kernel, psol_v2_kernel,
                    psol_p1_kernel, psol_p2_kernel,
                    rms_error_u[3], rms_error_v[3], rms_error_p[3], 
                    rms_diff_u[3], rms_diff_v[3], rms_diff_p[3], 
                    time_setup[3], time_eval[3], time_direct_eval[3]);
      eval_solution<9>(knots, comsol, 
                    psol_u1_kernel, psol_u2_kernel,
                    psol_v1_kernel, psol_v2_kernel,
                    psol_p1_kernel, psol_p2_kernel,
                    rms_error_u[4], rms_error_v[4], rms_error_p[4], 
                    rms_diff_u[4], rms_diff_v[4], rms_diff_p[4], 
                    time_setup[4], time_eval[4], time_direct_eval[4]);
      eval_solution<10>(knots, comsol, 
                    psol_u1_kernel, psol_u2_kernel,
                    psol_v1_kernel, psol_v2_kernel,
                    psol_p1_kernel, psol_p2_kernel,
                    rms_error_u[5], rms_error_v[5], rms_error_p[5], 
                    rms_diff_u[5], rms_diff_v[5], rms_diff_p[5], 
                    time_setup[5], time_eval[5], time_direct_eval[5]);




      

      file << std::setw(15) << nbucket 
           << std::setw(15) << rms_diff_u[0]
           << std::setw(15) << rms_diff_u[1]
           << std::setw(15) << rms_diff_u[2]
           << std::setw(15) << rms_diff_u[3]
           << std::setw(15) << rms_diff_u[4]
           << std::setw(15) << rms_diff_u[5]
           << std::setw(15) << rms_diff_v[0]
           << std::setw(15) << rms_diff_v[1]
           << std::setw(15) << rms_diff_v[2]
           << std::setw(15) << rms_diff_v[3]
           << std::setw(15) << rms_diff_v[4]
           << std::setw(15) << rms_diff_v[5]
           << std::setw(15) << rms_diff_p[0]
           << std::setw(15) << rms_diff_p[1]
           << std::setw(15) << rms_diff_p[2]
           << std::setw(15) << rms_diff_p[3]
           << std::setw(15) << rms_diff_p[4]
           << std::setw(15) << rms_diff_p[5]
           << std::setw(15) << rms_error_u[0]
           << std::setw(15) << rms_error_u[1]
           << std::setw(15) << rms_error_u[2]
           << std::setw(15) << rms_error_u[3]
           << std::setw(15) << rms_error_u[4]
           << std::setw(15) << rms_error_u[5]
           << std::setw(15) << rms_error_v[0]
           << std::setw(15) << rms_error_v[1]
           << std::setw(15) << rms_error_v[2]
           << std::setw(15) << rms_error_v[3]
           << std::setw(15) << rms_error_v[4]
           << std::setw(15) << rms_error_v[5]
           << std::setw(15) << rms_error_p[0]
           << std::setw(15) << rms_error_p[1]
           << std::setw(15) << rms_error_p[2]
           << std::setw(15) << rms_error_p[3]
           << std::setw(15) << rms_error_p[4]
           << std::setw(15) << rms_error_p[5]
           << std::setw(15) << time_setup[0]
           << std::setw(15) << time_setup[1]
           << std::setw(15) << time_setup[2]
           << std::setw(15) << time_setup[3]
           << std::setw(15) << time_setup[4]
           << std::setw(15) << time_setup[5]
           << std::setw(15) << time_eval[0]
           << std::setw(15) << time_eval[1]
           << std::setw(15) << time_eval[2]
           << std::setw(15) << time_eval[3]
           << std::setw(15) << time_eval[4]
           << std::setw(15) << time_eval[5]
           << std::setw(15) << time_direct_eval[0]
           << std::setw(15) << time_direct_eval[1]
           << std::setw(15) << time_direct_eval[2]
           << std::setw(15) << time_direct_eval[3]
           << std::setw(15) << time_direct_eval[4]
           << std::setw(15) << time_direct_eval[5]

           << std::endl;

    }
    file.close();
}
