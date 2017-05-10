#include "filter.h"
#include "solve_stokes_MAPS.h"
#include "setup_knots.h"

int main(int argc, char **argv) {

    unsigned int nout,max_iter_linear,restart_linear,nx;
    int fibre_number,seed;
    double fibre_resolution,fibre_radius,particle_rate,react_rate,D;
    double dt_aim,k,gamma,rf,c0,epsilon_strength,epsilon_falloff;
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
        ("particle_rate", po::value<double>(&particle_rate)->default_value(1000.0), "particle rate")
        ("periodic", po::value<bool>(&periodic)->default_value(false), "periodic in x")
        ("react_rate", po::value<double>(&react_rate)->default_value(0.5), "particle reaction rate")
        ("epsilon_strength", po::value<double>(&epsilon_strength)->default_value(2.0), "boundary clustering fall-off")
        ("epsilon_falloff", po::value<double>(&epsilon_falloff)->default_value(0.3), "boundary clustering fall-off")

        ("c0", po::value<double>(&c0)->default_value(0.1), "kernel constant")
        ("nx", po::value<unsigned int>(&nx)->default_value(20), "nx")
        ("fibre_resolution", po::value<double>(&fibre_resolution)->default_value(0.75), "number of knots around each fibre")
        ("seed", po::value<int>(&seed)->default_value(10), "seed")
        ("fibre_number", po::value<int>(&fibre_number)->default_value(3), "number of fibres")
        ("fibre_radius", po::value<double>(&fibre_radius)->default_value(0.05), "radius of fibres")
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

    const double up_mult = 1.5;
    const double c2 = std::pow(c0,2);
    const double c3 = std::pow(c0,3);
    const double c4 = std::pow(c0,4);
    const double mu = 1.0;
    const double Tf = 2.0;
    const double L = 1.0;
    const double delta = L/nx;
    const int ny = up_mult*L/delta;
    const double deltay = up_mult*L/ny;
    const double boundary_layer = delta/5;
    const double s = 1.1*delta;
    const int timesteps = Tf/dt_aim;
    const double dt = Tf/timesteps;
    const double2 domain_min(0,0);
    const double2 domain_max(L,up_mult*L+1e-10);
    const double2 ns_buffer(L/2,1e-10);

    std::poisson_distribution<int> poisson(particle_rate*dt);
    std::uniform_real_distribution<double> uniform(L/10.0,L-L/10.0);
    std::default_random_engine generator(seed);

    // particles
    {
        particles.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,bool2(false));

        std::cout << "added "<<particles.size()<<" particles"<<std::endl;
    }

    // fibres
    {
        typename ParticlesType::value_type p;
        for (int ii=0; ii<fibre_number; ++ii) {
            for (int jj=0; jj<fibre_number; ++jj) {
                const double dx = L/fibre_number;
                const double2 origin = double2(
                                            (ii+0.5)*dx,
                                            (jj+0.5)*dx+0.5*(up_mult-1)*L
                                            );
                get<position>(p) = origin;
                fibres.push_back(p);
            }
        }
        fibres.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,bool2(false));

        std::cout << "added "<<fibres.size()<<" fibres"<<std::endl;
    }

    setup_knots(knots, fibres, fibre_radius, fibre_resolution, nx, domain_min, domain_max, c0, k,epsilon_strength,epsilon_falloff,periodic);

    solve_stokes_MAPS(knots,max_iter_linear,restart_linear,solver_in,c0);

    typedef typename position::value_type const & const_position_reference;
    typedef typename KnotsType::const_reference const_knot_reference;
    typedef typename KnotsType::reference knot_reference;
    typedef typename ParticlesType::const_reference const_particles_reference;
    typedef typename ParticlesType::reference particles_reference;

    auto psol_u1_op = create_dense_operator(particles,knots,
            [](const_position_reference dx,
               const_particles_reference a,
               const_knot_reference b) {
            return psol_u1(dx,get<kernel_constant>(b));
            });

      auto psol_v1_op = create_dense_operator(particles,knots,
            [](const_position_reference dx,
               const_particles_reference a,
               const_knot_reference b) {
            return psol_v1(dx,get<kernel_constant>(b));
            });

      auto psol_p1_op = create_dense_operator(particles,knots,
            [](const_position_reference dx,
               const_particles_reference a,
               const_knot_reference b) {
            return psol_p1(dx,get<kernel_constant>(b));
            });

      auto psol_u2_op = create_dense_operator(particles,knots,
            [](const_position_reference dx,
               const_particles_reference a,
               const_knot_reference b) {
            return psol_u2(dx,get<kernel_constant>(b));
            });

      auto psol_v2_op = create_dense_operator(particles,knots,
            [](const_position_reference dx,
               const_particles_reference a,
               const_knot_reference b) {
            return psol_v2(dx,get<kernel_constant>(b));
            });

      auto psol_p2_op = create_dense_operator(particles,knots,
            [](const_position_reference dx,
               const_particles_reference a,
               const_knot_reference b) {
            return psol_p2(dx,get<kernel_constant>(b));
            });

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
    Accumulate<std::plus<double> > sum;
    AccumulateWithinDistance<std::plus<double2> > sumv(fibre_radius);
    Accumulate<std::plus<double3> > sumv3;
    Accumulate<Aboria::max<double> > max;
    max.set_init(0);
    AccumulateWithinDistance<std::bit_or<bool> > any(fibre_radius);
    any.set_init(false);
    VectorSymbolic<double,2> vector;
    VectorSymbolic<double,3> vector3;
    VectorSymbolic<double,9> matrix;
    Normal N;
    Uniform U;


    std::cout << "starting timesteps!"<<std::endl;
    const int timesteps_per_out = timesteps/nout;
    for (int ii=0; ii<timesteps; ++ii) {
        // add new particles
        const int new_n = poisson(generator);
        for (int jj=0; jj<new_n; ++jj) {
            ParticlesType::value_type p;
            get<position>(p) = double2(uniform(generator),up_mult*L);
            get<kernel_constant>(p) = c0;
            particles.push_back(p);
        }


        // diffusion with drift
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
        typedef Eigen::Map<vector_type> map_type;
        map_type(get<velocity_u>(particles).data(),particles.size()) = 
            psol_u1_op*map_type(get<alpha1>(knots).data(),knots.size())+
            psol_u2_op*map_type(get<alpha2>(knots).data(),knots.size());
        map_type(get<velocity_v>(particles).data(),particles.size()) = 
            psol_v1_op*map_type(get<alpha1>(knots).data(),knots.size())+
            psol_v2_op*map_type(get<alpha2>(knots).data(),knots.size());

        r[a] += std::sqrt(2.0*D*dt)*vector(N[a],N[a]) + dt*vector(vu[a],vv[a]);

        if (ii % timesteps_per_out == 0) {
            std::cout << "timestep "<<ii<<" of "<<timesteps<<std::endl;
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

    }

}
