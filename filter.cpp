#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40
#undef NDEBUG
#include "Aboria.h"
using namespace Aboria;

#include <boost/math/constants/constants.hpp>
#include "boost/program_options.hpp" 
namespace po = boost::program_options;
#include <chrono>
typedef std::chrono::system_clock Clock;

const double PI = boost::math::constants::pi<double>();
enum linear_solver {CG, BiCGSTAB, GMRES, DGMRES};

template<typename Kernel,typename VectorType>
void solve(Kernel &&kernel, VectorType &&result, VectorType &&source, size_t max_iter=10, size_t restart=10, linear_solver solver=CG) {
    typedef typename std::remove_reference<Kernel>::type::Scalar Scalar;
    switch (solver) {
        case CG: {
            Eigen::ConjugateGradient<
                typename std::remove_reference<Kernel>::type, 
                         Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<Scalar>> cg;
            cg.setMaxIterations(max_iter);
            cg.compute(kernel);
            result = cg.solveWithGuess(source,result);
            std::cout << "CG:    #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << " actual error: "<<(kernel*result - source).norm() / source.norm() << std::endl;
            break;
                 }
        case BiCGSTAB: {
            Eigen::BiCGSTAB<
                typename std::remove_reference<Kernel>::type, 
                     Eigen::DiagonalPreconditioner<Scalar>> bicg;
            bicg.setMaxIterations(max_iter);
            bicg.compute(kernel);
            result = bicg.solveWithGuess(source,result);
            //result = bicg.solve(source);
            std::cout << "BiCGSTAB:    #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << " actual error: "<<(kernel*result - source).norm() / source.norm() << std::endl;
            break;
               }
        case GMRES: {
            Eigen::GMRES<
                typename std::remove_reference<Kernel>::type, 
                    Eigen::DiagonalPreconditioner<Scalar>> gmres;
            gmres.set_restart(restart);
            gmres.setMaxIterations(max_iter);
            gmres.compute(kernel);
            result = gmres.solveWithGuess(source,result);
            std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << " actual error: "<<(kernel*result - source).norm() / source.norm()<<std::endl;
            break;
                    }
        case DGMRES: {
            Eigen::DGMRES<
                typename std::remove_reference<Kernel>::type, 
                    Eigen::DiagonalPreconditioner<Scalar>> dgmres;
            dgmres.set_restart(restart);
            dgmres.setMaxIterations(max_iter);
            dgmres.compute(kernel);
            result = dgmres.solveWithGuess(source,result);
            std::cout << "DGMRES:    #iterations: " << dgmres.iterations() << ", estimated error: " << dgmres.error() << " actual error: "<<(kernel*result - source).norm() / source.norm()<<std::endl;
            break;
                    }


    }


}

int main(int argc, char **argv) {

    unsigned int nout,max_iter_linear,restart_linear,nx;
    int fibre_resolution,fibre_number,seed;
    double fibre_radius,particle_rate,react_rate,D;
    double dt_aim,h0_factor,k,gamma,rf,c0,epsilon;
    unsigned int solver_in;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("max_iter_linear", po::value<unsigned int>(&max_iter_linear)->default_value(100), "maximum iterations for linear solve")
        ("restart_linear", po::value<unsigned int>(&restart_linear)->default_value(101), "iterations until restart for linear solve")
        ("linear_solver", po::value<unsigned int>(&solver_in)->default_value(2), "linear solver")
        ("nout", po::value<unsigned int>(&nout)->default_value(10), "number of output points")
        ("k", po::value<double>(&k)->default_value(0.1), "spring constant")
        ("D", po::value<double>(&D)->default_value(1.0), "diffusion constant")
        ("particle_rate", po::value<double>(&particle_rate)->default_value(1000.0), "particle rate")
        ("react_rate", po::value<double>(&react_rate)->default_value(0.5), "particle reaction rate")
        ("epsilon", po::value<double>(&epsilon)->default_value(10.0), "boundary clustering fall-off")
        ("c0", po::value<double>(&c0)->default_value(0.01), "kernel constant")
        ("nx", po::value<unsigned int>(&nx)->default_value(19), "nx")
        ("fibre_resolution", po::value<int>(&fibre_resolution)->default_value(10), "number of knots around each fibre")
        ("seed", po::value<int>(&seed)->default_value(10), "seed")
        ("fibre_number", po::value<int>(&fibre_number)->default_value(3), "number of fibres")
        ("fibre_radius", po::value<double>(&fibre_radius)->default_value(0.1), "radius of fibres")
        ("h0_factor", po::value<double>(&h0_factor)->default_value(4.0), "h0 factor")
        ("dt", po::value<double>(&dt_aim)->default_value(0.1), "timestep")
    ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);  

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    

    ABORIA_VARIABLE(inlet,uint8_t,"is inlet knot")
    ABORIA_VARIABLE(outlet,uint8_t,"is outlet knot")
    ABORIA_VARIABLE(boundary,uint8_t,"is boundary knot")
    ABORIA_VARIABLE(interior,uint8_t,"is interior knot")

    ABORIA_VARIABLE(velocity_u,double,"velocity u");
    ABORIA_VARIABLE(velocity_v,double,"velocity v");
    ABORIA_VARIABLE(pressure,double,"pressure");
    ABORIA_VARIABLE(alpha,double2,"alpha");
    
    ABORIA_VARIABLE(boundary_normal,double2,"unit normal vector to boundary")
    ABORIA_VARIABLE(kernel_constant,double,"kernel constant")

    typedef Particles<std::tuple<alpha,boundary,inlet,outlet,interior,velocity_u,velocity_v,pressure,kernel_constant,boundary_normal>,2> KnotsType;
    typedef Particles<std::tuple<kernel_constant>,2> ParticlesType;
    typedef position_d<2> position;
    KnotsType knots;
    ParticlesType particles;
    ParticlesType fibres;

    const double c2 = std::pow(c0,2);
    const double c3 = std::pow(c0,3);
    const double c4 = std::pow(c0,4);
    const double mu = 1.0;
    const double flow_rate = 1.0; 
    const double Tf = 2.0;
    const double L = 1.0;
    const double delta = L/nx;
    const double boundary_layer = delta/5;
    const double s = 1.1*delta;
    const double h0 = h0_factor*delta;
    const int timesteps = Tf/dt_aim;
    const double dt = Tf/timesteps;
    const double dt_adapt = (1.0/100.0)*PI/sqrt(2*k);
    const double2 domain_min(-L/2,0);
    const double2 domain_max(1.5*L,L+1e-10);

    std::poisson_distribution<int> poisson(particle_rate*dt);
    std::uniform_real_distribution<double> uniform(L/10.0,L-L/10.0);
    std::default_random_engine generator(seed);

    // particles 
    {
        particles.init_neighbour_search(domain_min,domain_max,fibre_radius+boundary_layer,bool2(false));

        std::cout << "added "<<particles.size()<<" particles"<<std::endl;
    }

    // fibres
    {
        typename ParticlesType::value_type p;
        for (int ii=0; ii<fibre_number; ++ii) {
            for (int jj=0; jj<fibre_number; ++jj) {
                const double dx = L/fibre_number;
                const double2 origin = double2((ii+0.5)*dx,(jj+0.5)*dx);
                get<position>(p) = origin;
                fibres.push_back(p);
            }
        }
        fibres.init_neighbour_search(domain_min,domain_max,fibre_radius+boundary_layer,bool2(false));

        std::cout << "added "<<fibres.size()<<" fibres"<<std::endl;
    }
 
    // knots
    {

        typename KnotsType::value_type p;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> uniform(0+delta/10.0,L-delta/10.0);
        const double volume = L*L-fibres.size()*PI*std::pow(fibre_radius,2);
        const int N = nx*nx*volume; 
        while (knots.size() < N) {
            get<position>(p) = double2(uniform(generator),uniform(generator));
            bool outside_fibres = true;
            for (auto tpl: fibres.get_neighbours(get<position>(p))) {
                if ((get<position>(p)-get<position>(std::get<0>(tpl))).norm() 
                        < fibre_radius) {
                    outside_fibres = false;
                }
            }
            if (outside_fibres) {
                get<boundary>(p) = false;
                get<inlet>(p) = false;
                get<outlet>(p) = false;
                get<interior>(p) = true;
                knots.push_back(p);
            }

        }
        std::cout << "added "<<knots.size()<<" interior knots"<<std::endl;

        // fibre boundaries
        for (int ii=0; ii<fibres.size(); ++ii) {
            const double2 origin = get<position>(fibres)[ii];
            const double dtheta = 2*PI/fibre_resolution;
            for (int kk=0; kk<fibre_resolution; ++kk) {
                get<position>(p) = origin + fibre_radius*double2(cos(kk*dtheta),sin(kk*dtheta));
                //get<boundary_tangent>(p) = double2(-sin(kk*dtheta),cos(kk*dtheta));
                get<boundary_normal>(p) = double2(cos(kk*dtheta),sin(kk*dtheta));
                get<boundary>(p) = true;
                get<inlet>(p) = false;
                get<outlet>(p) = false;
                get<interior>(p) = false;
                knots.push_back(p);
            }
        }

        const double dx = L/nx;
        for (int ii=0; ii<=nx; ++ii) {
            if ((ii > 0) && (ii<nx)) {
                // inlet
                get<position>(p) = double2(ii*dx,L);
                get<boundary>(p) = false;
                get<inlet>(p) = true;
                get<outlet>(p) = false;
                get<interior>(p) = false;
                get<boundary_normal>(p) = double2(0,-1);
                knots.push_back(p);

                // outlet 
                get<position>(p) = double2(ii*dx,0.0);
                get<boundary>(p) = false;
                get<inlet>(p) = false;
                get<outlet>(p) = true;
                get<interior>(p) = false;
                get<boundary_normal>(p) = double2(0,1);
                knots.push_back(p);
            }

            // boundary - left 
            get<position>(p) = double2(0.0,ii*dx);
            get<boundary>(p) = true;
            get<inlet>(p) = false;
            get<interior>(p) = false;
            get<outlet>(p) = false;
            get<boundary_normal>(p) = double2(1,0);
            knots.push_back(p);

            // boundary - right 
            get<position>(p) = double2(L,ii*dx);
            get<boundary>(p) = true;
            get<inlet>(p) = false;
            get<interior>(p) = false;
            get<outlet>(p) = false;
            get<boundary_normal>(p) = double2(-1,0);
            knots.push_back(p);
        }

        knots.init_neighbour_search(domain_min,domain_max,2*h0,bool2(false));
        std::cout << "added "<<knots.size()<<" knots with c0 = " <<c0<< std::endl;
    }
    
    vtkWriteGrid("init_knots",0,knots.get_grid(true));

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
    Symbol<alpha> al;

    Symbol<boundary_normal> normal;

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
    Accumulate<std::plus<double2> > sumv;
    Accumulate<std::plus<double3> > sumv3;
    Accumulate<Aboria::max<double> > max;
    max.set_init(0);
    Accumulate<std::bit_or<bool> > any;
    any.set_init(false);
    VectorSymbolic<double,2> vector;      
    VectorSymbolic<double,3> vector3;      
    VectorSymbolic<double,9> matrix;      
    Normal N;
    Uniform U;

    c[i] = c0;
    c[a] = c0;

    auto kernel_mq = deep_copy(
            sqrt(dot(dx,dx)+pow(c[i],2))
        );

    auto phi_sol = deep_copy(
            ((1.0/75.0)*sqrt(dot(dx,dx)+pow(c[i],2))*(4.0*pow(dot(dx,dx),2)+48.0*dot(dx,dx)*pow(c[i],2)-61.0*pow(c[i],4)) - pow(c[i],3)*log(c[i])*dot(dx,dx) - (1.0/5.0)*(5.0*dot(dx,dx)-2*pow(c[i],2))*pow(c[i],3)*log(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))/(12.0*mu)
            );

    auto phi_sol_dash_div_r = deep_copy(
            (4.0*pow(dot(dx,dx),3) + 36.0*pow(dot(dx,dx),2)*pow(c[i],2) + 39.0*dot(dx,dx)*pow(c[i],4) + 7.0*pow(c[i],6))/(180.0*pow(dot(dx,dx)+pow(c[i],2),1.5))
            - (5.0*pow(dot(dx,dx),2) + 3.0*dot(dx,dx)*pow(c[i],2) - 2.0*pow(c[i],4))*pow(c[i],3)/(60.0*pow(dot(dx,dx)+pow(c[i],2),1.5)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            - (1.0/6.0)*pow(c[i],3)*log(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))) - (1.0/6.0)*pow(c[i],3)*log(c[i])
            );

    

    auto phi_sol_dash_dash = deep_copy(
            (16.0*pow(dot(dx,dx),3)+84.0*pow(dot(dx,dx),2)*pow(c[i],2)+96.0*dot(dx,dx)*pow(c[i],4)+7*pow(c[i],6))/(180.0*pow(dot(dx,dx)+pow(c[i],2),1.5))
            - (20.0*pow(dot(dx,dx),2)+25.0*dot(dx,dx)*pow(c[i],2)-2.0*pow(c[i],4))*pow(c[i],3)/(60.0*pow(dot(dx,dx)+pow(c[i],2),1.5)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            + (5.0*dot(dx,dx)-2*pow(c[i],2))*pow(c[i],3)*dot(dx,dx)/(60.0*(dot(dx,dx)+pow(c[i],2))*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            - (1.0/6.0)*pow(c[i],3)*log(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))) - (1.0/6.0)*pow(c[i],3)*log(c[i])
            );
 
    auto phi_sol_dash_dash_dash = deep_copy(
            (76.0*pow(dot(dx,dx),2)+176.0*dot(dx,dx)*pow(c[i],2)+285.0*pow(c[i],4))*norm(dx)/(300.0*pow(dot(dx,dx)+pow(c[i],2),1.5))
            + (4.0*pow(dot(dx,dx),2)+48.0*dot(dx,dx)*pow(c[i],2)-61.0*pow(c[i],4))*pow(norm(dx),3)/(300.0*pow(dot(dx,dx)+pow(c[i],2),5.0/2.0))
            + (10.0*pow(dot(dx,dx),2)+15.0*dot(dx,dx)*pow(c[i],2)-2.0*pow(c[i],4))*pow(c[i],3)*norm(dx)/(20.0*pow(dot(dx,dx)+pow(c[i],2),2)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            - (5.0*dot(dx,dx)+22.0*pow(c[i],2))*pow(c[i],3)*norm(dx)/(20.0*pow(dot(dx,dx)+pow(c[i],2),1.5)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            + (-5.0*dot(dx,dx)+2*pow(c[i],2))*pow(c[i],3)*pow(norm(dx),3)/(20.0*pow(dot(dx,dx)+pow(c[i],2),5.0/2.0)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            + (-5.0*dot(dx,dx)+2*pow(c[i],2))*pow(c[i],3)*pow(norm(dx),3)/(30.0*pow(dot(dx,dx)+pow(c[i],2),1.5)*pow(c[i]+sqrt(dot(dx,dx)+pow(c[i],2)),3))
            );

    auto kernel_mq_test = deep_copy(
            sqrt(dot(dx,dx)+c2)
        );

    auto phi_sol_test = deep_copy(
            ((1.0/75.0)*sqrt(dot(dx,dx)+c2)*(4.0*pow(dot(dx,dx),2)+48.0*dot(dx,dx)*c2-61.0*c2*c2) - c3*log(c0)*dot(dx,dx) - (1.0/5.0)*(5.0*dot(dx,dx)-2*c2)*c3*log(c0+sqrt(dot(dx,dx)+c2)))/(12.0*mu)
            );

    auto phi_sol_dash_div_r_test = deep_copy(
            (4.0*pow(dot(dx,dx),3) + 36.0*pow(dot(dx,dx),2)*c2 + 39.0*dot(dx,dx)*c2*c2 + 7.0*c3*c3)/(180.0*pow(dot(dx,dx)+c2,1.5))
            - (5.0*pow(dot(dx,dx),2) + 3.0*dot(dx,dx)*c2 - 2*c2*c2)*c3/(60.0*pow(dot(dx,dx)+c2,1.5)*(c0+sqrt(dot(dx,dx)+c2)))
            - (1.0/6.0)*c3*log(c0+sqrt(dot(dx,dx)+c2)) - (1.0/6.0)*c3*log(c0)
            );

    auto phi_sol_dash_dash_test = deep_copy(
            (16.0*pow(dot(dx,dx),3)+84.0*pow(dot(dx,dx),2)*c2+96.0*dot(dx,dx)*c2*c2+7*c3*c3)/(180.0*pow(dot(dx,dx)+c2,1.5))
            - (20.0*pow(dot(dx,dx),2)+25.0*dot(dx,dx)*c2-2.0*c2*c2)*c3/(60.0*pow(dot(dx,dx)+c2,1.5)*(c0+sqrt(dot(dx,dx)+c2)))
            + (5.0*dot(dx,dx)-2.0*c2)*c3*dot(dx,dx)/(60.0*(dot(dx,dx)+c2)*(c0+sqrt(dot(dx,dx)+c2)))
            - (1.0/6.0)*c3*log(c0+sqrt(dot(dx,dx)+c2)) - (1.0/6.0)*c3*log(c0)
            );
 
    auto phi_sol_dash_dash_dash_test = deep_copy(
            (76.0*pow(dot(dx,dx),2)+176.0*dot(dx,dx)*c2+285.0*c2*c2)*norm(dx)/(300.0*pow(dot(dx,dx)+c2,1.5))
            + (4.0*pow(dot(dx,dx),2)+48.0*dot(dx,dx)*c2-61.0*c2*c2)*pow(norm(dx),3)/(300.0*pow(dot(dx,dx)+c2,5.0/2.0))
            + (10.0*pow(dot(dx,dx),2)+15.0*dot(dx,dx)*c2-2*c2*c2)*c3*norm(dx)/(20.0*pow(dot(dx,dx)+c2,2)*(c0+sqrt(dot(dx,dx)+c2)))
            - (5.0*dot(dx,dx)+22.0*c2)*c3*norm(dx)/(20.0*pow(dot(dx,dx)+c2,1.5)*(c0+sqrt(dot(dx,dx)+c2)))
            + (-5.0*dot(dx,dx)+2*c2)*c3*pow(norm(dx),3)/(20.0*pow(dot(dx,dx)+c2,5.0/2.0)*(c0+sqrt(dot(dx,dx)+c2)))
            + (-5.0*dot(dx,dx)+2.0*c2)*c3*pow(norm(dx),3)/(30.0*pow(dot(dx,dx)+c2,1.5)*pow(c0+sqrt(dot(dx,dx)+c2),3))
            );

    // i = 1, l = 1
    auto psol_u1_test = deep_copy(
            if_else(norm(dx)==0
            ,(1.0/mu)*(phi_sol_dash_div_r_test + phi_sol_dash_dash_test)
            ,(1.0/mu)*(phi_sol_dash_div_r_test*dx[0]*dx[0] + phi_sol_dash_dash_test*dx[1]*dx[1])/dot(dx,dx)
            )
            );

    auto psol_p1_test = deep_copy(
            if_else(norm(dx)==0
            ,-(phi_sol_dash_dash_dash_test)
            ,(phi_sol_dash_dash_dash_test + (phi_sol_dash_dash_test - phi_sol_dash_div_r_test)/norm(dx))*dx[0]/norm(dx)
            )
            );




    // i = 1, l = 1
    auto psol_u1 = deep_copy(
            if_else(norm(dx)==0
            ,(1.0/mu)*(phi_sol_dash_div_r + phi_sol_dash_dash)
            ,(1.0/mu)*(phi_sol_dash_div_r*dx[0]*dx[0] + phi_sol_dash_dash*dx[1]*dx[1])/dot(dx,dx)
            )
            );

    // i = 1, l = 2
    auto psol_u2 = deep_copy(
            if_else(norm(dx)==0
            ,(1.0/mu)*(phi_sol_dash_div_r - phi_sol_dash_dash)
            ,(1.0/mu)*(phi_sol_dash_div_r - phi_sol_dash_dash)*(dx[0]*dx[1])/dot(dx,dx)
            )
            );

    // i = 2, l = 1
    auto psol_v1 = deep_copy(
            if_else(norm(dx)==0
            ,(1.0/mu)*(phi_sol_dash_div_r - phi_sol_dash_dash)
            ,(1.0/mu)*(phi_sol_dash_div_r - phi_sol_dash_dash)*(dx[0]*dx[1])/dot(dx,dx)
            )
            );

    // i = 2, l = 2
    auto psol_v2 = deep_copy(
            if_else(norm(dx)==0
            ,(1.0/mu)*(phi_sol_dash_div_r+phi_sol_dash_dash)
            ,(1.0/mu)*(phi_sol_dash_div_r*dx[1]*dx[1] + phi_sol_dash_dash*dx[0]*dx[0])/dot(dx,dx)
            )
            );


    auto psol_p1 = deep_copy(
            if_else(norm(dx)==0
            ,-(phi_sol_dash_dash_dash)
            ,(phi_sol_dash_dash_dash + (phi_sol_dash_dash - phi_sol_dash_div_r)/norm(dx))*dx[0]/norm(dx)
            )
            );
    auto psol_p2 = deep_copy(
            if_else(norm(dx)==0
            ,-(phi_sol_dash_dash_dash)
            ,(phi_sol_dash_dash_dash + (phi_sol_dash_dash - phi_sol_dash_div_r)/norm(dx))*dx[1]/norm(dx)
            )
            );


    auto kernel_mq_pk = deep_copy(
            sqrt(dot(dpk,dpk)+pow(c[a],2))
        );

    auto phi_sol_pk = deep_copy(
            ((1.0/75.0)*sqrt(dot(dpk,dpk)+pow(c[a],2))*(4.0*pow(dot(dpk,dpk),2)+48.0*dot(dpk,dpk)*pow(c[a],2)-61.0*pow(c[a],4)) - pow(c[a],3)*log(c[a])*dot(dpk,dpk) - (1.0/5.0)*(5.0*dot(dpk,dpk)-2*pow(c[a],2))*pow(c[a],3)*log(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2))))/(12.0*mu)
            );

    auto phi_sol_dash_div_r_pk = deep_copy(
            (4.0*pow(dot(dpk,dpk),3) + 36.0*pow(dot(dpk,dpk),2)*pow(c[a],2) + 39.0*dot(dpk,dpk)*pow(c[a],4) + 7*pow(c[a],6))/(180.0*pow(dot(dpk,dpk)+pow(c[a],2),1.5))
            - (5.0*pow(dot(dpk,dpk),2) + 3*dot(dpk,dpk)*pow(c[a],2) - 2*pow(c[a],4))*pow(c[a],3)/(60.0*pow(dot(dpk,dpk)+pow(c[a],2),1.5)*(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2))))
            - (1.0/6.0)*pow(c[a],3)*log(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2))) - (1.0/6.0)*pow(c[a],3)*log(c[a])
            );

    auto phi_sol_dash_dash_pk = deep_copy(
            (16.0*pow(dot(dpk,dpk),3)+84*pow(dot(dpk,dpk),2)*pow(c[a],2)+96*dot(dpk,dpk)*pow(c[a],4)+7*pow(c[a],6))/(180.0*pow(dot(dpk,dpk)+pow(c[a],2),1.5))
            - (20*pow(dot(dpk,dpk),2)+25*dot(dpk,dpk)*pow(c[a],2)-2*pow(c[a],4))*pow(c[a],3)/(60.0*pow(dot(dpk,dpk)+pow(c[a],2),1.5)*(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2))))
            + (5*dot(dpk,dpk)-2*pow(c[a],2))*pow(c[a],3)*dot(dpk,dpk)/(60*(dot(dpk,dpk)+pow(c[a],2))*(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2))))
            - (1.0/6.0)*pow(c[a],3)*log(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2))) - (1.0/6.0)*pow(c[a],3)*log(c[a])
            );
 
    auto phi_sol_dash_dash_dash_pk = deep_copy(
            (76.0*pow(dot(dpk,dpk),2)+176*dot(dpk,dpk)*pow(c[a],2)+285*pow(c[a],4))*norm(dpk)/(300*pow(dot(dpk,dpk)+pow(c[a],2),1.5))
            + (4*pow(dot(dpk,dpk),2)+48*dot(dpk,dpk)*pow(c[a],2)-61*pow(c[a],4))*pow(norm(dpk),3)/(300*pow(dot(dpk,dpk)+pow(c[a],2),5.0/2.0))
            + (10*pow(dot(dpk,dpk),2)+15*dot(dpk,dpk)*pow(c[a],2)-2*pow(c[a],4))*pow(c[a],3)*norm(dpk)/(20*pow(dot(dpk,dpk)+pow(c[a],2),2)*(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2))))
            - (5*dot(dpk,dpk)+22*pow(c[a],2))*pow(c[a],3)*norm(dpk)/(20*pow(dot(dpk,dpk)+pow(c[a],2),1.5)*(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2))))
            + (-5*dot(dpk,dpk)+2*pow(c[a],2))*pow(c[a],3)*pow(norm(dpk),3)/(20*pow(dot(dpk,dpk)+pow(c[a],2),5.0/2.0)*(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2))))
            + (-5*dot(dpk,dpk)+2*pow(c[a],2))*pow(c[a],3)*pow(norm(dpk),3)/(30*pow(dot(dpk,dpk)+pow(c[a],2),1.5)*pow(c[a]+sqrt(dot(dpk,dpk)+pow(c[a],2)),3))
            );


    // i = 1, l = 1
    auto psol_u1_pk = deep_copy(
            if_else(norm(dpk)==0
            ,(1.0/mu)*(phi_sol_dash_div_r_pk + phi_sol_dash_dash_pk)
            ,(1.0/mu)*(phi_sol_dash_div_r_pk*dpk[0]*dpk[0] + phi_sol_dash_dash_pk*dpk[1]*dpk[1])/dot(dpk,dpk)
            )
            );

    // i = 1, l = 2
    auto psol_u2_pk = deep_copy(
            if_else(norm(dpk)==0
            ,(1.0/mu)*(phi_sol_dash_div_r_pk - phi_sol_dash_dash_pk)
            ,(1.0/mu)*(phi_sol_dash_div_r_pk - phi_sol_dash_dash_pk)*(dpk[0]*dpk[1])/dot(dpk,dpk)
            )
            );

    // i = 2, l = 1
    auto psol_v1_pk = deep_copy(
            if_else(norm(dpk)==0
            ,(1.0/mu)*(phi_sol_dash_div_r_pk - phi_sol_dash_dash_pk)
            ,(1.0/mu)*(phi_sol_dash_div_r_pk - phi_sol_dash_dash_pk)*(dpk[0]*dpk[1])/dot(dpk,dpk)
            )
            );

    // i = 2, l = 2
    auto psol_v2_pk = deep_copy(
            if_else(norm(dpk)==0
            ,(1.0/mu)*(phi_sol_dash_div_r_pk+phi_sol_dash_dash_pk)
            ,(1.0/mu)*(phi_sol_dash_div_r_pk*dpk[1]*dpk[1] + phi_sol_dash_dash_pk*dpk[0]*dpk[0])/dot(dpk,dpk)
            )
            );

    auto psol_p1_pk = deep_copy(
            if_else(norm(dpk)==0
            ,-(phi_sol_dash_dash_dash_pk)
            ,(phi_sol_dash_dash_dash_pk + (phi_sol_dash_dash_pk - phi_sol_dash_div_r_pk)/norm(dpk))*dpk[0]/norm(dpk)
            )
            );
    auto psol_p2_pk = deep_copy(
            if_else(norm(dpk)==0
            ,-(phi_sol_dash_dash_dash_pk)
            ,(phi_sol_dash_dash_dash_pk + (phi_sol_dash_dash_pk - phi_sol_dash_div_r_pk)/norm(dpk))*dpk[1]/norm(dpk)
            )
            );



    
    auto spring_force_kk = deep_copy(
        if_else(dot(dx,dx)==0
            ,vector(0,0)
            ,(-k*(s-norm(dx))/norm(dx))*dx
            )
        );

    auto spring_force_kf = deep_copy(
        if_else(dot(dkf,dkf)==0
            ,vector(0,0)
            ,if_else(norm(dkf)>fibre_radius+boundary_layer
                ,0.0
                //,delta*k*exp(epsilon*(fibre_radius+boundary_layer-norm(dkf)))
                ,-10*k*(fibre_radius+boundary_layer-norm(dkf))
                )*dkf/norm(dkf)
            )
        );


    auto spring_force_kb = deep_copy(
        vector(
                if_else(r[i][0] < boundary_layer
                ,10*k*(boundary_layer-r[i][0])
                ,if_else(r[i][0] > L-boundary_layer
                    ,10*k*(L-boundary_layer-r[i][0])
                    ,0.0
                    //,-delta*k*(exp(epsilon*(boundary_layer-r[i][0]))
                    //                      - exp(epsilon*(r[i][0]-L+boundary_layer)))
                    )
                )
                ,if_else(r[i][1] < boundary_layer 
                ,10*k*(boundary_layer-r[i][1])
                ,if_else(r[i][1] > L-boundary_layer
                    ,10*k*(L-boundary_layer-r[i][1])
                    ,0.0
                    )
                )
            )
        );
            

    // adapt knot locations
    for (int ii=0; ii<1000; ii++) {
        r[i] += dt_adapt*if_else(is_i[i]
                    ,sumv(j,norm(dx)<s,spring_force_kk)
                        + sumv(bf,true,spring_force_kf)
                        + spring_force_kb
                    ,vector(0,0)
                );
    }

    vtkWriteGrid("init_knots",1,knots.get_grid(true));


    //B1: u = 0 at inlet and b, p = 0 at outlet
    auto A11 = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,kernel_mq
               ,if_else(is_out[i]
                   ,psol_p1
                   ,psol_u1
                   )
               )
            );
    auto A11_test = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,kernel_mq_test
               ,if_else(is_out[i]
                   ,psol_p1_test
                   ,psol_u1_test
                   )
               )
            );
    auto A11_ne = deep_copy(
            if_else(is_i[i]
               ,kernel_mq
               ,if_else(is_out[i]
                   ,psol_p1
                   ,psol_u1
                   )
               )
            );


    auto A11_test_ne = deep_copy(
            if_else(is_i[i]
               ,kernel_mq_test
               ,if_else(is_out[i]
                   ,psol_p1_test
                   ,psol_u1_test
                   )
               )
            );



    auto A12 = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,0.0
               ,if_else(is_out[i]
                   ,psol_p2
                   ,psol_u2
                   )
               )
            );

    //B2: v = -flow_rate at inlet and v=0 at b, u = 0 at outlet
    auto A21 = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,0.0
               ,if_else(is_out[i]
                   ,psol_u1
                   ,psol_v1
                   )
               )
            );

    auto A22 = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,kernel_mq
               ,if_else(is_out[i]
                   ,psol_u2
                   ,psol_v2
                   )
               )
            );
    
    auto A = create_block_eigen_operator<2,2>(
            A11,A12,
            A21,A22
            );

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
    typedef Eigen::Map<vector_type> map_type;
    map_type eigen_pressure(get<pressure>(knots).data(),knots.size());
    map_type eigen_velocity_u(get<velocity_u>(knots).data(),knots.size());
    map_type eigen_velocity_v(get<velocity_u>(knots).data(),knots.size());
    pr[i] = 1.0;
    eigen_velocity_u = A11*eigen_pressure;
    eigen_velocity_v = A11_test*eigen_pressure;
    al[i] = sumv(j,true,vector(A11_test_ne*pr[j],A11_ne*pr[j]));
    for (int ii=0;ii<knots.size();++ii) {
        const double difference = 3.0*eigen_velocity_u[ii]-eigen_velocity_v[ii]-get<alpha>(knots)[ii][0]-get<alpha>(knots)[ii][1];
        if (difference != 0) {
            std::cout << "velocity_u = "<<eigen_velocity_u[ii]<<
                         "velocity_v = "<<eigen_velocity_v[ii]<<
                         "alpha = "<<get<alpha>(knots)[ii][0]<<
                         "alpha = "<<get<alpha>(knots)[ii][1]<<
                         "difference= "<<difference<<std::endl;
        }
    }

    vector_type source(2*knots.size());
    vector_type alphas(2*knots.size());


    for (int ii=0; ii<knots.size(); ++ii) {
        source(ii) = 0.0;
        if (get<inlet>(knots)[ii]) {
            source(knots.size()+ii) = -flow_rate;
        } else {
            source(knots.size()+ii) = 0.0;
        }
        
        alphas[ii] = 0.0;
        alphas[knots.size()+ii] = 0.0;
    }
    std::cout << "good!"<<std::endl;

    /*
    solve(A,alphas,source,
                max_iter_linear,restart_linear,(linear_solver)solver_in);
                */

    for (int ii=0; ii<knots.size(); ++ii) {
        get<alpha>(knots)[ii][0] = alphas[ii];
        get<alpha>(knots)[ii][1] = alphas[ii+knots.size()];
    }

    vu[i] = sum(j,true,psol_u1*al[j][0] + psol_u2*al[j][1]);
    vv[i] = sum(j,true,psol_v1*al[j][0] + psol_v2*al[j][1]);
    pr[i] = sum(j,true,psol_p1*al[j][0] + psol_p2*al[j][1]);

    vtkWriteGrid("MAPS",0,knots.get_grid(true));

    for (int ii=0; ii<timesteps; ++ii) {
        // add new particles
        const int new_n = poisson(generator);
        for (int jj=0; jj<new_n; ++jj) {
            ParticlesType::value_type p;
            get<position>(p) = double2(uniform(generator),L);
            get<kernel_constant>(p) = c0;
            particles.push_back(p);
        }
        std::cout << "finished adding particles"<<std::endl;


        // diffusion with drift
        r[a] += std::sqrt(2.0*D*dt)*vector(N,N);

        vtkWriteGrid("particles",ii,particles.get_grid(true));
            
        /*
            + dt*vector(
                    sum(j,true,psol_u1_pk*al[j][0] + psol_u2_pk*al[j][1]),
                    sum(j,true,psol_v1_pk*al[j][0] + psol_v2_pk*al[j][1])
                    );
                    */

        // react with fibres
        alive_[a] = !any(bf,norm(dpf) < fibre_radius,U<react_rate); 

        // react with side walls
        alive_[a] = !if_else(r[a][0] > L
                           ,U<react_rate
                           ,if_else(r[a][0] < 0
                               ,U<react_rate
                               ,false
                               )
                           );
        
        // reflect off fibres (if still alive) 
        r[a] += sumv(bf, norm(dpf)<fibre_radius
                        ,(fibre_radius/norm(dpf)-1)*dpf);
            
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

