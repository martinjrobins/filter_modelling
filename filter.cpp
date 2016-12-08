#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40
#include "Aboria.h"
using namespace Aboria;

#include <boost/math/constants/constants.hpp>
#include "boost/program_options.hpp" 
namespace po = boost::program_options;
#include <chrono>
typedef std::chrono::system_clock Clock;

const double PI = boost::math::constants::pi<double>();
enum linear_solver {CG, BiCGSTAB, GMRES};

template<typename Kernel,typename VectorType>
void solve(Kernel &&kernel, VectorType &&result, VectorType &&source, size_t max_iter=10, size_t restart=10, linear_solver solver=CG) {
    switch (solver) {
        case CG: {
            Eigen::ConjugateGradient<
                typename std::remove_reference<Kernel>::type, 
                         Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
            cg.setMaxIterations(max_iter);
            cg.compute(kernel);
            result = cg.solveWithGuess(source,result);
            std::cout << "CG:    #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
            break;
                 }
        case BiCGSTAB: {
            Eigen::BiCGSTAB<
                typename std::remove_reference<Kernel>::type, 
                     Eigen::DiagonalPreconditioner<double>> bicg;
            bicg.setMaxIterations(max_iter);
            bicg.compute(kernel);
            result = bicg.solveWithGuess(source,result);
            //result = bicg.solve(source);
            std::cout << "BiCGSTAB:    #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
            break;
               }
        case GMRES: {
            Eigen::GMRES<
                typename std::remove_reference<Kernel>::type, 
                    Eigen::DiagonalPreconditioner<double>> gmres;
            gmres.set_restart(restart);
            gmres.setMaxIterations(max_iter);
            gmres.compute(kernel);
            result = gmres.solveWithGuess(source,result);
            std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
            break;
                    }
    }

}

int main(int argc, char **argv) {

    unsigned int nout,max_iter_linear,restart_linear,nx;
    int fibre_resolution,fibre_number;
    double fibre_radius;
    double dt_aim,c0,h0_factor,k,gamma,rf;
    unsigned int solver_in;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("max_iter_linear", po::value<unsigned int>(&max_iter_linear)->default_value(100), "maximum iterations for linear solve")
        ("restart_linear", po::value<unsigned int>(&restart_linear)->default_value(20), "iterations until restart for linear solve")
        ("linear_solver", po::value<unsigned int>(&solver_in)->default_value(1), "linear solver")
        ("nout", po::value<unsigned int>(&nout)->default_value(10), "number of output points")
        ("k", po::value<double>(&k)->default_value(0.1), "spring constant")
        ("nx", po::value<unsigned int>(&nx)->default_value(19), "nx")
        ("fibre_resolution", po::value<int>(&fibre_resolution)->default_value(10), "number of knots around each fibre")
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
    ABORIA_VARIABLE(solution,double3,"solution")
    ABORIA_VARIABLE(solution_weights,double3,"solution weights")
    ABORIA_VARIABLE(boundary_tangent,double2,"unit tangent vector to boundary")
    ABORIA_VARIABLE(kernel_constant,double,"kernel constant")

    typedef Particles<std::tuple<boundary,inlet,outlet,interior,solution,solution_weights,kernel_constant,boundary_tangent>,2> KnotsType;
    typedef Particles<std::tuple<>,2> ParticlesType;
    typedef position_d<2> position;
    KnotsType knots;
    ParticlesType particles;
    ParticlesType fibres;

    const double mu = 1.0;
    const double flow_rate = 1.0; 
    const double Tf = 2.0;
    const double L = 1.0;
    const double delta = L/nx;
    const double boundary_layer = 0;
    const double s = 1.1*delta;
    const double h0 = h0_factor*delta;
    const int timesteps = Tf/dt_aim;
    const double dt = Tf/timesteps;
    const double dt_adapt = (1.0/100.0)*PI/sqrt(2*k);


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
        fibres.init_neighbour_search(double2(-0.1*L),double2(1.1*L),fibre_radius+boundary_layer,bool2(false));

        std::cout << "added "<<fibres.size()<<" fibres"<<std::endl;
    }
 
    // knots
    {

        typename KnotsType::value_type p;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> uniform(0,L);
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
                get<boundary_tangent>(p) = double2(-sin(kk*dtheta),cos(kk*dtheta));
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
                get<boundary_tangent>(p) = double2(1,0);
                knots.push_back(p);

                // outlet 
                get<position>(p) = double2(ii*dx,0.0);
                get<boundary>(p) = false;
                get<inlet>(p) = false;
                get<outlet>(p) = true;
                get<interior>(p) = false;
                get<boundary_tangent>(p) = double2(1,0);
                knots.push_back(p);
            }

            // boundary - left 
            get<position>(p) = double2(0.0,ii*dx);
            get<boundary>(p) = true;
            get<inlet>(p) = false;
            get<interior>(p) = false;
            get<outlet>(p) = false;
            get<boundary_tangent>(p) = double2(0,1);
            knots.push_back(p);

            // boundary - right 
            get<position>(p) = double2(L,ii*dx);
            get<boundary>(p) = true;
            get<inlet>(p) = false;
            get<outlet>(p) = false;
            get<interior>(p) = false;
            get<boundary_tangent>(p) = double2(0,1);
            knots.push_back(p);
        }

        knots.init_neighbour_search(double2(-0.1*L),double2(1.1*L),4*h0,bool2(false));
        std::cout << "added "<<knots.size()<<" knots with h0 = " <<h0<< std::endl;
    }
    
    vtkWriteGrid("init_knots",0,knots.get_grid(true));

    Symbol<boundary> is_b;
    Symbol<inlet> is_in;
    Symbol<outlet> is_out;
    Symbol<interior> is_i;
    Symbol<position> r;
    Symbol<kernel_constant> h;
    Symbol<solution> u;
    Symbol<boundary_tangent> tangent;
    Symbol<solution_weights> w;
    Label<0,KnotsType> i(knots);
    Label<1,KnotsType> j(knots);
    Label<0,ParticlesType> a(particles);
    Label<1,ParticlesType> b(particles);
    Label<0,ParticlesType> af(fibres);
    Label<1,ParticlesType> bf(fibres);
    auto dx = create_dx(i,j);
    auto dkf = create_dx(i,bf);
    auto dpf = create_dx(a,bf);
    Accumulate<std::plus<double> > sum;
    Accumulate<std::plus<double2> > sumv;
    Accumulate<std::plus<double3> > sumv3;
    Accumulate<Aboria::max<double> > max;
    max.set_init(0);
    VectorSymbolic<double,2> vector;      
    VectorSymbolic<double,3> vector3;      
    VectorSymbolic<double,9> matrix;      

    auto kernel = deep_copy(
        pow(h[i]+h[j]-norm(dx),4)*(16.0*(h[i]+h[j]) + 64.0*norm(dx))/pow(h[i]+h[j],5)
        );

    auto gradient = deep_copy(
        -320.0*dx*pow(h[i]+h[j]-norm(dx),3)/pow(h[i]+h[j],5)
        );

    auto gradientx = deep_copy(
        -320.0*dx[0]*pow(h[i]+h[j]-norm(dx),3)/pow(h[i]+h[j],5)
        );
    auto gradienty = deep_copy(
        -320.0*dx[1]*pow(h[i]+h[j]-norm(dx),3)/pow(h[i]+h[j],5)
        );

    auto laplace = deep_copy(
            pow(h[i]+h[j]-norm(dx),2)/pow(h[i]+h[j],5)*(128.0*(-2.5*(h[i]+h[j]) + 10.0*norm(dx)) - 320.0*(h[i]+h[j]-norm(dx)))
        );

    
    auto spring_force_kk = deep_copy(
        if_else(dot(dx,dx)==0
            ,0
            ,(-k*(s-norm(dx))/norm(dx))*dx
            )
        );

    auto spring_force_kf = deep_copy(
        if_else(dot(dkf,dkf)==0
            ,0
            ,(-10*k*(fibre_radius+boundary_layer-norm(dkf))/norm(dkf))*dkf
            )
        );

    auto spring_force_kb = deep_copy(
        10*k*(vector(
                if_else(r[i][0] < boundary_layer
                ,boundary_layer-r[i][0]
                ,if_else(r[i][0] > L-boundary_layer
                    ,L-boundary_layer-r[i][0]
                    ,0
                    )
                )
                ,if_else(r[i][1] < boundary_layer 
                ,boundary_layer-r[i][1]
                ,if_else(r[i][1] > L-boundary_layer
                    ,L-boundary_layer-r[i][1]
                    ,0
                    )
                ))
            )
        );
            

    // adapt knot locations
    for (int ii=0; ii<1000; ii++) {
        r[i] += dt_adapt*if_else(is_i[i]
                    ,sumv(j,norm(dx)<s,spring_force_kk)
                        + sumv(bf,norm(dkf)<fibre_radius+boundary_layer,spring_force_kf)
                        + spring_force_kb
                    ,vector(0,0)
                );
    }

    vtkWriteGrid("init_knots",1,knots.get_grid(true));

    h[i] = h0;

    //stokes - incompressible - newtonian
    //   mu * (du2/dx2 + du2/dy2) - dp/dx = 0
    //   mu * (dv2/dx2 + dv2/dy2) - dp/dy = 0
    //   du/dx + dv/dy = 0
    //
    //bc - inlet - constant (u,v) = (0,-flow_rate)
    //   - outlet - constant (u,v) = (0,-flow_rate)
    //   - boundary - constant (u,v) = (0,0)
    auto stokes = create_eigen_operator(i,j,
            if_else(is_i[i]
                ,matrix(
                     mu*laplace,0,-gradientx
                    ,0,mu*laplace,-gradienty
                    ,gradientx,gradienty,0
                    )
                ,matrix(
                    kernel,0,0
                    ,0,kernel,0
                    ,0,0,dot(gradient,tangent[i])
                )
            )
            ,norm(dx) < h[i]+h[j] 
         );

    u[i] = if_else(is_i[i] || is_b[i]
                ,vector3(0,0,0)
                ,vector3(0,-flow_rate,0));

    typedef Eigen::Matrix<double3,Eigen::Dynamic,1> vector_type; 
    typedef Eigen::Map<vector_type> map_type;
    solve(stokes,map_type(get<solution_weights>(knots).data(),knots.size()),
            map_type(get<solution>(knots).data(),knots.size()),
            max_iter_linear,restart_linear,(linear_solver)solver_in);

    // calculate solution
    u[i] = sumv3(j,norm(dx)<h[i]+h[j],kernel*w[j]);

    // write knots to file
    vtkWriteGrid("flow",0,knots.get_grid(true));
    
}

