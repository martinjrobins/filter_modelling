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
    int fibre_resolution,fibre_number;
    double fibre_radius;
    double dt_aim,c0,h0_factor,k,gamma,rf;
    unsigned int solver_in;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("max_iter_linear", po::value<unsigned int>(&max_iter_linear)->default_value(100), "maximum iterations for linear solve")
        ("restart_linear", po::value<unsigned int>(&restart_linear)->default_value(101), "iterations until restart for linear solve")
        ("linear_solver", po::value<unsigned int>(&solver_in)->default_value(2), "linear solver")
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

    ABORIA_VARIABLE(velocity_u,double,"velocity u");
    ABORIA_VARIABLE(velocity_u_dx,double,"velocity u dx");
    ABORIA_VARIABLE(divergence,double,"divergenceu dx");
    ABORIA_VARIABLE(velocity_v,double,"velocity v");
    ABORIA_VARIABLE(vorticity,double,"vorticity");
    
    ABORIA_VARIABLE(velocity_u_old,double,"old velocity u");
    ABORIA_VARIABLE(velocity_v_old,double,"old velocity v");
    ABORIA_VARIABLE(vorticity_old,double,"old vorticity");
    
    ABORIA_VARIABLE(velocity_u_weights,double,"velocity u weights");
    ABORIA_VARIABLE(velocity_v_weights,double,"velocity v weights");
    ABORIA_VARIABLE(vorticity_weights,double,"vorticity weights");

    ABORIA_VARIABLE(boundary_normal,double2,"unit normal vector to boundary")
    ABORIA_VARIABLE(kernel_constant,double,"kernel constant")

    //typedef Particles<std::tuple<boundary,inlet,outlet,interior,solution,solution_weights,kernel_constant,boundary_normal>,2> KnotsType;
    typedef Particles<std::tuple<boundary,inlet,outlet,interior,velocity_u,velocity_v,vorticity,velocity_u_old,velocity_v_old,vorticity_old,velocity_u_weights,velocity_v_weights,vorticity_weights,kernel_constant,boundary_normal,velocity_u_dx,divergence>,2> KnotsType;
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
    const double boundary_layer = delta/5;
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
                //get<boundary_tangent>(p) = double2(1,0);
                get<boundary_normal>(p) = double2(0,-1);
                knots.push_back(p);

                // outlet 
                get<position>(p) = double2(ii*dx,0.0);
                get<boundary>(p) = false;
                get<inlet>(p) = false;
                get<outlet>(p) = true;
                get<interior>(p) = false;
                //get<boundary_tangent>(p) = double2(1,0);
                get<boundary_normal>(p) = double2(0,1);
                knots.push_back(p);
            }

            // boundary - left 
            get<position>(p) = double2(0.0,ii*dx);
            get<boundary>(p) = true;
            get<inlet>(p) = false;
            get<interior>(p) = false;
            get<outlet>(p) = false;
            //get<boundary_tangent>(p) = double2(0,1);
            get<boundary_normal>(p) = double2(1,0);
            knots.push_back(p);

            // boundary - right 
            get<position>(p) = double2(L,ii*dx);
            get<boundary>(p) = true;
            get<inlet>(p) = false;
            get<interior>(p) = false;
            get<outlet>(p) = false;
            //get<boundary_tangent>(p) = double2(0,1);
            get<boundary_normal>(p) = double2(-1,0);
            knots.push_back(p);
        }

        knots.init_neighbour_search(double2(-0.1*L),double2(1.1*L),2*h0,bool2(false));
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
    //Symbol<solution_weights> w;

    Symbol<velocity_u> vu;
    Symbol<velocity_u_dx> vudx;
    Symbol<divergence> div;
    Symbol<velocity_v> vv;
    Symbol<vorticity> w;
    Symbol<velocity_u_old> vu_old;
    Symbol<velocity_v_old> vv_old;
    Symbol<vorticity_old> w_old;
    Symbol<velocity_u_weights> vuw;
    Symbol<velocity_v_weights> vvw;
    Symbol<vorticity_weights> ww;

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
        320.0*dx*pow(h[i]+h[j]-norm(dx),3)/pow(h[i]+h[j],5)
        );

    auto gradientx = deep_copy(
        320.0*dx[0]*pow(h[i]+h[j]-norm(dx),3)/pow(h[i]+h[j],5)
        );
    auto gradienty = deep_copy(
        320.0*dx[1]*pow(h[i]+h[j]-norm(dx),3)/pow(h[i]+h[j],5)
        );

    auto gradientx2 = deep_copy(

            );

    auto gradienty2 = deep_copy(

            );

    auto gradientxy = deep_copy(
            if_else(norm(dx)==0
                ,0.0
                ,960.0*dx[0]*dx[1]*pow(h[i]+h[j]-norm(dx),2)/(pow(h[i]+h[j],5)*norm(dx))
                )

            );




    auto laplace = deep_copy(
            pow(h[i]+h[j]-norm(dx),2)/pow(h[i]+h[j],5)*(128.0*(-2.5*(h[i]+h[j]) + 10.0*norm(dx)) - 320.0*(h[i]+h[j]-norm(dx)))
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
    //vorticity formulation
    //   du2/dx2 + du2/dy2 = -dwdy
    //   dv2/dx2 + dv2/dy2 = dwdx
    //   dw2/dx2 + dw2/dy2 = 0
    //
    //bc - inlet - constant u = 0 v = -flow_rate p = 0
    //   - outlet - dudy = 0 dudx = 0 dvdy = 0
    //   - boundary - constant u,v = 0 n dot nabla p = 0 
    /*
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
                        ,0,0,dot(gradient,normal[i])
                )
            )
            ,norm(dx) < h[i]+h[j] 
         );
    */

    /*

    auto UU = create_eigen_operator(i,j,
            if_else(is_i[i] || is_out[i]
                ,mu*laplace
                ,kernel
                )
            ,norm(dx) < h[i]+h[j] 
         );

    auto UP = create_eigen_operator(i,j,
            if_else(is_i[i] || is_out[i]
                ,-gradientx
                ,0.0
                )
            ,norm(dx) < h[i]+h[j] 
         );

    auto VV = create_eigen_operator(i,j,
            if_else(is_i[i] || is_out[i]
                ,mu*laplace
                ,kernel
                )
            ,norm(dx) < h[i]+h[j] 
         );

    auto VP = create_eigen_operator(i,j,
            if_else(is_i[i] || is_out[i]
                ,-gradienty
                ,0.0
                )
            ,norm(dx) < h[i]+h[j] 
         );

    auto PU = create_eigen_operator(i,j,
            if_else(is_i[i] || is_out[i]
                ,gradientx
                ,0.0
                )
            ,norm(dx) < h[i]+h[j] 
         );

    auto PV = create_eigen_operator(i,j,
            if_else(is_i[i] || is_out[i]
                ,gradienty
                ,0.0
                )
            ,norm(dx) < h[i]+h[j] 
         );

    auto PP = create_eigen_operator(i,j,
            if_else(is_i[i] || is_out[i]
                ,0.0
                ,kernel
                )
            ,norm(dx) < h[i]+h[j] 
         );

    auto Zero = create_eigen_operator(i,j, 0.);

    auto stokes = create_block_eigen_operator<3,3>( UU,Zero,UP
                                                   ,Zero,VV,VP
                                                   ,PU,PV,PP);

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type; 
    matrix_type A(3*knots.size(),3*knots.size());
    for (int ii=0; ii<3*knots.size(); ++ii) {
        int num_non_zero_row = 0;
        int num_non_zero_col = 0;
        for (int jj=0; jj<3*knots.size(); ++jj) {
            if (stokes.coeff(ii,jj) != 0) {
                num_non_zero_row++;
            }
            if (stokes.coeff(jj,ii) != 0) {
                num_non_zero_col++;
            }

            if (ii != jj) {
                assert((get<position>(knots)[ii] - get<position>(knots)[jj]).norm()>0.0001);
            }
            A(ii,jj) = stokes.coeff(ii,jj);
        }
        assert(num_non_zero_row > 0);
        assert(num_non_zero_col > 0);
    }


    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
    typedef Eigen::Map<vector_type> map_type;
    vector_type w_eigen(3*knots.size());
    w_eigen.setZero(3*knots.size());
    vector_type u_eigen(3*knots.size());
    u_eigen.setZero(3*knots.size());

    for (int ii=0; ii<knots.size(); ++ii) {
        if (get<interior>(knots)[ii] || get<boundary>(knots)[ii] || get<outlet>(knots)[ii]) {
            u_eigen[ii] = 0;
            u_eigen[knots.size()+ii] = 0;
            u_eigen[2*knots.size()+ii] = 0;
        } else {
            u_eigen[ii] = 0;
            u_eigen[knots.size()+ii] = -1.0;
            u_eigen[2*knots.size()+ii] = 0.0;
        }
    }

    w_eigen = A.colPivHouseholderQr().solve(u_eigen);
    //w_eigen = A.llt().solve(u_eigen);
    //
    solve(stokes,w_eigen,
        u_eigen,
        max_iter_linear,restart_linear,(linear_solver)solver_in);
    double relative_error = (A*w_eigen - u_eigen).norm() / u_eigen.norm();               
    std::cout << "The relative error is:\n" << relative_error << std::endl;
        */

    ww[i] = 0.0;
    w[i] = 0.0;
    vu[i] = 0.0;
    vv[i] = 0.0;
    
    double tol = 1;
    int ii = -1;
    while (tol > 1e-4) {
        ii++;
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        typedef Eigen::Map<vector_type> map_type;

        //save old values
        vu_old[i] = vu[i];
        vv_old[i] = vv[i];
        w_old[i] = w[i];

        // adapt knot locations
        //  potential = gradientx u * gradienty v - gradienty u*gradientx v
        //  potential = sqrt(gradientx u ^2  +  gradient y u ^2) + 
        //              sqrt(gradientx v ^2  +  gradient y v ^2)
        //  f = gradient potential
        //    = (gx2 u * gy v + gx u * gxy v - gxy u * gx v + gy u*gx2 v,
        //       gxy u * gy v + gx y * gy2 v - gy2 u * gx v + gy u*gxy v)
        for (int ii=0; ii<100; ii++) {
            r[i] += dt_adapt*if_else(is_i[i]
                        ,sumv(j,norm(dx)<s,spring_force_kk)
                            + sumv(bf,norm(dkf)<fibre_radius+boundary_layer,spring_force_kf)
                            + spring_force_kb
                            + 
                        ,vector(0,0)
                    );
        }

        //vorticity formulation
        //   du2/dx2 + du2/dy2 = -dwdy
        //   dv2/dx2 + dv2/dy2 = dwdx
        //   dw2/dx2 + dw2/dy2 = 0

        //solve for u and v
        auto A = create_eigen_operator(i,j,
                if_else(is_i[i]
                   ,laplace
                   //,if_else(is_out[i]
                   //    ,gradienty
                       ,kernel
                   //)
                )
                ,norm(dx) < h[i]+h[j] 
            );


        vu[i] = if_else(is_i[i]
                    ,-sum(j,norm(dx)<h[i]+h[j],gradienty*ww[j])
                    ,0.0
                );
        vv[i] = if_else(is_i[i]
                    ,sum(j,norm(dx)<h[i]+h[j],gradientx*ww[j])
                    ,if_else(is_in[i] || is_out[i]
                        //,-flow_rate*(1-pow(r[i][0]-0.5*L,2)/std::pow(0.5*L,2))
                        ,-flow_rate
                        ,0.0
                     )
                );

        

        solve(A,map_type(get<velocity_u_weights>(knots).data(),knots.size()),
                map_type(get<velocity_u>(knots).data(),knots.size()),
                max_iter_linear,restart_linear,(linear_solver)solver_in);

        solve(A,map_type(get<velocity_v_weights>(knots).data(),knots.size()),
                map_type(get<velocity_v>(knots).data(),knots.size()),
                max_iter_linear,restart_linear,(linear_solver)solver_in);

        vu[i] = sum(j,norm(dx)<h[i]+h[j],kernel*vuw[j]);
        vv[i] = sum(j,norm(dx)<h[i]+h[j],kernel*vvw[j]);
        vudx[i] = sum(j,norm(dx)<h[i]+h[j],gradientx*vuw[j]);
        div[i] = sum(j,norm(dx)<h[i]+h[j],gradientx*vuw[j] + gradienty*vvw[j]);


        //solve for w
        auto B = create_eigen_operator(i,j,
                if_else(is_i[i]
                   ,laplace
                   ,kernel
                )
                ,norm(dx) < h[i]+h[j] 
            );

        //calculate vorticity boundary conditions
        w[i] = if_else(is_i[i]
                    ,0.0 
                    ,sum(j,norm(dx)<h[i]+h[j],gradientx*vvw[j]) 
                        - sum(j,norm(dx)<h[i]+h[j],gradienty*vuw[j]) 
                );

        //solve for w
        solve(B,map_type(get<vorticity_weights>(knots).data(),knots.size()),
                map_type(get<vorticity>(knots).data(),knots.size()),
                max_iter_linear,restart_linear,(linear_solver)solver_in);

        //calculate solution
        w[i] = sum(j,norm(dx)<h[i]+h[j],kernel*ww[j]);
        
        //calculate tol
        tol = eval(sum(i,true,pow(vu[i]-vu_old[i],2)+pow(vv[i]-vv_old[i],2)+pow(w[i]-w_old[i],2)));

        vtkWriteGrid("iteration",ii,knots.get_grid(true));
        std::cout << "finished iteration "<<ii<<", tol = "<<tol<<std::endl;
    }


    /*
    for (int ii=0; ii<knots.size(); ++ii) {
        get<solution_weights>(knots)[ii][0] = w_eigen[ii];
        get<solution_weights>(knots)[ii][1] = w_eigen[knots.size()+ii];
        get<solution_weights>(knots)[ii][2] = w_eigen[2*knots.size()+ii];

        //get<solution>(knots)[ii][0] = u_eigen[ii];
        //get<solution>(knots)[ii][1] = u_eigen[knots.size()+ii];
        //get<solution>(knots)[ii][2] = u_eigen[2*knots.size()+ii];
    }
    
    // calculate solution
    u[i] = sumv3(j,norm(dx)<h[i]+h[j],kernel*w[j]);

    // write knots to file
    vtkWriteGrid("flow",0,knots.get_grid(true));
    */
    
}

