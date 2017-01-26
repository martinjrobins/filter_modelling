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
                         Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
            cg.setMaxIterations(max_iter);
            cg.compute(kernel);
            result = cg.solveWithGuess(source,result);
            std::cout << "CG:    #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << " actual error: "<<(kernel*result - source).norm() / source.norm() << std::endl;
            break;
                 }
        case BiCGSTAB: {
            Eigen::BiCGSTAB<
                typename std::remove_reference<Kernel>::type, 
                     Eigen::IdentityPreconditioner> bicg;
            bicg.setMaxIterations(max_iter);
            bicg.compute(kernel);
            result = bicg.solveWithGuess(source,result);
            //result = bicg.solve(source);
            std::cout << "BiCGSTAB:    #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << " actual error: "<<(kernel*result - source).norm() / source.norm() << std::endl;
            break;
               }
        case GMRES: {
                        /*
            Eigen::GMRES<
                typename std::remove_reference<Kernel>::type, 
                    Eigen::DiagonalPreconditioner<Scalar>> gmres;
                    */
            Eigen::GMRES<
                typename std::remove_reference<Kernel>::type, 
                    Eigen::IdentityPreconditioner> gmres;

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
                    Eigen::IdentityPreconditioner> dgmres;
            dgmres.set_restart(restart);
            dgmres.setMaxIterations(max_iter);
            dgmres.compute(kernel);
            result = dgmres.solveWithGuess(source,result);
            std::cout << "DGMRES:    #iterations: " << dgmres.iterations() << ", estimated error: " << dgmres.error() << " actual error: "<<(kernel*result - source).norm() / source.norm()<<std::endl;
            break;
                    }


    }


}

ABORIA_VARIABLE(inlet,uint8_t,"is inlet knot")
ABORIA_VARIABLE(outlet,uint8_t,"is outlet knot")
ABORIA_VARIABLE(boundary,uint8_t,"is boundary knot")
ABORIA_VARIABLE(interior,uint8_t,"is interior knot")
ABORIA_VARIABLE(velocity_u,double,"velocity u");
ABORIA_VARIABLE(velocity_v,double,"velocity v");
ABORIA_VARIABLE(pressure,double,"pressure");
ABORIA_VARIABLE(dvelocity_u,double,"error velocity u");
ABORIA_VARIABLE(dvelocity_v,double,"error velocity v");
ABORIA_VARIABLE(dpressure,double,"error comsol pressure");
ABORIA_VARIABLE(alpha,double2,"alpha");
ABORIA_VARIABLE(kernel_constant,double,"kernel constant")



template <typename Particles> 
void read_data_files(Particles &particles) {

    typedef typename Particles::position position;

    std::cout << "reading files..." << std::endl;
    std::ifstream pressure_file("five layer_square/[p]_five_layer_square.txt" );
    std::ifstream vel_horz_file("five layer_square/[vel_horz]_five_layer_square.txt" );
    std::ifstream vel_vert_file("five layer_square/[vel_vert]_five_layer_square.txt" );
    std::string line;
    while ( pressure_file.good() ) {
        double2 pos,velocity;
        double pressure,dummy;
        std::getline(pressure_file, line);
        std::istringstream buffer(line);
        buffer >> pos[0];
        buffer >> pos[1];
        buffer >> pressure;

        std::getline(vel_horz_file, line);
        buffer >> dummy;
        buffer >> dummy;
        buffer >> velocity[0];

        std::getline(vel_vert_file, line);
        buffer >> dummy;
        buffer >> dummy;
        buffer >> velocity[1];

        typename Particles::value_type p;
        get<position>(p) = pos;
        get<dvelocity_u>(p) = velocity[0];
        get<dvelocity_v>(p) = velocity[1];
        get<dpressure>(p) = pressure;
        particles.push_back(p);
    }
    std::cout << "done reading files"<< std::endl;
}

template <typename I, typename J>
auto gen_spring(I &i, J &j, const double &k, const double &s) {
    auto dx = create_dx(i,j);
    VectorSymbolic<double,2> vector;      
    return deep_copy(
            if_else(dot(dx,dx)==0
            ,vector(0,0)
            ,(-k*(s-norm(dx))/norm(dx))*dx
            )
            );
}

template<typename KnotsType, typename FibresType>
void setup_knots(KnotsType &knots, FibresType &fibres, const double fibre_radius, const double fibre_resolution_factor, const double nx, double2 domain_min, double2 domain_max, const double c0, const double k, const double epsilon_strength, const double epsilon_falloff) {

    std::cout << "setup knots..." << std::endl;

    typedef typename KnotsType::position position;
    const double L = domain_max[0] - domain_min[0];
    const double delta = L/nx;
    const double fibre_resolution = delta*fibre_resolution_factor;
    const double s = 1.1*delta;

    const double boundary_layer = delta/5;
    const double volume = ((domain_max-domain_min).prod() 
                            - fibres.size()*PI*std::pow(fibre_radius,2))
                            /(domain_max-domain_min).prod();
    const int N = nx*nx*(domain_max[1]-domain_min[1])*volume/(domain_max[0]-domain_min[0]); 
    const double2 ns_buffer(L/2,1e-10);

    const double dt_adapt = (1.0/100.0)*PI/sqrt(2*k);

    typename KnotsType::value_type p;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniformx(domain_min[0]+delta/10.0,domain_max[0]-delta/10.0);
    std::uniform_real_distribution<double> uniformy(domain_min[1]+delta/10.0,domain_max[1]-delta/10.0);
    while (knots.size() < N) {
        get<position>(p) = double2(uniformx(generator),uniformy(generator));
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
            get<kernel_constant>(p) = c0;
            knots.push_back(p);
        }

    }
    std::cout << "added "<<knots.size()<<" interior knots"<<std::endl;

    // fibre boundaries
    for (int ii=0; ii<fibres.size(); ++ii) {
        const double2 origin = get<position>(fibres)[ii];
        const double dtheta_aim = fibre_resolution/fibre_radius;
        const int n = std::ceil(2*PI/dtheta_aim);
        const double dtheta = 2*PI/n;
        for (int kk=0; kk<n; ++kk) {
            get<position>(p) = origin + fibre_radius*double2(cos(kk*dtheta),sin(kk*dtheta));
            get<boundary>(p) = true;
            get<inlet>(p) = false;
            get<outlet>(p) = false;
            get<interior>(p) = false;
            get<kernel_constant>(p) = c0;
            knots.push_back(p);
        }
    }

    for (int ii=1; ii<nx; ++ii) {
        const double dx = L/nx;
        // inlet
        get<position>(p) = double2(domain_min[0] + ii*dx,domain_max[1]);
        get<boundary>(p) = false;
        get<inlet>(p) = true;
        get<outlet>(p) = false;
        get<interior>(p) = false;
        get<kernel_constant>(p) = c0;
        knots.push_back(p);

        // outlet 
        get<position>(p) = double2(domain_min[0] + ii*dx,domain_min[1]);
        get<boundary>(p) = false;
        get<inlet>(p) = false;
        get<outlet>(p) = true;
        get<interior>(p) = false;
        get<kernel_constant>(p) = c0;
        knots.push_back(p);
    }
    const int ny = (domain_max[1]-domain_min[1])/delta;
    const double deltay = (domain_max[1]-domain_min[1])/ny;
    for (int ii=0; ii<=ny; ++ii) {
        // boundary - left 
        get<position>(p) = double2(domain_min[0],domain_min[1]+ii*deltay);
        get<boundary>(p) = true;
        get<inlet>(p) = false;
        get<interior>(p) = false;
        get<outlet>(p) = false;
        get<kernel_constant>(p) = c0;
        knots.push_back(p);

        // boundary - right 
        get<position>(p) = double2(domain_max[0],domain_min[1]+ii*deltay);
        get<boundary>(p) = true;
        get<inlet>(p) = false;
        get<interior>(p) = false;
        get<outlet>(p) = false;
        get<kernel_constant>(p) = c0;
        knots.push_back(p);
    }

    knots.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,L/5,bool2(false));
    std::cout << "added "<<knots.size()<<" knots with c0 = " <<c0<< std::endl;
    
    vtkWriteGrid("init_knots",0,knots.get_grid(true));

    Symbol<boundary> is_b;
    Symbol<inlet> is_in;
    Symbol<outlet> is_out;
    Symbol<interior> is_i;
    Symbol<position> r;
    Symbol<alive> alive_;
    Symbol<kernel_constant> c;

    Label<0,KnotsType> i(knots);
    Label<1,KnotsType> j(knots);
    Label<1,FibresType> bf(fibres);
    auto dx = create_dx(i,j);
    auto dkf = create_dx(i,bf);
    Accumulate<std::plus<double2> > sumv;
    VectorSymbolic<double,2> vector;      
    
    auto spring_force_kk = gen_spring(i,j,k,s);

    auto spring_force_kf = deep_copy(
        if_else(dot(dkf,dkf)==0
            ,vector(0,0)
            ,if_else(norm(dkf)>fibre_radius+boundary_layer
                //,0.0
                ,epsilon_strength*delta*exp(1.0/epsilon_falloff*(fibre_radius+boundary_layer-norm(dkf)))
                ,-10*k*(fibre_radius+boundary_layer-norm(dkf))
                )*dkf/norm(dkf)
            )
        );


    auto spring_force_kb = deep_copy(
        vector(
                if_else(r[i][0] < domain_min[0]+boundary_layer
                ,10*k*(domain_min[0]+boundary_layer-r[i][0])
                ,if_else(r[i][0] > domain_max[0]-boundary_layer
                    ,10*k*(domain_max[0]-boundary_layer-r[i][0])
                    //,0.0
                    ,-epsilon_strength*delta*(exp(1.0/epsilon_falloff*(domain_min[0]+boundary_layer-r[i][0]))
                                          - exp(1.0/delta*(r[i][0]-domain_max[0]+boundary_layer)))
                    )
                )
                ,if_else(r[i][1] < domain_min[1]+boundary_layer 
                ,10*k*(domain_min[1]+boundary_layer-r[i][1])
                ,if_else(r[i][1] > domain_max[1]-boundary_layer
                    ,10*k*(domain_max[1]-boundary_layer-r[i][1])
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

    std::cout << "finshed adapting!"<<std::endl;

    vtkWriteGrid("init_knots",1,knots.get_grid(true));
}





template <typename I, typename J, typename C>
auto gen_kernel_mq(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    return deep_copy(
        sqrt(dot(dx,dx)+pow(c[i],2))
        );
}

template <typename I, typename J, typename C>
auto gen_phi_sol(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    const double mu = 1.0;

    return deep_copy(
        ((1.0/75.0)*sqrt(dot(dx,dx)+pow(c[i],2))*(4.0*pow(dot(dx,dx),2)+48.0*dot(dx,dx)*pow(c[i],2)-61.0*pow(c[i],4)) - pow(c[i],3)*log(c[i])*dot(dx,dx) - (1.0/5.0)*(5.0*dot(dx,dx)-2*pow(c[i],2))*pow(c[i],3)*log(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))/(12.0*mu)
        );
}

template <typename I, typename J, typename C>
auto gen_phi_sol_dash_div_r(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    return deep_copy(
        (4.0*pow(dot(dx,dx),3) + 36.0*pow(dot(dx,dx),2)*pow(c[i],2) + 39.0*dot(dx,dx)*pow(c[i],4) + 7.0*pow(c[i],6))/(180.0*pow(dot(dx,dx)+pow(c[i],2),1.5))
            - (5.0*pow(dot(dx,dx),2) + 3.0*dot(dx,dx)*pow(c[i],2) - 2.0*pow(c[i],4))*pow(c[i],3)/(60.0*pow(dot(dx,dx)+pow(c[i],2),1.5)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            - (1.0/6.0)*pow(c[i],3)*log(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))) - (1.0/6.0)*pow(c[i],3)*log(c[i])
            );
}


template <typename I, typename J, typename C>
auto gen_phi_sol_dash_dash(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    return deep_copy(
        (16.0*pow(dot(dx,dx),3)+84.0*pow(dot(dx,dx),2)*pow(c[i],2)+96.0*dot(dx,dx)*pow(c[i],4)+7*pow(c[i],6))/(180.0*pow(dot(dx,dx)+pow(c[i],2),1.5))
            - (20.0*pow(dot(dx,dx),2)+25.0*dot(dx,dx)*pow(c[i],2)-2.0*pow(c[i],4))*pow(c[i],3)/(60.0*pow(dot(dx,dx)+pow(c[i],2),1.5)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            + (5.0*dot(dx,dx)-2*pow(c[i],2))*pow(c[i],3)*dot(dx,dx)/(60.0*(dot(dx,dx)+pow(c[i],2))*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            - (1.0/6.0)*pow(c[i],3)*log(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))) - (1.0/6.0)*pow(c[i],3)*log(c[i])
            );
}

 
template <typename I, typename J, typename C>
auto gen_phi_sol_dash_dash_dash(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    return deep_copy(
        (76.0*pow(dot(dx,dx),2)+176.0*dot(dx,dx)*pow(c[i],2)+285.0*pow(c[i],4))*norm(dx)/(300.0*pow(dot(dx,dx)+pow(c[i],2),1.5))
            + (4.0*pow(dot(dx,dx),2)+48.0*dot(dx,dx)*pow(c[i],2)-61.0*pow(c[i],4))*pow(norm(dx),3)/(300.0*pow(dot(dx,dx)+pow(c[i],2),5.0/2.0))
            + (10.0*pow(dot(dx,dx),2)+15.0*dot(dx,dx)*pow(c[i],2)-2.0*pow(c[i],4))*pow(c[i],3)*norm(dx)/(20.0*pow(dot(dx,dx)+pow(c[i],2),2)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            - (5.0*dot(dx,dx)+22.0*pow(c[i],2))*pow(c[i],3)*norm(dx)/(20.0*pow(dot(dx,dx)+pow(c[i],2),1.5)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            + (-5.0*dot(dx,dx)+2*pow(c[i],2))*pow(c[i],3)*pow(norm(dx),3)/(20.0*pow(dot(dx,dx)+pow(c[i],2),5.0/2.0)*(c[i]+sqrt(dot(dx,dx)+pow(c[i],2))))
            + (-5.0*dot(dx,dx)+2*pow(c[i],2))*pow(c[i],3)*pow(norm(dx),3)/(30.0*pow(dot(dx,dx)+pow(c[i],2),1.5)*pow(c[i]+sqrt(dot(dx,dx)+pow(c[i],2)),3))
            );
}


template <typename I, typename J, typename C>
auto gen_psol_u1(I &i, J &j, C &c) { 
    auto dx = create_dx(i,j);
    auto phi_sol_dash_div_r = gen_phi_sol_dash_div_r(i,j,c);
    auto phi_sol_dash_dash = gen_phi_sol_dash_dash(i,j,c);
    auto phi_sol_dash_dash_dash = gen_phi_sol_dash_dash_dash(i,j,c);

    const double mu = 1.0;
    // i = 1, l = 1
    return deep_copy(
        if_else(norm(dx)==0
            ,(1.0/mu)*(phi_sol_dash_div_r + phi_sol_dash_dash)
            ,(1.0/mu)*(phi_sol_dash_div_r*dx[0]*dx[0] + phi_sol_dash_dash*dx[1]*dx[1])/dot(dx,dx)
            )
        );
}

template <typename I, typename J, typename C>
auto gen_psol_u2(I &i, J &j, C &c) { 
    auto dx = create_dx(i,j);
    auto phi_sol_dash_div_r = gen_phi_sol_dash_div_r(i,j,c);
    auto phi_sol_dash_dash = gen_phi_sol_dash_dash(i,j,c);
    auto phi_sol_dash_dash_dash = gen_phi_sol_dash_dash_dash(i,j,c);

    const double mu = 1.0;
    // i = 1, l = 2
    return deep_copy(
        if_else(norm(dx)==0
            ,(1.0/mu)*(phi_sol_dash_div_r - phi_sol_dash_dash)
            ,(1.0/mu)*(phi_sol_dash_div_r - phi_sol_dash_dash)*(dx[0]*dx[1])/dot(dx,dx)
            )
        );
}

template <typename I, typename J, typename C>
auto gen_psol_v1(I &i, J &j, C &c) { 
    auto dx = create_dx(i,j);
    auto phi_sol_dash_div_r = gen_phi_sol_dash_div_r(i,j,c);
    auto phi_sol_dash_dash = gen_phi_sol_dash_dash(i,j,c);
    auto phi_sol_dash_dash_dash = gen_phi_sol_dash_dash_dash(i,j,c);

    const double mu = 1.0;
    // i = 2, l = 1
    return deep_copy(
        if_else(norm(dx)==0
            ,(1.0/mu)*(phi_sol_dash_div_r - phi_sol_dash_dash)
            ,(1.0/mu)*(phi_sol_dash_div_r - phi_sol_dash_dash)*(dx[0]*dx[1])/dot(dx,dx)
            )
        );
}

template <typename I, typename J, typename C>
auto gen_psol_v2(I &i, J &j, C &c) { 
    auto dx = create_dx(i,j);
    auto phi_sol_dash_div_r = gen_phi_sol_dash_div_r(i,j,c);
    auto phi_sol_dash_dash = gen_phi_sol_dash_dash(i,j,c);
    auto phi_sol_dash_dash_dash = gen_phi_sol_dash_dash_dash(i,j,c);

    const double mu = 1.0;
    return deep_copy(
        if_else(norm(dx)==0
            ,(1.0/mu)*(phi_sol_dash_div_r+phi_sol_dash_dash)
            ,(1.0/mu)*(phi_sol_dash_div_r*dx[1]*dx[1] + phi_sol_dash_dash*dx[0]*dx[0])/dot(dx,dx)
            )
        );
}

template <typename I, typename J, typename C>
auto gen_psol_p1(I &i, J &j, C &c) { 
    auto dx = create_dx(i,j);
    auto phi_sol_dash_div_r = gen_phi_sol_dash_div_r(i,j,c);
    auto phi_sol_dash_dash = gen_phi_sol_dash_dash(i,j,c);
    auto phi_sol_dash_dash_dash = gen_phi_sol_dash_dash_dash(i,j,c);

    return deep_copy(
        if_else(norm(dx)==0
            ,-(phi_sol_dash_dash_dash)
            ,(phi_sol_dash_dash_dash + (phi_sol_dash_dash - phi_sol_dash_div_r)/norm(dx))*dx[0]/norm(dx)
            )
        );
}

template <typename I, typename J, typename C>
auto gen_psol_p2(I &i, J &j, C &c) { 
    auto dx = create_dx(i,j);
    auto phi_sol_dash_div_r = gen_phi_sol_dash_div_r(i,j,c);
    auto phi_sol_dash_dash = gen_phi_sol_dash_dash(i,j,c);
    auto phi_sol_dash_dash_dash = gen_phi_sol_dash_dash_dash(i,j,c);

    return deep_copy(
        if_else(norm(dx)==0
            ,-(phi_sol_dash_dash_dash)
            ,(phi_sol_dash_dash_dash + (phi_sol_dash_dash - phi_sol_dash_div_r)/norm(dx))*dx[1]/norm(dx)
            )
        );
}



template <typename KnotsType>
void solve_stokes(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in) {
    std::cout << "solving stokes..."<<knots.size()<<std::endl;

    const double flow_rate = 1.0; 

    typedef typename KnotsType::position position;

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

    Label<0,KnotsType> i(knots);
    Label<1,KnotsType> j(knots);
    auto dx = create_dx(i,j);
    Accumulate<std::plus<double> > sum;
    Accumulate<std::plus<double2> > sumv;
    Accumulate<std::plus<double3> > sumv3;
    VectorSymbolic<double,2> vector;      

    auto kernel_mq = gen_kernel_mq(i,j,c);
    auto psol_p1 = gen_psol_p1(i,j,c);
    auto psol_p2 = gen_psol_p2(i,j,c);
    auto psol_u1 = gen_psol_u1(i,j,c);
    auto psol_u2 = gen_psol_u2(i,j,c);
    auto psol_v1 = gen_psol_v1(i,j,c);
    auto psol_v2 = gen_psol_v2(i,j,c);

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


    std::cout << "setup equations..."<<std::endl;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type; 
    typedef Eigen::Map<vector_type> map_type;
    map_type eigen_pressure(get<pressure>(knots).data(),knots.size());
    map_type eigen_velocity_u(get<velocity_u>(knots).data(),knots.size());
    map_type eigen_velocity_v(get<velocity_u>(knots).data(),knots.size());
    

    vector_type source(2*knots.size());
    vector_type alphas(2*knots.size());
    vector_type alphas_test(2*knots.size());
    matrix_type A_eigen(2*knots.size(),2*knots.size());
    std::cout << "assemble matrix..."<<std::endl;
    A.assemble(A_eigen);

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

    std::cout << "solve ..."<<std::endl;
    solve(A_eigen,alphas,source,
                max_iter_linear,restart_linear,(linear_solver)solver_in);


    std::cout << "done solve..."<<std::endl;
    for (int ii=0; ii<knots.size(); ++ii) {
        get<alpha>(knots)[ii][0] = alphas[ii];
        get<alpha>(knots)[ii][1] = alphas[ii+knots.size()];
    }

    vu[i] = sum(j,true,psol_u1*al[j][0] + psol_u2*al[j][1]);
    vv[i] = sum(j,true,psol_v1*al[j][0] + psol_v2*al[j][1]);
    pr[i] = sum(j,true,psol_p1*al[j][0] + psol_p2*al[j][1]);

    vtkWriteGrid("MAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
}


