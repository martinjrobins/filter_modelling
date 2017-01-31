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
ABORIA_VARIABLE(alpha,double3,"alpha");
ABORIA_VARIABLE(kernel_constant,double,"kernel constant")

typedef Particles<std::tuple<alpha,boundary,inlet,outlet,interior,velocity_u,velocity_v,pressure,kernel_constant>,2> KnotsType;
typedef Particles<std::tuple<dvelocity_u,dvelocity_v,dpressure,velocity_u,velocity_v,pressure,kernel_constant>,2> ComsolType;
typedef Particles<std::tuple<kernel_constant>,2> ParticlesType;
typedef position_d<2> position;


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





template <typename I, typename J, typename C>
auto gen_kernel_mq(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    return deep_copy(
        sqrt(dot(dx,dx)+pow(c[j],2))
        );
}

template <typename I, typename J, typename C>
auto gen_phi_sol(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    const double mu = 1.0;

    return deep_copy(
        ((1.0/75.0)*sqrt(dot(dx,dx)+pow(c[j],2))*(4.0*pow(dot(dx,dx),2)+48.0*dot(dx,dx)*pow(c[j],2)-61.0*pow(c[j],4)) - pow(c[j],3)*log(c[j])*dot(dx,dx) - (1.0/5.0)*(5.0*dot(dx,dx)-2*pow(c[j],2))*pow(c[j],3)*log(c[j]+sqrt(dot(dx,dx)+pow(c[j],2))))/(12.0*mu)
        );
}

template <typename I, typename J, typename C>
auto gen_phi_sol_dash_div_r(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    return deep_copy(
        (4.0*pow(dot(dx,dx),3) + 36.0*pow(dot(dx,dx),2)*pow(c[j],2) + 39.0*dot(dx,dx)*pow(c[j],4) + 7.0*pow(c[j],6))/(180.0*pow(dot(dx,dx)+pow(c[j],2),1.5))
            - (5.0*pow(dot(dx,dx),2) + 3.0*dot(dx,dx)*pow(c[j],2) - 2.0*pow(c[j],4))*pow(c[j],3)/(60.0*pow(dot(dx,dx)+pow(c[j],2),1.5)*(c[j]+sqrt(dot(dx,dx)+pow(c[j],2))))
            - (1.0/6.0)*pow(c[j],3)*log(c[j]+sqrt(dot(dx,dx)+pow(c[j],2))) - (1.0/6.0)*pow(c[j],3)*log(c[j])
            );
}


template <typename I, typename J, typename C>
auto gen_phi_sol_dash_dash(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    return deep_copy(
        (16.0*pow(dot(dx,dx),3)+84.0*pow(dot(dx,dx),2)*pow(c[j],2)+96.0*dot(dx,dx)*pow(c[j],4)+7*pow(c[j],6))/(180.0*pow(dot(dx,dx)+pow(c[j],2),1.5))
            - (20.0*pow(dot(dx,dx),2)+25.0*dot(dx,dx)*pow(c[j],2)-2.0*pow(c[j],4))*pow(c[j],3)/(60.0*pow(dot(dx,dx)+pow(c[j],2),1.5)*(c[j]+sqrt(dot(dx,dx)+pow(c[j],2))))
            + (5.0*dot(dx,dx)-2*pow(c[j],2))*pow(c[j],3)*dot(dx,dx)/(60.0*(dot(dx,dx)+pow(c[j],2))*(c[j]+sqrt(dot(dx,dx)+pow(c[j],2))))
            - (1.0/6.0)*pow(c[j],3)*log(c[j]+sqrt(dot(dx,dx)+pow(c[j],2))) - (1.0/6.0)*pow(c[j],3)*log(c[j])
            );
}

 
template <typename I, typename J, typename C>
auto gen_phi_sol_dash_dash_dash(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    return deep_copy(
        (76.0*pow(dot(dx,dx),2)+176.0*dot(dx,dx)*pow(c[j],2)+285.0*pow(c[j],4))*norm(dx)/(300.0*pow(dot(dx,dx)+pow(c[j],2),1.5))
            + (4.0*pow(dot(dx,dx),2)+48.0*dot(dx,dx)*pow(c[j],2)-61.0*pow(c[j],4))*pow(norm(dx),3)/(300.0*pow(dot(dx,dx)+pow(c[j],2),5.0/2.0))
            + (10.0*pow(dot(dx,dx),2)+15.0*dot(dx,dx)*pow(c[j],2)-2.0*pow(c[j],4))*pow(c[j],3)*norm(dx)/(20.0*pow(dot(dx,dx)+pow(c[j],2),2)*(c[j]+sqrt(dot(dx,dx)+pow(c[j],2))))
            - (5.0*dot(dx,dx)+22.0*pow(c[j],2))*pow(c[j],3)*norm(dx)/(20.0*pow(dot(dx,dx)+pow(c[j],2),1.5)*(c[j]+sqrt(dot(dx,dx)+pow(c[j],2))))
            + (-5.0*dot(dx,dx)+2*pow(c[j],2))*pow(c[j],3)*pow(norm(dx),3)/(20.0*pow(dot(dx,dx)+pow(c[j],2),5.0/2.0)*(c[j]+sqrt(dot(dx,dx)+pow(c[j],2))))
            + (-5.0*dot(dx,dx)+2*pow(c[j],2))*pow(c[j],3)*pow(norm(dx),3)/(30.0*pow(dot(dx,dx)+pow(c[j],2),1.5)*pow(c[j]+sqrt(dot(dx,dx)+pow(c[j],2)),3))
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



void read_data_files(ComsolType &particles);
void setup_knots(KnotsType &knots, ParticlesType &fibres, const double fibre_radius, const double fibre_resolution_factor, const double nx, double2 domain_min, double2 domain_max, const double c0, const double k, const double epsilon_strength, const double epsilon_falloff);
void calculate_c(KnotsType &knots, double c0, const double nx, double2 domain_min, double2 domain_max);
void solve_stokes_MAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in);
void solve_stokes_LMAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double delta);
    
