#pragma once

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

typedef Eigen::Vector2d eigen_vector;
ABORIA_VARIABLE(inlet,uint8_t,"is_inlet_knot")
ABORIA_VARIABLE(outlet,uint8_t,"is_outlet_knot")
ABORIA_VARIABLE(boundary,uint8_t,"is_boundary_knot")
ABORIA_VARIABLE(target,uint8_t,"is_target_knot")
ABORIA_VARIABLE(interior,uint8_t,"is_interior_knot")
ABORIA_VARIABLE(velocity,eigen_vector,"velocity");
ABORIA_VARIABLE(stokes_velocity,eigen_vector,"stokes_velocity");
ABORIA_VARIABLE(electro_velocity,eigen_vector,"electro_velocity");
ABORIA_VARIABLE(velocity_dudx,double,"velocity_dudx");
ABORIA_VARIABLE(velocity_dudy,double,"velocity_dudy");
ABORIA_VARIABLE(velocity_dvdx,double,"velocity_dvdx");
ABORIA_VARIABLE(velocity_dvdy,double,"velocity_dvdy");
ABORIA_VARIABLE(pressure,double,"pressure");
ABORIA_VARIABLE(dvelocity,eigen_vector,"error_velocity_u");
ABORIA_VARIABLE(dpressure,double,"error_comsol_pressure");
ABORIA_VARIABLE(alpha1,double,"alpha_1");
ABORIA_VARIABLE(alpha2,double,"alpha_2");
ABORIA_VARIABLE(count,int,"count");
ABORIA_VARIABLE(angle,double,"angle");
ABORIA_VARIABLE(charge,double,"charge");
ABORIA_VARIABLE(kernel_constant,double,"kernel_constant")
ABORIA_VARIABLE(point_a,vdouble2,"start point of element")
ABORIA_VARIABLE(point_b,vdouble2,"end point of element")
ABORIA_VARIABLE(traction,eigen_vector,"traction of element")
ABORIA_VARIABLE(knotpot,double,"knot potential")
ABORIA_VARIABLE(gradknotpot,eigen_vector,"gradient knot potential")
ABORIA_VARIABLE(gradP,double,"gradient_of_potential")
ABORIA_VARIABLE(normal,eigen_vector,"normal")

typedef Particles<std::tuple<boundary,alpha1,alpha2,boundary,target,inlet,outlet,interior,velocity,knotpot,gradknotpot>,2,std::vector,nanoflann_adaptor> KnotsType;
typedef Particles<std::tuple<boundary,dvelocity,dpressure,velocity,pressure,kernel_constant>,2> ComsolType;
typedef Particles<std::tuple<velocity,stokes_velocity,electro_velocity,angle,count,charge>,2> ParticlesType;
typedef Particles<std::tuple<point_a,point_b,traction,boundary,normal,gradP,charge>,2> ElementsType;
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


void read_data_files(ComsolType &particles);


