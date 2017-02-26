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

ABORIA_VARIABLE(inlet,uint8_t,"is inlet knot")
ABORIA_VARIABLE(outlet,uint8_t,"is outlet knot")
ABORIA_VARIABLE(boundary,uint8_t,"is boundary knot")
ABORIA_VARIABLE(target,uint8_t,"is target knot")
ABORIA_VARIABLE(interior,uint8_t,"is interior knot")
ABORIA_VARIABLE(velocity_u,double,"velocity u");
ABORIA_VARIABLE(velocity_v,double,"velocity v");
ABORIA_VARIABLE(velocity_dudx,double,"velocity dudx");
ABORIA_VARIABLE(velocity_dudy,double,"velocity dudy");
ABORIA_VARIABLE(velocity_dvdx,double,"velocity dvdx");
ABORIA_VARIABLE(velocity_dvdy,double,"velocity dvdy");
ABORIA_VARIABLE(pressure,double,"pressure");
ABORIA_VARIABLE(dvelocity_u,double,"error velocity u");
ABORIA_VARIABLE(dvelocity_v,double,"error velocity v");
ABORIA_VARIABLE(dpressure,double,"error comsol pressure");
ABORIA_VARIABLE(alpha1,double,"alpha 1");
ABORIA_VARIABLE(alpha2,double,"alpha 2");
ABORIA_VARIABLE(kernel_constant,double,"kernel constant")

typedef Particles<std::tuple<alpha1,alpha2,boundary,target,inlet,outlet,interior,velocity_u,velocity_v,velocity_dudx,velocity_dudy,velocity_dvdx,velocity_dvdy,pressure,kernel_constant>,2> KnotsType;
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


void read_data_files(ComsolType &particles);
