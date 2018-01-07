#pragma once
#include "filter.h"
#include "setup_knots.h"
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/quadrature/gauss.hpp>


double solve_laplace_BEM(KnotsType &knots, ElementsType& elements);
 
//https://books.google.co.uk/books?id=iXTLBQAAQBAJ&pg=PA207&lpg=PA207&dq=laplace+green%27s+function+singly+periodic&source=bl&ots=UwwLqXz4xb&sig=GvO-3wMVN_wkcg3tm7GmOkfcxTw&hl=en&sa=X&ved=0ahUKEwjBlJi828LYAhUIJMAKHZN6BpoQ6AEIODAC#v=onepage&q=laplace%20green's%20function%20singly%20periodic&f=false 
auto make_laplace_SLP_2d1p(const vdouble2 domain_min, const vdouble2 domain_max, const bool self) {
    typedef Eigen::Matrix<double,1,1> mat1x1;

    const vdouble2 box_size = domain_max-domain_min;
    const double k = 2*PI/box_size[0];

    // x = x - Llamba
    auto integrate_real_space = [=](
            const vdouble2& dx_a, 
            const vdouble2& dx_b,
            const bool self) {

        mat1x1 result;
        const vdouble2 d = dx_b - dx_a;

        auto A = [=](const vdouble2& x) {
            return 0.5*std::log(2*(std::cosh(k*x[1]) - std::cos(k*x[0])));
        };

        auto fSxx = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            return A(x);
        };

        auto fSxx_minus_ST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return A(x)
                 - std::log(std::sqrt(r2));
        };

        const double h = d.norm();

        if (self) {
            const vdouble2 dn = d/h;
            const double C = std::log(h) - 1.69314718055995;
            result(0,0) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSxx_minus_ST, 0.0, 1.0)
                + h*C;
        } else {
            result(0,0) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSxx, 0.0, 1.0);

        }

        return result;
    };


    auto Aeval = [=](auto a, auto b) {

        //std::cout << "Aeval:  a = "<<get<position>(a)<<" b = "<<get<position>(b)<< " b1 = "<<get<point_a>(b)<<" b2 = "<<get<point_b>(b)<<std::endl;
        mat1x1 result = mat1x1::Zero();

        const vdouble2 dx_a = get<point_a>(b)-get<position>(a);
        const vdouble2 dx_b = get<point_b>(b)-get<position>(a);

        result = integrate_real_space(dx_a,dx_b,self && get<id>(a)==get<id>(b));
        //result = integrate_real_space(dx_a,dx_b, dx_a[0]<0 != dx_b[0]<0);
        
        return result;
    };
    
    return Aeval;
}

auto make_laplace_DLP_2d1p(const vdouble2 domain_min, const vdouble2 domain_max, const bool self) {

    typedef Eigen::Vector2d mat2x1;
    const vdouble2 box_size = domain_max-domain_min;
    const double k = 2*PI/box_size[0];

    // x = x - Llamba
    auto integrate_real_space = [=](
            const vdouble2& dx_a, 
            const vdouble2& dx_b,
            const bool self) {

        const vdouble2 d = dx_b - dx_a;

        auto A = [=](const vdouble2& x) {
            return 0.5*std::log(2*(std::cosh(k*x[1]) - std::cos(k*x[0])));
        };

        auto Ay = [=](const vdouble2& x) {
            return -0.5*std::sinh(k*x[1])/(std::cos(k*x[0])-std::cosh(k*x[1]));
        };

        auto Ax = [=](const vdouble2& x) {
            return -0.5*std::sin(k*x[0])/(std::cos(k*x[0])-std::cosh(k*x[1]));
        };

        auto Tx = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            return k*Ax(x);
        };

        auto Ty = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            return k*Ay(x);
        };

        auto Tx_minus_ST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return k*Ax(x)
                        - x[0]/r2;
        };

        auto Ty_minus_ST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return k*Ay(x)
                        - x[1]/r2;
        };

        mat2x1 result;
        const double h = d.norm();

        if (self) {
            const vdouble2 dn = d/h;
            result(0) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(Tx_minus_ST, 0.0, 1.0); 
            result(1) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(Ty_minus_ST, 0.0, 1.0); 
        } else {
            const vdouble2 dn = d/h;
            result(0) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(Tx, 0.0, 1.0); 
            result(1) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(Ty, 0.0, 1.0); 
        }

        return result;
    };

    auto Aeval = [=](auto a, auto b) {
        mat2x1 result = mat2x1::Zero();
        //std::cout << "Aeval:  a = "<<get<position>(a)<<" b = "<<get<position>(b)<< " b1 = "<<get<point_a>(b)<<" b2 = "<<get<point_b>(b)<<std::endl;

        const vdouble2 dx_a = get<point_a>(b)-get<position>(a);
        const vdouble2 dx_b = get<point_b>(b)-get<position>(a);


        result = integrate_real_space(dx_a,dx_b,self && get<id>(a) == get<id>(b));
        //std::cout << "Aeval: boundary = "<<get<boundary>(b)<<" dx_a = "<<dx_a<<"  dx_b = "<<dx_b<< std::endl;
        //std::cout << "result = "<<result << std::endl;
        return result.dot(get<normal>(b));
    };

    return Aeval;

}

