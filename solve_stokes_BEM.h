#pragma once
#include "filter.h"
#include "setup_knots.h"
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/quadrature/gauss.hpp>


void setup_elements(ElementsType& elements, ElementsType& boundary, ParticlesType& fibres, vdouble2 domain_min, vdouble2 domain_max, const unsigned int nx, const double fibre_radius );
 
double solve_stokes_BEM(KnotsType &knots, ElementsType& elements, ElementsType& boundary, const double alpha, const int nlambda, const int nmu);
 
 
auto make_greens_kernel(const double alpha, const int nlambda, const int nmu, const double h, const vdouble2 domain_min, const vdouble2 domain_max, const bool self) {
    typedef Eigen::Matrix<double,2,2> mat2x2;


    const vdouble2 box_size = domain_max-domain_min;
    const double tau = box_size.prod();
    const double inv4alpha = 1.0/(4.0*alpha);
    const double inv4sqrtalpha = 1.0/(4.0*std::sqrt(alpha));
    const double sqrtpialpha = std::sqrt(PI*alpha);

    auto integrate_wave_space = [=](const vdouble2& K,
            const vdouble2& dx_a,
            const vdouble2& dx_b) {
        mat2x2 result;
        const double k2 = K.squaredNorm();
        const vdouble2 xd = dx_b - dx_a;
        const double sin_plus_sin = std::sin(dx_b[0]*K[0]+dx_b[1]*K[1]) 
            - std::sin(dx_a[0]*K[0] + dx_a[1]*K[1]);
        const double xd_dot_K = xd[0]*K[0] + xd[1]*K[1];
        double C;
        if (xd_dot_K == 0.0) {
            C = 4.0*PI*h*(1 + alpha*k2)*std::exp(-alpha*k2)/(tau*std::pow(k2,2));
        } else {
            C = 4.0*PI*h*(1 + alpha*k2)*std::exp(-alpha*k2)*sin_plus_sin/(tau*std::pow(k2,2)*xd_dot_K);
        }
        result(0,0) = C*(-k2 + K[0]*K[0]);
        result(0,1) = C*(      K[0]*K[1]);
        result(1,0) = result(0,1);
        result(1,1) = C*(-k2 + K[1]*K[1]);
        return result;
    };

    auto integrate_real_space0 = [=](
            const vdouble2& dx_a, 
            const vdouble2& dx_b) {
        mat2x2 result;
        vdouble2 d = dx_b - dx_a;
        //const double h = d.norm();
        d /= h;
        const double E1 = boost::math::expint(1,std::pow(h,2)*0.25*inv4alpha);
        const double erfh = boost::math::erf(h*inv4sqrtalpha);
        const double C1 = 0.5*h*E1;
        const double C2 = 2*sqrtpialpha*erfh;
        result(0,0) = C1 + C2*d[0]*d[0];
        result(0,1) =      C2*d[0]*d[1];
        result(1,0) = result(0,1);
        result(1,1) = C1 + C2*d[1]*d[1];
        return result;
    };

    // x = x - Llamba
    auto integrate_real_space = [=](
            const vdouble2& dx_a, 
            const vdouble2& dx_b) {

        mat2x2 result;
        const vdouble2 d = dx_b - dx_a;

        auto fE1 = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            const double E1 = boost::math::expint(1,r2*inv4alpha);
            return E1;
        };

        auto f00 = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            const double exp = std::exp(-r2*inv4alpha);
            return (x[0]*x[0]/r2 - 1)*exp;
        };

        auto f11 = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            const double exp = std::exp(-r2*inv4alpha);
            return (x[1]*x[1]/r2 - 1)*exp;
        };

        auto f01 = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            const double exp = std::exp(-r2*inv4alpha);
            return (x[0]*x[1]/r2    )*exp;
        };


        if (dx_b.squaredNorm() > 0.2*0.2*box_size[0]*box_size[0]) {
            const double E1 = boost::math::quadrature::gauss<double, 1>::
                integrate(fE1, 0.0, 1.0);

            result(0,0) = 0.5*E1 + boost::math::quadrature::gauss<double, 1>::
                integrate(f00, 0.0, 1.0);
            result(0,1) =          boost::math::quadrature::gauss<double, 1>::
                integrate(f01, 0.0, 1.0);
            result(1,0) = result(0,1);
            result(1,1) = 0.5*E1 + boost::math::quadrature::gauss<double, 1>::
                integrate(f11, 0.0, 1.0);
            result *= d.norm();

        } else {
            const double E1 = boost::math::quadrature::gauss<double, 7>::
                integrate(fE1, 0.0, 1.0);

            result(0,0) = 0.5*E1 + boost::math::quadrature::gauss<double, 7>::
                integrate(f00, 0.0, 1.0);
            result(0,1) =          boost::math::quadrature::gauss<double, 7>::
                integrate(f01, 0.0, 1.0);
            result(1,0) = result(0,1);
            result(1,1) = 0.5*E1 + boost::math::quadrature::gauss<double, 7>::
                integrate(f11, 0.0, 1.0);
            result *= d.norm();
        }
        return result;
    };


    auto Aeval = [=](auto a, auto b) {

        //std::cout << "Aeval:  a = "<<get<position>(a)<<" b = "<<get<position>(b)<< " b1 = "<<get<point_a>(b)<<" b2 = "<<get<point_b>(b)<<std::endl;
        mat2x2 result = mat2x2::Zero();

        const vdouble2 dx_a = get<point_a>(b)-get<position>(a);
        const vdouble2 dx_b = get<point_b>(b)-get<position>(a);
        //std::cout << "Aeval: dx_a = "<<dx_a<<"  dx_b = "<<dx_b<< std::endl;

        for (int lambda1 = -nlambda; lambda1 <= nlambda; ++lambda1) {
            for (int lambda2 = -nlambda; lambda2 <= nlambda; ++lambda2) {
                const vdouble2 L = box_size*vdouble2(lambda1,lambda2);
                if (self && lambda1 == 0 && lambda2 == 0 && get<id>(a) == get<id>(b)) {
                    result += integrate_real_space0(dx_a,dx_b);
                } else {
                    result += integrate_real_space(dx_a-L,dx_b-L);
                }
                //std::cout << "result for L = "<<L<<" = "<<result << std::endl;
            }
        }

        for (int mu1 = -nmu; mu1 <= nmu; ++mu1) {
            for (int mu2 = -nmu; mu2 <= nmu; ++mu2) {
                if (mu2 == 0 && mu1 == 0) continue;
                const vdouble2 K = (2*PI/tau)*vdouble2(box_size[1]*mu1,
                        box_size[0]*mu2); 
                result += integrate_wave_space(K,dx_a,dx_b);
                //std::cout << "result for K = "<<K<<" = "<<result << std::endl;
            }
        }
        //std::cout << "result = "<<result << std::endl;
        return result;
    };
    
    return Aeval;
}

auto make_greens_kernel_2d1p(const double alpha, const int nlambda, const int nmu, const double h, const vdouble2 domain_min, const vdouble2 domain_max, const bool self) {
    typedef Eigen::Matrix<double,2,2> mat2x2;

    const vdouble2 box_size = domain_max-domain_min;
    const double k = 2*PI/box_size[0];

    // x = x - Llamba
    auto integrate_real_space = [=](
            const vdouble2& dx_a, 
            const vdouble2& dx_b,
            const bool self) {

        mat2x2 result;
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

        auto fSxx = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            return A(x) + k*x[1]*Ay(x) - 1.0;
        };

        auto fSxy = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            return -k*x[1]*Ax(x);
        };

        auto fSyy = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            return A(x) - k*x[1]*Ay(x);
        };


        auto fSxx_minus_ST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return A(x) + k*x[1]*Ay(x) - 1.0
                 - (std::log(std::sqrt(r2)) - x[0]*x[0]/r2);
        };

        auto fSxy_minus_ST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return -k*x[1]*Ax(x)
                 - (0.0                     - x[0]*x[1]/r2);
        };


        auto fSyy_minus_ST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return A(x) - k*x[1]*Ay(x)
                 - (std::log(std::sqrt(r2)) - x[1]*x[1]/r2);
        };

        auto fSxxST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return (std::log(std::sqrt(r2)) - x[0]*x[0]/r2);
        };

        auto fSxyST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return  (0.0                     - x[0]*x[1]/r2);
        };


        auto fSyyST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return (std::log(std::sqrt(r2)) - x[1]*x[1]/r2);
        };


        const double h = d.norm();

        if (self) {
            const vdouble2 dn = d/h;
            const double C = std::log(h) - 1.69314718055995;
            result(0,0) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSxx_minus_ST, 0.0, 1.0)
                + h*(
                   //-2*dn[0]*dn[0] - dn[1]*dn[1] + C
                   -dn[0]*dn[0] + C
                  );
            result(0,1) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSxy_minus_ST, 0.0, 1.0)
                + h*(
                   -dn[0]*dn[1]
                  );
            result(1,0) = result(0,1);
            result(1,1) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSyy_minus_ST, 0.0, 1.0)
                + h*(
                   //-dn[0]*dn[0] - 2*dn[1]*dn[1] + C
                   -dn[1]*dn[1] + C
                  );
        } else {

            result(0,0) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSxx, 0.0, 1.0);
            result(0,1) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSxy, 0.0, 1.0);
            result(1,0) = result(0,1);
            result(1,1) = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSyy, 0.0, 1.0);

            /*
            if (std::abs(0.5*(dx_a[0]+dx_b[0])) < 0.01) {

                mat2x2 result2;
                result2(0,0) = h*boost::math::quadrature::gauss<double, 8>::
                    integrate(fSxxST, 0.0, 1.0);
                result2(0,1) = h*boost::math::quadrature::gauss<double, 8>::
                    integrate(fSxyST, 0.0, 1.0);
                result2(1,0) = result2(0,1);
                result2(1,1) = h*boost::math::quadrature::gauss<double, 8>::
                    integrate(fSyyST, 0.0, 1.0);

                std::cout << " dx_a = "<<dx_a<<"  dx_b = "<<dx_b<< std::endl;
                std::cout << "result = "<<result << std::endl;
                std::cout << "result2 = "<<result2 << std::endl;
            }
            */
        }

        return result;
    };


    auto Aeval = [=](auto a, auto b) {

        //std::cout << "Aeval:  a = "<<get<position>(a)<<" b = "<<get<position>(b)<< " b1 = "<<get<point_a>(b)<<" b2 = "<<get<point_b>(b)<<std::endl;
        mat2x2 result = mat2x2::Zero();

        const vdouble2 dx_a = get<point_a>(b)-get<position>(a);
        const vdouble2 dx_b = get<point_b>(b)-get<position>(a);

        result = integrate_real_space(dx_a,dx_b,self && get<id>(a)==get<id>(b));
        //result = integrate_real_space(dx_a,dx_b, dx_a[0]<0 != dx_b[0]<0);
        
        return result;
    };
    
    return Aeval;
}

auto make_boundary_layer_kernel_2d1p(const double mu, const vdouble2 domain_min, const vdouble2 domain_max, const bool self) {

    typedef Eigen::Matrix<double,2,1> mat2x1;
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

        auto Axx = [=](const vdouble2& x) {
            const double coskx = std::cos(k*x[0]);
            const double coshky = std::cosh(k*x[1]);

            const double T1 = 0.25*std::pow(std::sin(k*x[0]),2)*coshky;
            const double T2 = 0.25*coskx*std::pow(coshky,2);
            const double T3 = 0.25*coskx;
            const double T4 = -0.5*coshky;
            const double T5 = std::pow(coskx-coshky,3);
            return -(T1+T2+T3+T4)/T5;
        };

        auto Axy = [=](const vdouble2& x) {
            const double coskx = std::cos(k*x[0]);
            const double coshky = std::cosh(k*x[1]);
            const double T5 = std::pow(coskx-coshky,2);
            return -0.5*std::sin(k*x[0])*std::sinh(k*x[1])/T5;
        };

        auto Tyxy = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            return -k*x[1]*Axy(x);
        };

        auto Tyyy = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            return k*x[1]*Axx(x) + Ay(x);
        };

        auto Tyxy_minus_ST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return -k*x[1]*Axy(x)
                        - 2*x[1]*x[0]*x[1]/std::pow(r2,2);
        };

        auto Tyyy_minus_ST = [=](const double t) { 
            const vdouble2 x = d*t + dx_a;
            const double r2 = x.squaredNorm();
            return k*x[1]*Axx(x) + Ay(x)
                    - 2*x[1]*x[1]*x[1]/std::pow(r2,2);
        };

        mat2x1 result;
        const double h = d.norm();

        if (self) {
            const vdouble2 dn = d/h;
            result(0) = 2*mu*h*k*boost::math::quadrature::gauss<double, 8>::
                integrate(Tyxy_minus_ST, 0.0, 1.0); 
            result(1) = 2*mu*h*k*boost::math::quadrature::gauss<double, 8>::
                integrate(Tyyy_minus_ST, 0.0, 1.0); 
        } else {
            const vdouble2 dn = d/h;
            result(0) = 2*mu*h*k*boost::math::quadrature::gauss<double, 8>::
                integrate(Tyxy, 0.0, 1.0); 
            result(1) = 2*mu*h*k*boost::math::quadrature::gauss<double, 8>::
                integrate(Tyyy, 0.0, 1.0); 
        }
        
        return result;
    };

    auto Aeval = [=](auto a, auto b) {
        mat2x1 result = mat2x1::Zero();
        //std::cout << "Aeval:  a = "<<get<position>(a)<<" b = "<<get<position>(b)<< " b1 = "<<get<point_a>(b)<<" b2 = "<<get<point_b>(b)<<std::endl;

        const vdouble2 dx_a = get<point_a>(b)-get<position>(a);
        const vdouble2 dx_b = get<point_b>(b)-get<position>(a);


        if (get<boundary>(b)) {
            result = integrate_real_space(dx_a,dx_b,self && get<id>(a) == get<id>(b));
        }
        //std::cout << "Aeval: boundary = "<<get<boundary>(b)<<" dx_a = "<<dx_a<<"  dx_b = "<<dx_b<< std::endl;
        //std::cout << "result = "<<result << std::endl;
        return result;
    };

    return Aeval;

}

