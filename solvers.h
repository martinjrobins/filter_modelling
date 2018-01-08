#include "filter.h"
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/quadrature/gauss.hpp>


//https://books.google.co.uk/books?id=iXTLBQAAQBAJ&pg=PA207&lpg=PA207&dq=laplace+green%27s+function+singly+periodic&source=bl&ots=UwwLqXz4xb&sig=GvO-3wMVN_wkcg3tm7GmOkfcxTw&hl=en&sa=X&ved=0ahUKEwjBlJi828LYAhUIJMAKHZN6BpoQ6AEIODAC#v=onepage&q=laplace%20green's%20function%20singly%20periodic&f=false 
auto make_laplace_SLP_2d1p(const vdouble2 domain_min, const vdouble2 domain_max, const bool self) {

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

        double result;
        if (self) {
            const vdouble2 dn = d/h;
            const double C = std::log(h) - 1.69314718055995;
            result = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSxx_minus_ST, 0.0, 1.0)
                + h*C;
        } else {
            result = h*boost::math::quadrature::gauss<double, 8>::
                integrate(fSxx, 0.0, 1.0);

        }

        return result;
    };


    auto Aeval = [=](auto a, auto b) {

        //std::cout << "Aeval:  a = "<<get<position>(a)<<" b = "<<get<position>(b)<< " b1 = "<<get<point_a>(b)<<" b2 = "<<get<point_b>(b)<<std::endl;

        const vdouble2 dx_a = get<point_a>(b)-get<position>(a);
        const vdouble2 dx_b = get<point_b>(b)-get<position>(a);

        return integrate_real_space(dx_a,dx_b,self && get<id>(a)==get<id>(b));
        //result = integrate_real_space(dx_a,dx_b, dx_a[0]<0 != dx_b[0]<0);
        
    };
    
    return Aeval;
}

auto make_laplace_gradSLP_2d1p(const vdouble2 domain_min, const vdouble2 domain_max, const bool self) {

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
        //std::cout << "Aeval:  a = "<<get<position>(a)<<" b = "<<get<position>(b)<< " b1 = "<<get<point_a>(b)<<" b2 = "<<get<point_b>(b)<<std::endl;

        const vdouble2 dx_a = get<point_a>(b)-get<position>(a);
        const vdouble2 dx_b = get<point_b>(b)-get<position>(a);



        /*
        std::cout << "Aeval:  k*dx_a = "<<k*dx_a<<"  k*dx_b = "<<k*dx_b<< std::endl;
        std::cout << "position a = "<<get<position>(a) << std::endl;
        std::cout << "pointa = "<<get<point_a>(b)<<" pointb = "<<get<point_b>(b) << std::endl;
        std::cout << "result = "<<integrate_real_space(dx_a,dx_b,self && get<id>(a) == get<id>(b))<< std::endl;
        */
        return integrate_real_space(dx_a,dx_b,self && get<id>(a) == get<id>(b));
    };

    return Aeval;

}

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
    std::cout << "make greens kernel 2d1p k = "<<k << std::endl;

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

        mat2x2 result = mat2x2::Zero();

        const vdouble2 dx_a = get<point_a>(b)-get<position>(a);
        const vdouble2 dx_b = get<point_b>(b)-get<position>(a);

        result = integrate_real_space(dx_a,dx_b,self && get<id>(a)==get<id>(b));

        //std::cout << "Aeval:  a = "<<get<position>(a)<<" b = "<<get<position>(b)<< " b1 = "<<get<point_a>(b)<<" b2 = "<<get<point_b>(b)<<std::endl;
        //std::cout << "result = "<< result << std::endl;
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

void setup_elements(ElementsType& elements, ElementsType& boundarye, ParticlesType& fibres, vdouble2 domain_min, vdouble2 domain_max, const unsigned int nx, const double fibre_radius ) {
    const double L = domain_max[0] - domain_min[0];
    const double Ly = domain_max[1] - domain_min[1];
    const double dtheta = 2*PI/nx;

    elements.resize(fibres.size()*nx);
    vdouble2 ns_buffer(0,L/10.0);

    std::cout << "setup elements: domain_min = "<<domain_min<<" domain_max = "<<domain_max<<" fibre_radius = "<<fibre_radius<<" nx = "<<nx << std::endl;

    

    for (int ii=0; ii<fibres.size(); ++ii) {
        const vdouble2 origin = get<position>(fibres)[ii];
        bool outside,started;
        started = false;
        for (int kk=0; kk<nx; ++kk) {
            ElementsType::reference p = elements[ii*nx+kk];
            get<position>(p) = origin + fibre_radius*vdouble2(std::cos(kk*dtheta),std::sin(kk*dtheta));
            get<point_a>(p) = origin + fibre_radius*vdouble2(std::cos((kk-0.5)*dtheta),std::sin((kk-0.5)*dtheta));
            get<point_b>(p) = origin + fibre_radius*vdouble2(std::cos((kk+0.5)*dtheta),std::sin((kk+0.5)*dtheta));
            get<boundary>(p) = false;
            get<normal>(p) = eigen_vector(std::cos(kk*dtheta),std::sin(kk*dtheta));
                
            for (int i = 0; i < 2; ++i) {
                if (get<position>(p)[i] < domain_min[i]) {
                    get<position>(p)[i] += domain_max[i] - domain_min[i];
                    get<point_a>(p)[i] += domain_max[i] - domain_min[i];
                    get<point_b>(p)[i] += domain_max[i] - domain_min[i];
                } else if ((get<position>(p)[i] >= domain_max[i])) {
                    get<position>(p)[i] -= domain_max[i] - domain_min[i];
                    get<point_a>(p)[i] -= domain_max[i] - domain_min[i];
                    get<point_b>(p)[i] -= domain_max[i] - domain_min[i];
                }
            }
            //std::cout << "element "<<ii*nx+kk<<" has p = "<<get<position>(p)<<" p1 = "<<get<point_a>(p)<<" and p2 = "<<get<point_b>(p) << std::endl;
        }
    }

    /*
    const double dx_aim = (get<position>(elements)[0]-get<position>(elements)[1]).norm();
    const int n_inlet = std::ceil(L/dx_aim);
    const double dx = L/n_inlet;
    //elements.resize(fibres.size()*nx + 2*n_inlet);
    boundarye.resize(2*n_inlet);

    for (int i = 0; i < n_inlet; ++i) {
        //ElementsType::reference p = elements[fibres.size()*nx + 2*i];
        ElementsType::reference p = boundarye[2*i];
        get<position>(p) = vdouble2((i+0.5)*dx,domain_min[1]-2*(domain_max[1]-domain_min[1]));
        get<point_a>(p) = vdouble2((i)*dx,domain_min[1]-2*(domain_max[1]-domain_min[1]));
        get<point_b>(p) = vdouble2((i+1)*dx,domain_min[1]-2*(domain_max[1]-domain_min[1]));
        get<boundary>(p) = true;

        //ElementsType::reference p2 = elements[fibres.size()*nx + 2*i+1];
        ElementsType::reference p2 = boundarye[2*i+1];
        get<position>(p2) = vdouble2((i+0.5)*dx,domain_max[1]+2*(domain_max[1]-domain_min[1]));
        get<point_a>(p2) = vdouble2((i)*dx,domain_max[1]+2*(domain_max[1]-domain_min[1]));
        get<point_b>(p2) = vdouble2((i+1)*dx,domain_max[1]+2*(domain_max[1]-domain_min[1]));
        get<boundary>(p2) = true;
        
    }
    */

    elements.init_neighbour_search(domain_min,domain_max,vbool2(true,false));
    //boundarye.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,vbool2(true,false));
}



double solve_stokes_BEM(KnotsType &knots, ElementsType& elements, ElementsType& boundarye, const double alpha, const int nlambda, const int nmu) {
    std::cout << "solving stokes with BEM..."<<elements.size()<<std::endl;

    const double flow_rate = 1.0;
    const double mu = 1.0;

    typedef typename KnotsType::position position;
    typedef typename position::value_type const & const_position_reference;
    typedef typename KnotsType::const_reference const_particle_reference;

    const vdouble2 min = elements.get_min();
    const vdouble2 max = elements.get_max();
    const double h = (get<point_b>(elements)[0] - get<point_a>(elements)[0]).norm();


    auto kernel = make_greens_kernel_2d1p(alpha,nlambda,nmu,h,min,max,true);
    auto A = create_dense_operator(elements,elements,kernel);
    /*
    auto Abl = create_dense_operator(elements,boundarye,
            make_boundary_layer_kernel_2d1p(mu,min,max,false));
            */

    std::cout << "setup equations..."<<std::endl;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Map<vector_type> map_type;

    const size_t N = elements.size();
    //const size_t Nb = boundarye.size();

    vector_type source(2*N);
    vector_type alphas(2*N);
    //vector_type vel(Nb);
    matrix_type A_eigen(2*N,2*N);

    std::cout << "assemble matrix..."<<std::endl;
    //c[i] = c0;
    A.assemble(A_eigen);
    

    for (int ii=0; ii<N; ++ii) {
        source(2*ii) = 0;
        source(2*ii+1) = flow_rate*4*PI*mu;
    }

    std::cout << "solve w BEM ..."<<std::endl;
    alphas = A_eigen.householderQr().solve(source);
    //alphas = A_eigen.colPivHouseholderQr().solve(source);
    double relative_error = (A_eigen*alphas - source).norm() / source.norm();
    std::cout << "The relative error is:\n" << relative_error << std::endl;
    
    std::cout << "done solve..."<<std::endl;

    for (int ii=0; ii<N; ++ii) {
        get<traction>(elements)[ii][0] = alphas[2*ii];
        get<traction>(elements)[ii][1] = alphas[2*ii+1];
    }

    vtkWriteGrid("BEMelements",0,elements.get_grid(true));

    std::cout << "assemble knot matrix..."<<std::endl;
    matrix_type A_knot(2*knots.size(),2*N);
    auto kernel2 = make_greens_kernel_2d1p(alpha,nlambda,nmu,h,min,max,false);
    auto Aknots = create_dense_operator(knots,elements,kernel2);

    vector_type svel = Aknots*alphas/(4*PI*mu);

    for (int ii=0; ii<knots.size(); ++ii) {
        get<velocity>(knots)[ii][0] = svel[2*ii];
        get<velocity>(knots)[ii][1] = svel[2*ii+1]-flow_rate;
    }

    vtkWriteGrid("BEMknots",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
    return relative_error;
}

double solve_laplace_BEM(KnotsType &knots, ElementsType& elements, const double fibre_charge) {
    std::cout << "solving stokes with BEM..."<<elements.size()<<std::endl;

    const double pot = fibre_charge;
    const double mu = 1.0;

    typedef typename KnotsType::position position;
    typedef typename position::value_type const & const_position_reference;
    typedef typename KnotsType::const_reference const_particle_reference;

    const vdouble2 min = elements.get_min();
    const vdouble2 max = elements.get_max();
    const double h = (get<point_b>(elements)[0] - get<point_a>(elements)[0]).norm();


    auto Aslp = create_dense_operator(elements,elements,
            make_laplace_SLP_2d1p(min,max,true));

    std::cout << "setup equations..."<<std::endl;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Map<vector_type> map_type;

    const size_t N = elements.size();

    vector_type potential(N);
    vector_type source(N);
    vector_type alphas(N);
    matrix_type A_eigen(N,N);

    std::cout << "assemble matrix..."<<std::endl;
    Aslp.assemble(A_eigen);

    potential = vector_type::Ones(N)*pot;
    source = -2*PI*potential;

    std::cout << "solve w BEM ..."<<std::endl;
    alphas = A_eigen.householderQr().solve(source);
    //alphas = A_eigen.colPivHouseholderQr().solve(source);
    double relative_error = (A_eigen*alphas - source).norm() / source.norm();
    std::cout << "The relative error is:\n" << relative_error << std::endl;
    
    std::cout << "done solve..."<<std::endl;

    for (int ii=0; ii<N; ++ii) {
        get<gradP>(elements)[ii] = alphas[ii];
    }

    vtkWriteGrid("BEMelements",0,elements.get_grid(true));

    std::cout << "assemble knot matrix..."<<std::endl;
    auto AknotsSLP = create_dense_operator(knots,elements,
            make_laplace_SLP_2d1p(min,max,false));
    auto AknotsDLP = create_dense_operator(knots,elements,
            make_laplace_gradSLP_2d1p(min,max,false));

    AknotsSLP.get_first_kernel().evaluate(get<knotpot>(knots),get<gradP>(elements));
    AknotsDLP.get_first_kernel().evaluate(get<gradknotpot>(knots),get<gradP>(elements));
    std::transform(std::begin(get<gradknotpot>(knots)),
                   std::end(get<gradknotpot>(knots)),
                    std::begin(get<gradknotpot>(knots)),
                   [](auto& i) { return i/(-2*PI); });
    std::transform(std::begin(get<knotpot>(knots)),
                   std::end(get<knotpot>(knots)),
                   std::begin(get<knotpot>(knots)),
                   [](auto& i) { return i/(-2*PI); });
    
    vtkWriteGrid("BEMknots",0,knots.get_grid(true));

    std::cout << "done solving laplace"<<std::endl;
    return relative_error;
}
