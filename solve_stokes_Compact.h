#pragma once
#include "filter.h"

    double kernel(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        const double r2 = r*r;
        const double r3 = r2*r;
        const double r4 = r2*r2;
        return std::pow(1.0-r,10)*(429*r4 + 450*r3 + 210*r2 + 50*r + 5); 
    }

    double kernel_x(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        return -26.0*dx[0]*std::pow(invh,2)*std::pow(r-1.0,9)*(5.0 + 3.0*r*(15.0 + r*(53.0+77.0*r))); 
    }
    double kernel_y(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        return -26.0*dx[1]*std::pow(invh,2)*std::pow(r-1.0,9)*(5.0 + 3.0*r*(15.0 + r*(53.0+77.0*r))); 
    }
    double kernel_xx(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        return 26.0*std::pow(invh,2)*std::pow(r-1.0,8)*((-1.0 + r)*(5.0+3.0*r*(15.0+r*(53.0+77.0*r))) 
                                    + 132.0*(1.0+r*(8.0+21.0*r))*dx[0]*dx[0]*invh*invh); 
    }

    double kernel_yy(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        return 26.0*std::pow(invh,2)*std::pow(r-1.0,8)*((-1.0 + r)*(5.0+3.0*r*(15.0+r*(53.0+77.0*r))) 
                                    + 132.0*(1.0+r*(8.0+21.0*r))*dx[1]*dx[1]*invh*invh); 
    }
    double kernel_xy(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        return 3432.0*dx[0]*dx[1]*std::pow(invh,4)*std::pow(r-1.0,8)*(1.0 + r*(8.0 + 21.0*r)); 
    }

    double laplace_xx(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        const double r2 = r*r;
        const double r3 = r2*r;
        const double r4 = r2*r2;
        return 6864.0*std::pow(invh,4)*std::pow(r-1.0,6)*(2.0+12.0*r-3.0*r2-158.0*r3+147.0*r4
                                +30.0*(-3.0+r*(-18.0+49.0*r))*dx[0]*dx[0]*invh*invh);
    };
    double laplace_yy(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        const double r2 = r*r;
        const double r3 = r2*r;
        const double r4 = r2*r2;
        return 6864.0*std::pow(invh,4)*std::pow(r-1.0,6)*(2.0+12.0*r-3.0*r2-158.0*r3+147.0*r4
                                +30.0*(-3.0+r*(-18.0+49.0*r))*dx[1]*dx[1]*invh*invh);
    };
    double laplace_xy(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        return 205920.0*dx[0]*dx[1]*std::pow(invh,6)*std::pow(r-1.0,6)*(-3.0+r*(-18.0+49.0*r));
    }
    double laplace2_xx(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        const double r2 = r*r;
        const double r3 = r2*r;
        const double r4 = r2*r2;
        return 2471040.0*std::pow(invh,6)*std::pow(r-1.0,4)*(-1.0-4.0*r+46.0*r2-90.0*r3+49.0*r4
                            +14.0*(8.0+r*(-31.0+28.0*r))*dx[0]*dx[0]*invh*invh);
    };
    double laplace2_yy(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        const double r2 = r*r;
        const double r3 = r2*r;
        const double r4 = r2*r2;
        return 2471040.0*std::pow(invh,6)*std::pow(r-1.0,4)*(-1.0-4.0*r+46.0*r2-90.0*r3+49.0*r4
                            +14.0*(8.0+r*(-31.0+28.0*r))*dx[1]*dx[1]*invh*invh);
    };
    double laplace2_xy(const double2& dx, const double invh) {
        const double r = dx.norm()*invh;
        if (r > 1) return 0.0;
        return 34594560.0*dx[0]*dx[1]*std::pow(invh,8)*std::pow(r-1.0,4)*(8.0+r*(-31.0+28.0*r));
    };


double solve_stokes_Compact(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0);
