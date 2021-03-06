
template <typename Dx, typename C>
double kernel_mq(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    return std::sqrt(dx2 + std::pow(c,2));
}


template <typename Dx, typename C>
double phi_sol(const Dx& dx, const C& c) {
        const double c2 = c*c;
        const double c3 = c2*c;
        const double dx2 = dx.squaredNorm();
        const double mu = 1;
        return ((1.0/75.0)*std::sqrt(dx2+c2)*(4.0*dx2*dx2+48.0*dx2*c2-61.0*c2*c2) - c3*log(c)*dx2 - (1.0/5.0)*(5.0*dx2-2*c2)*c3*log(c+std::sqrt(dx2+c2)))/(12.0*mu);
}

/*
template <typename Dx, typename C>
double phi_sol_dash_div_r(const Dx& dx, const A& a, const B& b) {
        const double c = get<kernel_constant>(b);
        const double c2 = c*c;
        const double c3 = c2*c;
        const double dx2 = dx.squaredNorm();
        return (4.0*dx2*dx2*dx2 + 36.0*dx2*dx2*c2 + 39.0*dx2*c2*c2 + 7.0*c3*c3)/(180.0*std::pow(dx2+c2,1.5))
            - (5.0*dx2*dx2 + 3.0*dx2*c2 - 2.0*c2*c2)*c3/(60.0*std::pow(dx2+c2,1.5)*(c+std::sqrt(dx2+c2)))
            - (1.0/6.0)*c3*std::log(c+std::sqrt(dx2+c2)) - (1.0/6.0)*c3*std::log(c);

}
*/

template <typename Dx, typename C>
double phi_sol_dash_div_r(const Dx& dx, const C& c) {
        const double c2 = std::pow(c,2);
        const double c3 = std::pow(c,3);
        const double c4 = std::pow(c,4);
        const double dx2 = dx.squaredNorm();
        const double dx4 = std::pow(dx2,2);
        const double logc = std::log(c);
        const double c_p_mq = std::sqrt(dx2+c2);
        const double logc_p_mq = std::log(c+c_p_mq);
        return (-30*c4*logc - 30*c4*logc_p_mq + 7*c4 - 30*c3*c_p_mq*logc - 30*c3*c_p_mq*logc_p_mq + 13*c3*c_p_mq + 32*c2*dx2 + 4*c*dx2*c_p_mq + 4*dx4)/(180*(c+c_p_mq));

}

template <typename Dx, typename C>
double phi_sol_dash_div_r_fast(const Dx& dx, const C& c) {
        const double c2 = c*c;
        const double c3 = c*c2;
        const double c4 = c2*c2;
        const double dx2 = dx.squaredNorm();
        const double dx4 = dx2*dx2;
        const double logc = std::log(c);
        const double c_p_mq = std::sqrt(dx2+c2);
        const double logc_p_mq = std::log(c+c_p_mq);
        return (-30*c4*logc - 30*c4*logc_p_mq + 7*c4 - 30*c3*c_p_mq*logc - 30*c3*c_p_mq*logc_p_mq + 13*c3*c_p_mq + 32*c2*dx2 + 4*c*dx2*c_p_mq + 4*dx4)/(180*(c+c_p_mq));

}


/*
template <typename Dx, typename C>
double phi_sol_dash_dash(const Dx& dx, const A& a, const B& b) {
        const double c = get<kernel_constant>(b);
        const double c2 = c*c;
        const double c3 = c2*c;
        const double dx2 = dx.squaredNorm();
        const double c_plus_sqrt = c + std::sqrt(dx2+c2);
        return (16.0*dx2*dx2*dx2+84.0*dx2*dx2*c2+96.0*dx2*c2*c2+7.0*c3*c3)/(180.0*std::pow(dx2+c2,1.5))
            - (20.0*dx2*dx2+25.0*dx2*c2-2.0*c2*c2)*c3/(60.0*std::pow(dx2+c2,1.5)*c_plus_sqrt)
            + (5.0*dx2-2.0*c2)*c3*dx2/(60.0*(dx2+c2)*c_plus_sqrt*c_plus_sqrt)
            - (1.0/6.0)*c3*std::log(c_plus_sqrt) - (1.0/6.0)*c3*std::log(c);
}
*/

template <typename Dx, typename C>
double phi_sol_dash_dash(const Dx& dx, const C& c) {
        const double c2 = std::pow(c,2);
        const double c3 = std::pow(c,3);
        const double c4 = std::pow(c,4);
        const double c5 = std::pow(c,5);
        const double dx2 = dx.squaredNorm();
        const double dx4 = std::pow(dx2,2);
        const double logc = std::log(c);
        const double c_p_mq = std::sqrt(dx2+c2);
        const double logc_p_mq = std::log(c+c_p_mq);
        return (-60*c5*logc - 60*c5*logc_p_mq + 20*c5 - 60*c4*c_p_mq*logc - 60*c4*c_p_mq*logc_p_mq + 20*c4*c_p_mq - 30*c3*dx2*logc - 30*c3*dx2*logc_p_mq + 91*c3*dx2 + 84*c2*dx2*c_p_mq + 32*c*dx4 + 16*dx4*c_p_mq)/(180*(2*c2+2*c*c_p_mq+dx2));

}

template <typename Dx, typename C>
double phi_sol_dash_dash_fast(const Dx& dx, const C& c) {
        const double c2 = c*c;
        const double c3 = c*c2;
        const double c4 = c2*c2;
        const double c5 = c2*c3;
        const double dx2 = dx.squaredNorm();
        const double dx4 = dx2*dx2;
        const double logc = std::log(c);
        const double c_p_mq = std::sqrt(dx2+c2);
        const double logc_p_mq = std::log(c+c_p_mq);
        return (-60*c5*logc - 60*c5*logc_p_mq + 20*c5 - 60*c4*c_p_mq*logc - 60*c4*c_p_mq*logc_p_mq + 20*c4*c_p_mq - 30*c3*dx2*logc - 30*c3*dx2*logc_p_mq + 91*c3*dx2 + 84*c2*dx2*c_p_mq + 32*c*dx4 + 16*dx4*c_p_mq)/(180*(2*c2+2*c*c_p_mq+dx2));

}


template <typename Dx, typename C>
double phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(const Dx& dx, const C& c) {
        const double c2 = std::pow(c,2);
        const double c3 = std::pow(c,3);
        const double dx2 = dx.squaredNorm();
        const double nx = std::sqrt(dx2);
        return -nx*(7.0*c3 + 8.0*c2*std::sqrt(dx2+c2) + 4.0*dx2*c + 2.0*dx2*std::sqrt(dx2+c2))/(60.0*c2 + 60.0*c*std::sqrt(dx2+c2) + 30.0*dx2);
}

template <typename Dx, typename C>
double phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r_fast(const Dx& dx, const C& c) {
        const double c2 = c*c;
        const double c3 = c2*c;
        const double dx2 = dx.squaredNorm();
        const double nx = std::sqrt(dx2);
        return -nx*(7.0*c3 + 8.0*c2*std::sqrt(dx2+c2) + 4.0*dx2*c + 2.0*dx2*std::sqrt(dx2+c2))/(60.0*c2 + 60.0*c*std::sqrt(dx2+c2) + 30.0*dx2);
}

/*
template <typename Dx, typename C>
double phi_sol_dash_dash_dash(const Dx& dx, const A& a, const B& b) {
        const double c = get<kernel_constant>(b);
        const double c2 = c*c;
        const double c3 = c2*c;
        const double dx2 = dx.squaredNorm();
        const double ndx = std::sqrt(dx2);
        const double c_plus_sqrt = c + std::sqrt(dx2+c2);

        return (76.0*dx2*dx2+176.0*dx2*c2+285.0*c2*c2)*ndx/(300.0*std::pow(dx2+c2,1.5))
            + (4.0*dx2*dx2+48.0*dx2*c2-61.0*c2*c2)*dx2*ndx/(300.0*std::pow(dx2+c2,5.0/2.0))
            + (10.0*dx2*dx2+15.0*dx2*c2-2.0*c2*c2)*c3*ndx/(20.0*(dx2+c2)*(dx2+c2)*c_plus_sqrt)
            - (5.0*dx2+22.0*c2)*c3*ndx/(20.0*std::pow(dx2+c2,1.5)*c_plus_sqrt)
            + (-5.0*dx2+2*c2)*c3*ndx*dx2/(20.0*std::pow(dx2+c2,5.0/2.0)*c_plus_sqrt)
            + (-5.0*dx2+2*c2)*c3*ndx*dx2/(30.0*std::pow(dx2+c2,1.5)*c_plus_sqrt*c_plus_sqrt*c_plus_sqrt);
}
*/

template <typename Dx, typename C>
double phi_sol_dash_dash_dash(const Dx& dx, const C& c) {
        const double c2 = std::pow(c,2);
        const double c3 = std::pow(c,3);
        const double c4 = std::pow(c,4);
        const double c5 = std::pow(c,5);
        const double dx2 = dx.squaredNorm();
        const double dx4 = std::pow(dx2,2);
        const double r = std::sqrt(dx2);
        const double logc = std::log(c);
        const double c_p_mq = std::sqrt(dx2+c2);
        const double logc_p_mq = std::log(c+c_p_mq);
        return r*(45*c5 + 45*c4*c_p_mq + 69*c3*dx2 + 46*c2*dx2*c_p_mq + 24*c*dx4 + 8*dx4*c_p_mq)/(30*(4*c4 + 4*c3*c_p_mq + 5*c2*dx2 + 3*c*dx2*c_p_mq + dx4));
}

template <typename Dx, typename C>
double phi_sol_dash_dash_dash_fast(const Dx& dx, const C& c) {
        const double c2 = c*c;
        const double c3 = c*c2;
        const double c4 = c2*c2;
        const double c5 = c2*c3;
        const double dx2 = dx.squaredNorm();
        const double dx4 = dx2*dx2;
        const double r = std::sqrt(dx2);
        const double logc = std::log(c);
        const double c_p_mq = std::sqrt(dx2+c2);
        const double logc_p_mq = std::log(c+c_p_mq);
        return r*(45*c5 + 45*c4*c_p_mq + 69*c3*dx2 + 46*c2*dx2*c_p_mq + 24*c*dx4 + 8*dx4*c_p_mq)/(30*(4*c4 + 4*c3*c_p_mq + 5*c2*dx2 + 3*c*dx2*c_p_mq + dx4));
}


template <typename Dx, typename C>
double psol_u1(const Dx& dx, const C& c) {
    const double vphi_sol_dash_div_r = phi_sol_dash_div_r(dx,c);
    const double vphi_sol_dash_dash = phi_sol_dash_dash(dx,c);
    const double dx2 = dx.squaredNorm();

    const double mu = 1.0;
    // i = 1, l = 1
    if (dx2==0) {
        return (0.5/mu)*(vphi_sol_dash_div_r + vphi_sol_dash_dash);
    } else {
        return (1.0/mu)*(vphi_sol_dash_div_r*dx[0]*dx[0] + vphi_sol_dash_dash*dx[1]*dx[1])/dx2;
    }
}

template <typename Dx, typename C>
double psol_u2(const Dx& dx, const C& c) {
    const double vphi_sol_dash_div_r = phi_sol_dash_div_r(dx,c);
    const double vphi_sol_dash_dash = phi_sol_dash_dash(dx,c);
    const double dx2 = dx.squaredNorm();

    const double mu = 1.0;
    // i = 1, l = 2
    if (dx2==0) {
        return (1.0/mu)*(vphi_sol_dash_div_r - vphi_sol_dash_dash);
    } else {
        return (1.0/mu)*(vphi_sol_dash_div_r - vphi_sol_dash_dash)*(dx[0]*dx[1])/dx2;
    }
}

template <typename Dx, typename C>
double psol_v1(const Dx& dx, const C& c) {
    const double vphi_sol_dash_div_r = phi_sol_dash_div_r(dx,c);
    const double vphi_sol_dash_dash = phi_sol_dash_dash(dx,c);
    const double dx2 = dx.squaredNorm();

    const double mu = 1.0;
    // i = 2, l = 1
    if (dx2==0) {
        return (1.0/mu)*(vphi_sol_dash_div_r - vphi_sol_dash_dash);
    } else {
        return (1.0/mu)*(vphi_sol_dash_div_r - vphi_sol_dash_dash)*(dx[0]*dx[1])/dx2;
    }
}

template <typename Dx, typename C>
double psol_v2(const Dx& dx, const C& c) {
    const double vphi_sol_dash_div_r = phi_sol_dash_div_r(dx,c);
    const double vphi_sol_dash_dash = phi_sol_dash_dash(dx,c);
    const double dx2 = dx.squaredNorm();

    const double mu = 1.0;
    // i = 1, l = 1
    if (dx2==0) {
        return (0.5/mu)*(vphi_sol_dash_div_r + vphi_sol_dash_dash);
    } else {
        return (1.0/mu)*(vphi_sol_dash_div_r*dx[1]*dx[1] + vphi_sol_dash_dash*dx[0]*dx[0])/dx2;
    }
}



template <typename Dx, typename C>
double psol_du1dx(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double x = dx[0];
    const double x2 = x*x;
    const double x4 = x2*x2;
    const double y = dx[1];
    const double y2 = y*y;
    const double y4 = y2*y2;
    const double mu = 1.0;
    if (dx2==0) {
        return phi_sol_dash_dash_dash(dx,c); 
    } else {
        return -x*(phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c)*(2*y2-x2)/(ndx*dx2) + phi_sol_dash_dash_dash(dx,c)*y2/(dx2*ndx));
    }

}

template <typename Dx, typename C>
double psol_du2dx(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double x = dx[0];
    const double x2 = x*x;
    const double x4 = x2*x2;
    const double y = dx[1];
    const double y2 = y*y;
    const double y4 = y2*y2;
    const double mu = 1.0;
    if (dx2==0) {
        return -phi_sol_dash_dash_dash(dx,c); 
    } else {
        return -y*(phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c)*(y2-2*x2)/(ndx*dx2)- phi_sol_dash_dash_dash(dx,c)*x2/(dx2*ndx));
    }

}

template <typename Dx, typename C>
double psol_dv1dx(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double x = dx[0];
    const double x2 = x*x;
    const double x4 = x2*x2;
    const double y = dx[1];
    const double y2 = y*y;
    const double y4 = y2*y2;
    const double mu = 1.0;
    if (dx2==0) {
        return -phi_sol_dash_dash_dash(dx,c); 
    } else {
        return -y*(phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c)*(y2-2*x2)/(ndx*dx2)- phi_sol_dash_dash_dash(dx,c)*x2/(dx2*ndx));
    }

}

template <typename Dx, typename C>
double psol_dv2dx(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double x = dx[0];
    const double x2 = x*x;
    const double x4 = x2*x2;
    const double y = dx[1];
    const double y2 = y*y;
    const double y4 = y2*y2;
    const double mu = 1.0;
    if (dx2==0) {
        return phi_sol_dash_dash_dash(dx,c);
    } else {
        return -x*(-phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c)*3*y2/(ndx*dx2) + phi_sol_dash_dash_dash(dx,c)*x2/(dx2*ndx));
    }

}

template <typename Dx, typename C>
double psol_du1dy(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double x = dx[0];
    const double x2 = x*x;
    const double x4 = x2*x2;
    const double y = dx[1];
    const double y2 = y*y;
    const double y4 = y2*y2;
    const double mu = 1.0;
    if (dx2==0) {
        return phi_sol_dash_dash_dash(dx,c); 
    } else {
        return -y*(-phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c)*3*x2/(ndx*dx2) + phi_sol_dash_dash_dash(dx,c)*y2/(dx2*ndx));
    }

}


template <typename Dx, typename C>
double psol_du2dy(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double x = dx[0];
    const double x2 = x*x;
    const double x4 = x2*x2;
    const double y = dx[1];
    const double y2 = y*y;
    const double y4 = y2*y2;
    const double mu = 1.0;
    if (dx2==0) {
        return -phi_sol_dash_dash_dash(dx,c); 
    } else {
        return x*(phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c)*(2*y2-x2)/(ndx*dx2) + phi_sol_dash_dash_dash(dx,c)*y2/(dx2*ndx));
    }
}

template <typename Dx, typename C>
double psol_dv1dy(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double x = dx[0];
    const double x2 = x*x;
    const double x4 = x2*x2;
    const double y = dx[1];
    const double y2 = y*y;
    const double y4 = y2*y2;
    const double mu = 1.0;
    if (dx2==0) {
        return -phi_sol_dash_dash_dash(dx,c); 
    } else {
        return x*(phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c)*(2*y2-x2)/(ndx*dx2) + phi_sol_dash_dash_dash(dx,c)*y2/(dx2*ndx));
    }
}

template <typename Dx, typename C>
double psol_dv2dy(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double x = dx[0];
    const double x2 = x*x;
    const double x4 = x2*x2;
    const double y = dx[1];
    const double y2 = y*y;
    const double y4 = y2*y2;
    const double mu = 1.0;
    if (dx2==0) {
        return phi_sol_dash_dash_dash(dx,c); 
    } else {
        return -y*(phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c)*(2*x2-y2)/(ndx*dx2) + phi_sol_dash_dash_dash(dx,c)*x2/(dx2*ndx));
    }

}

template <typename Dx, typename C>
double psol_p1(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double mu = 1.0;
    if (dx2==0) {
        return -phi_sol_dash_dash_dash(dx,c);
    } else {
        return (phi_sol_dash_dash_dash(dx,c) - phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c))*dx[0]/ndx;
        //return (phi_sol_dash_dash_dash(dx,a,b))*dx[0]/ndx;
    }
}

template <typename Dx, typename C>
double psol_p2(const Dx& dx, const C& c) {
    const double dx2 = dx.squaredNorm();
    const double ndx = std::sqrt(dx2);
    const double mu = 1.0;
    if (dx2==0) {
        return -phi_sol_dash_dash_dash(dx,c);
    } else {
        return (phi_sol_dash_dash_dash(dx,c) - phi_sol_dash_div_r2_minus_phi_sol_dash_dash_div_r(dx,c))*dx[1]/ndx;
        //return (phi_sol_dash_dash_dash(dx,a,b))*dx[1]/ndx;
    }
}


