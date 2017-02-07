template <typename I, typename J, typename C>
auto gen_kernel_mq(I &i, J &j, C &c) {
    auto dx = create_dx(i,j);
    return deep_copy(
        sqrt(dot(dx,dx)+c[i]*c[i])
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

    const double mu = 1.0;
    return deep_copy(
            //-(1.0/mu)*phi_sol_dash_dash_dash
        if_else(norm(dx)==0
            ,-(phi_sol_dash_dash_dash)
            //,0.0
            ,(phi_sol_dash_dash_dash + (phi_sol_dash_dash - phi_sol_dash_div_r)/norm(dx))*dx[0]/norm(dx)
            //,(phi_sol_dash_dash_dash)*dx[0]/norm(dx)
            )
        );
}

template <typename I, typename J, typename C>
auto gen_psol_p2(I &i, J &j, C &c) { 
    auto dx = create_dx(i,j);
    auto phi_sol_dash_div_r = gen_phi_sol_dash_div_r(i,j,c);
    auto phi_sol_dash_dash = gen_phi_sol_dash_dash(i,j,c);
    auto phi_sol_dash_dash_dash = gen_phi_sol_dash_dash_dash(i,j,c);

    const double mu = 1.0;
    return deep_copy(
            //-(1.0/mu)*phi_sol_dash_dash_dash
        if_else(norm(dx)==0
            ,-(phi_sol_dash_dash_dash)
            //,0.0
            ,(phi_sol_dash_dash_dash + (phi_sol_dash_dash - phi_sol_dash_div_r)/norm(dx))*dx[1]/norm(dx)
            //,(phi_sol_dash_dash_dash)*dx[1]/norm(dx)
            )
        );
}

