#include "filter.h"
void read_data_files(ComsolType &particles) {

    typedef typename ComsolType::position position;

    std::cout << "reading files..." << std::endl;
    std::ifstream pressure_file("five layer_square/[p]_five_layer_square.txt" );
    std::ifstream vel_horz_file("five layer_square/[vel_horz]_five_layer_square.txt" );
    std::ifstream vel_vert_file("five layer_square/[vel_vert]_five_layer_square.txt" );
    std::string line;
    for (int i=0; i<8; ++i) {
        std::getline(pressure_file, line);
        std::getline(vel_horz_file, line);
        std::getline(vel_vert_file, line);
    }

    int i=0;
    while ( pressure_file.good() ) {
        double2 pos,velocity;
        double pressure,dummy;
        std::getline(pressure_file, line);
        std::istringstream buffer(line);
        buffer >> pos[0];
        buffer >> pos[1];
        buffer >> pressure;
        buffer.clear();

        std::getline(vel_horz_file, line);
        buffer.str(line);
        buffer >> dummy;
        buffer >> dummy;
        buffer >> velocity[0];
        buffer.clear();

        std::getline(vel_vert_file, line);
        buffer.str(line);
        buffer >> dummy;
        buffer >> dummy;
        buffer >> velocity[1];
        buffer.clear();

        typename ComsolType::value_type p;
        get<position>(p) = pos;
        get<dvelocity_u>(p) = velocity[0];
        get<dvelocity_v>(p) = velocity[1];
        get<dpressure>(p) = pressure;
        if (i++ % 10 == 0) {
            particles.push_back(p);
            //std::cout << "position = "<<pos<<std::endl;
        }
    }
    std::cout << "done reading files"<< std::endl;
}

void setup_knots(KnotsType &knots, ParticlesType &fibres, const double fibre_radius, const double fibre_resolution_factor, const double nx, double2 domain_min, double2 domain_max, const double c0, const double k, const double epsilon_strength, const double epsilon_falloff) {

    std::cout << "setup knots..." << std::endl;

    typedef typename KnotsType::position position;
    const double L = domain_max[0] - domain_min[0];
    const double delta = L/nx;
    const double fibre_resolution = delta*fibre_resolution_factor;
    const double s = 2.0*delta;

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
        for (auto tpl: box_search(fibres.get_query(),get<position>(p))) {
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
    Label<1,ParticlesType> bf(fibres);
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

    std::cout << "finshed adapting. Writing to vtk file init_knots"<<std::endl;

    vtkWriteGrid("init_knots",1,knots.get_grid(true));

    std::cout << "finshed writing to file."<<std::endl;
}

void calculate_c(KnotsType &knots, double c0, const double nx, double2 domain_min, double2 domain_max) {
    std::cout << "calculate c..."<<knots.size()<<std::endl;

    typedef typename KnotsType::position position;

    const double L = domain_max[0] - domain_min[0];
    const double delta = L/nx;
    Symbol<kernel_constant> c;

    Label<0,KnotsType> i(knots);
    Label<1,KnotsType> j(knots);
    auto dx = create_dx(i,j);
    Accumulate<std::plus<double> > sum;

    const double mult = 2.0;
    c[i] = mult*delta;
    const double scale = 7.0/(64.0*PI);
    auto kernel_wendland = deep_copy(
        scale*pow(2.0-norm(dx)/c[i],4)*(1.0+2*norm(dx)/c[i])/(c[i]*c[i])
        );

    for (int ii=0; ii<10; ++ii) {
        c[i] = mult*sqrt(1.0/sum(j,norm(dx)<2*c[i],kernel_wendland));
    }

    c[i] = c0*c[i]/mult;

    vtkWriteGrid("init_knots",2,knots.get_grid(true));
    std::cout << "done calculate c."<<std::endl;
}


 

void solve_stokes_MAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in) {
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
    typedef Eigen::SparseMatrix<double> sparse_matrix_type; 
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

    std::cout << "solve w MAPS ..."<<std::endl;
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

    /*
    auto kernel_wendland = deep_copy(
        pow(c[j]-norm(dx),4)
            *(16.0*h[j] 
                + 64.0*norm(dx))/pow(c[j],5)
        );


    auto B = create_eigen_operator(i,j,
               ,kernel_wendland
               ,norm(dx) < h[j] 
            );
    
    source.resize(knots.size());
    alphas.resize(knots.size());
    sparse_matrix_type B_eigen(knots.size(),knots.size());
    std::cout << "assemble matrix..."<<std::endl;
    B.assemble(B_eigen);

    //
    // interpolate u
    //
    for (int ii=0; ii<knots.size(); ++ii) {
        source[ii] = vu[i];
        alphas[ii] = 0.0;
    }
    std::cout << "solve interpolant u ..."<<std::endl;
    solve(B_eigen,alphas,source,
                max_iter_linear,restart_linear,(linear_solver)solver_in);
    std::cout << "done solve..."<<std::endl;

    for (int ii=0; ii<knots.size(); ++ii) {
        get<alpha>(knots)[ii][0] = alphas[ii];
    }

    //
    // interpolate v
    //
    for (int ii=0; ii<knots.size(); ++ii) {
        source[ii] = vv[i];
        alphas[ii] = 0.0;
    }
    std::cout << "solve interpolant v ..."<<std::endl;
    solve(B_eigen,alphas,source,
                max_iter_linear,restart_linear,(linear_solver)solver_in);
    std::cout << "done solve..."<<std::endl;

    for (int ii=0; ii<knots.size(); ++ii) {
        get<alpha>(knots)[ii][1] = alphas[ii];
    }


    //
    // interpolate p
    //
    for (int ii=0; ii<knots.size(); ++ii) {
        source[ii] = pr[i];
        alphas[ii] = 0.0;
    }
    std::cout << "solve interpolant p ..."<<std::endl;
    solve(B_eigen,alphas,source,
                max_iter_linear,restart_linear,(linear_solver)solver_in);
    std::cout << "done solve..."<<std::endl;

    for (int ii=0; ii<knots.size(); ++ii) {
        get<alpha>(knots)[ii][2] = alphas[ii];
    }

    vu[i] = sum(j,true,psol_u1*al[j][0] + psol_u2*al[j][1]);
    vv[i] = sum(j,true,psol_v1*al[j][0] + psol_v2*al[j][1]);
    pr[i] = sum(j,true,psol_p1*al[j][0] + psol_p2*al[j][1]);
    */
    
    vtkWriteGrid("MAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
}

void solve_stokes_LMAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double delta) {
    std::cout << "solving stokes..."<<knots.size()<<std::endl;

    const double flow_rate = 1.0; 
    const double region_factor = 3;

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

            , norm(dx) < delta*region_factor
            );
   
    auto A12 = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,0.0
               ,if_else(is_out[i]
                   ,psol_p2
                   ,psol_u2
                   )
               )
            , norm(dx) < delta*region_factor
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
            , norm(dx) < delta*region_factor
            );

    auto A22 = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,kernel_mq
               ,if_else(is_out[i]
                   ,psol_u2
                   ,psol_v2
                   )
               )
            , norm(dx) < delta*region_factor
            );
    
    auto A = create_block_eigen_operator<2,2>(
            A11,A12,
            A21,A22
            );


    std::cout << "setup equations..."<<std::endl;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type; 
    typedef Eigen::SparseMatrix<double> sparse_matrix_type; 
    typedef Eigen::Map<vector_type> map_type;
    map_type eigen_pressure(get<pressure>(knots).data(),knots.size());
    map_type eigen_velocity_u(get<velocity_u>(knots).data(),knots.size());
    map_type eigen_velocity_v(get<velocity_u>(knots).data(),knots.size());
    

    vector_type source(2*knots.size());
    vector_type alphas(2*knots.size());
    vector_type alphas_test(2*knots.size());
    sparse_matrix_type A_eigen(2*knots.size(),2*knots.size());
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

    std::cout << "solve w LMAPS ..."<<std::endl;
    solve(A_eigen,alphas,source,
                max_iter_linear,restart_linear,(linear_solver)solver_in);


    std::cout << "done solve..."<<std::endl;

    for (int ii=0; ii<knots.size(); ++ii) {
        get<alpha>(knots)[ii][0] = alphas[ii];
        get<alpha>(knots)[ii][1] = alphas[ii+knots.size()];
    }

    vu[i] = sum(j,norm(dx)<delta*region_factor,psol_u1*al[j][0] + psol_u2*al[j][1]);
    vv[i] = sum(j,norm(dx)<delta*region_factor,psol_v1*al[j][0] + psol_v2*al[j][1]);
    pr[i] = sum(j,norm(dx)<delta*region_factor,psol_p1*al[j][0] + psol_p2*al[j][1]);
    
    vtkWriteGrid("MAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
}

