#include "solve_stokes_MAPS.h"

void solve_stokes_MAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0) {
    std::cout << "solving stokes..."<<knots.size()<<std::endl;

    const double flow_rate = 1.0;

    typedef typename KnotsType::position position;
    typedef typename position::value_type const & const_position_reference;
    typedef typename KnotsType::const_reference const_particle_reference;

    //----B1: u = 0 at inlet and b, p = 0 at outlet
    //B1: p = 0 at inlet and u=0 at b, dudy = 0 at outlet
    auto A11 = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    if (get<interior>(a)) {
                        return kernel_mq(dx,get<kernel_constant>(b));
                    } else if (get<outlet>(a)) {
                        //return psol_du1dy(dx,a,b);
                        return psol_p1(dx,get<kernel_constant>(b));
                    } else if (get<inlet>(a)) {
                        //return psol_p1(dx,get<kernel_constant>(b));
                        return psol_u1(dx,get<kernel_constant>(b));
                    } else { // b
                        return psol_u1(dx,get<kernel_constant>(b));
                    }
                    });

    auto A12 = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    if (get<interior>(a)) {
                        return 0.0;
                    } else if (get<outlet>(a)) {
                        //return psol_du2dy(dx,get<kernel_constant>(b));
                        return psol_p2(dx,get<kernel_constant>(b));
                    } else if (get<inlet>(a)) {
                        //return psol_p2(dx,get<kernel_constant>(b));
                        return psol_u2(dx,get<kernel_constant>(b));
                    } else { // b
                        return psol_u2(dx,get<kernel_constant>(b));
                    }
                    });


    //-----B2: v = -flow_rate at inlet and v=0 at b, u = 0 at outlet
    //B2: v = -flow_rate at inlet and v=0 at b, dvdy = 0 at outlet
    auto A21 = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    if (get<interior>(a)) {
                        return 0.0;
                    } else if (get<outlet>(a)) {
                        return psol_u1(dx,get<kernel_constant>(b));
                        //return psol_dv1dy(dx,get<kernel_constant>(b));
                    } else { // inlet or b
                        return psol_v1(dx,get<kernel_constant>(b));
                    }
                    });
        
    auto A22 = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    if (get<interior>(a)) {
                        return kernel_mq(dx,get<kernel_constant>(b));
                    } else if (get<outlet>(a)) {
                        return psol_u2(dx,get<kernel_constant>(b));
                        //return psol_dv2dy(dx,get<kernel_constant>(b));
                    } else { // inlet or b
                        return psol_v2(dx,get<kernel_constant>(b));
                    }
                    });

    auto A = create_block_operator<2,2>(
            A11,A12,
            A21,A22
            );

    std::cout << "setup equations..."<<std::endl;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Map<vector_type> map_type;

    const size_t N = knots.size();

    vector_type source(2*N);
    vector_type alphas(2*N);
    //vector_type alphas_test(2*N);
    matrix_type A_eigen(2*N,2*N);
    std::cout << "assemble matrix..."<<std::endl;
    //c[i] = c0;
    A.assemble(A_eigen);

    for (int ii=0; ii<N; ++ii) {
        source(ii) = 0.0;
        if (get<inlet>(knots)[ii]) {
            source(N+ii) = -flow_rate;
        } else {
            source(N+ii) = 0.0;
        }
    }

    std::cout << "solve w MAPS ..."<<std::endl;
    alphas = A_eigen.colPivHouseholderQr().solve(source);
    double relative_error = (A_eigen*alphas - source).norm() / source.norm();
    cout << "The relative error is:\n" << relative_error << endl;
    
    //alphas.setZero(2*N);
    
   // solve(A_eigen,alphas,source,
   //             max_iter_linear,restart_linear,(linear_solver)solver_in);


    std::cout << "done solve..."<<std::endl;

    auto psol_u1_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_u1(dx,get<kernel_constant>(b));
            });

      auto psol_v1_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_v1(dx,get<kernel_constant>(b));
            });

      auto psol_p1_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_p1(dx,get<kernel_constant>(b));
            });

      auto psol_u2_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_u2(dx,get<kernel_constant>(b));
            });

      auto psol_v2_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_v2(dx,get<kernel_constant>(b));
            });

      auto psol_p2_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_p2(dx,get<kernel_constant>(b));
            });

    auto psol_du1dx_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_du1dx(dx,get<kernel_constant>(b));
            });

    auto psol_du2dx_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_du2dx(dx,get<kernel_constant>(b));
            });

    auto psol_du1dy_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_du1dy(dx,get<kernel_constant>(b));
            });

    auto psol_du2dy_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_du2dy(dx,get<kernel_constant>(b));
            });

    auto psol_dv1dx_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_dv1dx(dx,get<kernel_constant>(b));
            });

    auto psol_dv2dx_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_dv2dx(dx,get<kernel_constant>(b));
            });

    auto psol_dv1dy_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_dv1dy(dx,get<kernel_constant>(b));
            });

    auto psol_dv2dy_op = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_dv2dy(dx,get<kernel_constant>(b));
            });


    map_type(get<alpha1>(knots).data(),N) = alphas.head(N);
    map_type(get<alpha2>(knots).data(),N) = alphas.tail(N);
    map_type(get<velocity_u>(knots).data(),N) = 
        psol_u1_op*alphas.head(N) + 
        psol_u2_op*alphas.tail(N);
    map_type(get<velocity_v>(knots).data(),N) = 
        psol_v1_op*alphas.head(N) + 
        psol_v2_op*alphas.tail(N);
    map_type(get<pressure>(knots).data(),N) = 
        psol_p1_op*alphas.head(N) + 
        psol_p2_op*alphas.tail(N);
    map_type(get<velocity_dudx>(knots).data(),N) = 
        psol_du1dx_op*alphas.head(N) + 
        psol_du2dx_op*alphas.tail(N);
    map_type(get<velocity_dudy>(knots).data(),N) = 
        psol_du1dy_op*alphas.head(N) + 
        psol_du2dy_op*alphas.tail(N);
    map_type(get<velocity_dvdx>(knots).data(),N) = 
        psol_dv1dx_op*alphas.head(N) + 
        psol_dv2dx_op*alphas.tail(N);
    map_type(get<velocity_dvdy>(knots).data(),N) = 
        psol_dv1dy_op*alphas.head(N) + 
        psol_dv2dy_op*alphas.tail(N);



    vtkWriteGrid("MAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
}
