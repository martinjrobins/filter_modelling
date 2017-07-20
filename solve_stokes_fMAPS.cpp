#include "solve_stokes_fMAPS.h"
#ifdef HAVE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif


template <unsigned int Ncheb>
void solve_stokes_fMAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0) {
    std::cout << "solving stokes..."<<knots.size()<<" with c0 = "<<c0<<std::endl;

    const double flow_rate = 1.0;
    const size_t N = knots.size();

    typedef typename KnotsType::position position;
    typedef typename KnotsType::const_reference const_knot_reference;
    typedef typename position::value_type const & const_position_reference;
    typedef typename position::value_type const & const_position_reference;

    KnotsType interior_knots;
    KnotsType inlet_boundary_knots;
    KnotsType outlet_knots;
    KnotsType side_knots;

    for (int i = 0; i < N; ++i) {
        if (get<interior>(knots[i])) {
            interior_knots.push_back(knots[i]);
        } else if (get<outlet>(knots[i])) {
            outlet_knots.push_back(knots[i]);
        } else if (get<inlet>(knots[i])) {
            inlet_boundary_knots.push_back(knots[i]);
        } else if (get<boundary>(knots[i])) {
            inlet_boundary_knots.push_back(knots[i]);
        } else { // b
            side_knots.push_back(knots[i]);
        }
    }

    int ii = 0;
    for (int i = 0; i < interior_knots.size(); ++i) {
        knots[ii] = interior_knots[i];
        ++ii;
    }
    for (int i = 0; i < inlet_boundary_knots.size(); ++i) {
        knots[ii] = inlet_boundary_knots[i];
        ++ii;
    }
    for (int i = 0; i < outlet_knots.size(); ++i) {
        knots[ii] = outlet_knots[i];
        ++ii;
    }
    for (int i = 0; i < side_knots.size(); ++i) {
        knots[ii] = side_knots[i];
        ++ii;
    }
    
    interior_knots.init_neighbour_search(knots.get_min(),knots.get_max(),knots.get_periodic());
    inlet_boundary_knots.init_neighbour_search(knots.get_min(),knots.get_max(),knots.get_periodic());
    outlet_knots.init_neighbour_search(knots.get_min(),knots.get_max(),knots.get_periodic());
    side_knots.init_neighbour_search(knots.get_min(),knots.get_max(),knots.get_periodic());

    //----B1: u = 0 at inlet and b, p = 0 at outlet
    //B1: p = 0 at inlet and u=0 at b, dudy = 0 at outlet
    auto A11_interior = create_h2_operator<Ncheb>(interior_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return kernel_mq(dx,c0);
                    });
    /*
    auto A11_interior = create_dense_operator(interior_knots,knots,
            [&](const_position_reference dx,
               const_knot_reference a,
               const_knot_reference b) {
                        return kernel_mq(dx,c0);
                    });
                    */

    auto A11_outlet = create_h2_operator<Ncheb>(outlet_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_p1(dx,c0);
                    });

    auto A11_inlet_boundary = create_h2_operator<Ncheb>(inlet_boundary_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u1(dx,c0);
                    });

    auto A11_side = create_h2_operator<Ncheb>(side_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u1(dx,c0);
                    });


    auto A12_interior = create_zero_operator(interior_knots,knots);

    auto A12_outlet = create_h2_operator<Ncheb>(outlet_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_p2(dx,c0);
                    });

    auto A12_inlet_boundary = create_h2_operator<Ncheb>(inlet_boundary_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u2(dx,c0);
                    
                    });

    auto A12_side = create_h2_operator<Ncheb>(side_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u2(dx,c0);
                    });


    //-----B2: v = -flow_rate at inlet and v=0 at b, u = 0 at outlet
    //B2: v = -flow_rate at inlet and v=0 at b, dvdy = 0 at outlet
    auto A21_interior = create_zero_operator(interior_knots,knots);

    auto A21_outlet = create_h2_operator<Ncheb>(outlet_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u1(dx,c0);
                    });

    auto A21_inlet_boundary = create_h2_operator<Ncheb>(inlet_boundary_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_v1(dx,c0);
                    });

    auto A21_side = create_h2_operator<Ncheb>(side_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_dv1dx(dx,c0);
                    });

    auto A22_interior = create_h2_operator<Ncheb>(interior_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return kernel_mq(dx,c0);
                    });
    /*
    auto A22_interior = create_dense_operator(interior_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return kernel_mq(dx,c0);
                    });
                    */


   auto A22_outlet = create_h2_operator<Ncheb>(outlet_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u2(dx,c0);
                    });

    auto A22_inlet_boundary = create_h2_operator<Ncheb>(inlet_boundary_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_v2(dx,c0);
                    });

    auto A22_side = create_h2_operator<Ncheb>(side_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_dv2dx(dx,c0);
                    });


    auto A = create_block_operator<8,2>(
            A11_interior,A12_interior,
            A21_interior,A22_interior,
            A11_inlet_boundary,A12_inlet_boundary,
            A21_inlet_boundary,A22_inlet_boundary,
            A11_outlet,A12_outlet,
            A21_outlet,A22_outlet,
            A11_side,A12_side,
            A21_side,A22_side
            );

    /*
    ii = 0;
    for (int i = 0; i < interior_knots.size(); ++i) {
        std::cout << "id check "<<get<position>(interior_knots)[i]<<" "<<get<position>(knots)[ii]<<std::endl;
        ++ii;
    }
    */


    std::cout << "A size = ("<<A.rows()<<","<<A.cols()<<")"<<std::endl;
    std::cout << "setup equations..."<<std::endl;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Map<vector_type> map_type;


    vector_type source(2*N);
    vector_type alphas(2*N);
    matrix_type A_eigen(2*N,2*N);
    A.assemble(A_eigen);
    //std::cout << "A_eigen size = ("<<A_eigen.rows()<<","<<A_eigen.cols()<<")"<<std::endl;


    for (int ii=0; ii<N; ++ii) {
        source(ii) = 0.0;
        if (get<inlet>(knots)[ii]) {
            source(N+ii) = -flow_rate;
        } else {
            source(N+ii) = 0.0;
        }
    }

    
    std::cout << "solve w fMAPS ..."<<std::endl;
    Eigen::GMRES<decltype(A),
                    Eigen::IdentityPreconditioner> gmres;
    gmres.set_restart(2*N+1);
    gmres.setMaxIterations(2*N);
    gmres.compute(A);
    alphas = gmres.solve(source);
    //alphas = A_eigen.colPivHouseholderQr().solve(source);
    double relative_error = (A_eigen*alphas - source).norm() / source.norm();
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << " gmres error: "<<(A*alphas - source).norm() / source.norm()<<" actual error: "<<(A_eigen*alphas - source).norm() / source.norm()<<std::endl;

    
    std::cout << "done solve..."<<std::endl;

    map_type(get<alpha1>(knots).data(),N) = alphas.head(N);
    map_type(get<alpha2>(knots).data(),N) = alphas.tail(N);

    auto psol_u1_op = create_h2_operator<Ncheb>(knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_u1(dx,c0);
            });

      auto psol_v1_op = create_h2_operator<Ncheb>(knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_v1(dx,c0);
            });

      auto psol_p1_op = create_h2_operator<Ncheb>(knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_p1(dx,c0);
            });

      auto psol_u2_op = create_h2_operator<Ncheb>(knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_u2(dx,c0);
            });

      auto psol_v2_op = create_h2_operator<Ncheb>(knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_v2(dx,c0);
            });

      auto psol_p2_op = create_h2_operator<Ncheb>(knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_p2(dx,c0);
            });

      /*
    auto psol_du1dx_op = create_h2_operator<Ncheb>(knots,knots,
            n,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_du1dx(dx,c0);
            });

    auto psol_du2dx_op = create_h2_operator<Ncheb>(knots,knots,
            n,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_du2dx(dx,c0);
            });

    auto psol_du1dy_op = create_h2_operator<Ncheb>(knots,knots,
            n,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_du1dy(dx,c0);
            });

    auto psol_du2dy_op = create_h2_operator<Ncheb>(knots,knots,
            n,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_du2dy(dx,c0);
            });

    auto psol_dv1dx_op = create_h2_operator<Ncheb>(knots,knots,
            n,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_dv1dx(dx,c0);
            });

    auto psol_dv2dx_op = create_h2_operator<Ncheb>(knots,knots,
            n,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_dv2dx(dx,c0);
            });

    auto psol_dv1dy_op = create_h2_operator<Ncheb>(knots,knots,
            n,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_dv1dy(dx,c0);
            });

    auto psol_dv2dy_op = create_h2_operator<Ncheb>(knots,knots,
            n,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
            return psol_dv2dy(dx,c0);
            });
            */


    map_type(get<velocity_u>(knots).data(),N) = 
        psol_u1_op*alphas.head(N) + 
        psol_u2_op*alphas.tail(N);
    map_type(get<velocity_v>(knots).data(),N) = 
        psol_v1_op*alphas.head(N) + 
        psol_v2_op*alphas.tail(N);
    map_type(get<pressure>(knots).data(),N) = 
        psol_p1_op*alphas.head(N) + 
        psol_p2_op*alphas.tail(N);
    /*
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
        */

    vtkWriteGrid("fMAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
    return relative_error;
}
