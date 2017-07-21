#pragma once
#include "filter.h"
#include "particular_MAPS.h"

template <unsigned int Ncheb>
double solve_stokes_fMAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0) {
    std::cout << "solving stokes..."<<knots.size()<<" with c0 = "<<c0<<std::endl;

    const double flow_rate = 1.0;
    const size_t N = knots.size();

    typedef typename KnotsType::position position;
    typedef typename KnotsType::const_reference const_knot_reference;
    typedef typename position::value_type const & const_position_reference;
    typedef typename position::value_type const & const_position_reference;
    typedef typename KnotsType::const_reference const_particle_reference;

    KnotsType interior_knots;
    KnotsType inlet_boundary_knots;
    KnotsType outlet_knots;
    KnotsType side_knots;

    for (int i = 0; i < N; ++i) {
        if (get<interior>(knots[i])) {
            interior_knots.push_back(knots[i]);
            get<id>(interior_knots[interior_knots.size()-1]) = get<id>(knots)[i];
        } else if (get<outlet>(knots[i])) {
            outlet_knots.push_back(knots[i]);
            get<id>(outlet_knots[outlet_knots.size()-1]) = get<id>(knots)[i];
        } else if (get<inlet>(knots[i])) {
            inlet_boundary_knots.push_back(knots[i]);
            get<id>(inlet_boundary_knots[inlet_boundary_knots.size()-1]) = get<id>(knots)[i];
        } else if (get<boundary>(knots[i])) {
            inlet_boundary_knots.push_back(knots[i]);
            get<id>(inlet_boundary_knots[inlet_boundary_knots.size()-1]) = get<id>(knots)[i];
        } else { // b
            side_knots.push_back(knots[i]);
            get<id>(side_knots[side_knots.size()-1]) = get<id>(knots)[i];
        }
    }
    
    interior_knots.init_neighbour_search(knots.get_min(),knots.get_max(),knots.get_periodic());
    inlet_boundary_knots.init_neighbour_search(knots.get_min(),knots.get_max(),knots.get_periodic());
    outlet_knots.init_neighbour_search(knots.get_min(),knots.get_max(),knots.get_periodic());
    side_knots.init_neighbour_search(knots.get_min(),knots.get_max(),knots.get_periodic());

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_int_type;
    vector_int_type other_ids(N);
    int ii = 0;
    for (int i = 0; i < interior_knots.size(); ++i) {
        other_ids(get<id>(interior_knots)[i]) = ii; 
        ++ii;
    }
    for (int i = 0; i < inlet_boundary_knots.size(); ++i) {
        other_ids(get<id>(inlet_boundary_knots)[i]) = ii; 
        ++ii;
    }
    for (int i = 0; i < outlet_knots.size(); ++i) {
        other_ids(get<id>(outlet_knots)[i]) = ii; 
        ++ii;
    }
    for (int i = 0; i < side_knots.size(); ++i) {
        other_ids(get<id>(side_knots)[i]) = ii; 
        ++ii;
    }

    /*
    //----B1: u = 0 at inlet and b, p = 0 at outlet
    //B1: p = 0 at inlet and u=0 at b, dudy = 0 at outlet
    auto A11_interior = create_chebyshev_operator(interior_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return kernel_mq(dx,c0);
                    });
    
    auto A11_outlet = create_chebyshev_operator(outlet_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_p1(dx,c0);
                    });

    auto A11_inlet_boundary = create_chebyshev_operator(inlet_boundary_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u1(dx,c0);
                    });

    auto A11_side = create_chebyshev_operator(side_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u1(dx,c0);
                    });

    

    auto A12_interior = create_zero_operator(interior_knots,knots);

    auto A12_outlet = create_chebyshev_operator(outlet_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_p2(dx,c0);
                    });

    auto A12_inlet_boundary = create_chebyshev_operator(inlet_boundary_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u2(dx,c0);
                    
                    });

    auto A12_side = create_chebyshev_operator(side_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u2(dx,c0);
                    });

    

    //-----B2: v = -flow_rate at inlet and v=0 at b, u = 0 at outlet
    //B2: v = -flow_rate at inlet and v=0 at b, dvdy = 0 at outlet
    auto A21_interior = create_zero_operator(interior_knots,knots);

    auto A21_outlet = create_chebyshev_operator(outlet_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u1(dx,c0);
                    });

    auto A21_inlet_boundary = create_chebyshev_operator(inlet_boundary_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_v1(dx,c0);
                    });

    auto A21_side = create_chebyshev_operator(side_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_dv1dx(dx,c0);
                    });

    auto A22_interior = create_chebyshev_operator(interior_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return kernel_mq(dx,c0);
                    });
    

   auto A22_outlet = create_chebyshev_operator(outlet_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u2(dx,c0);
                    });

    auto A22_inlet_boundary = create_chebyshev_operator(inlet_boundary_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_v2(dx,c0);
                    });

    auto A22_side = create_chebyshev_operator(side_knots,knots,
            Ncheb,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_dv2dx(dx,c0);
                    });
    */


    

    //----B1: u = 0 at inlet and b, p = 0 at outlet
    //B1: p = 0 at inlet and u=0 at b, dudy = 0 at outlet
    auto A11_interior = create_h2_operator<Ncheb>(interior_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return kernel_mq(dx,c0);
                    });
   
    auto A11_outlet = create_matrix_operator(outlet_knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                        return psol_p1(dx,c0);
                    });

    auto A11_inlet_boundary = create_h2_operator<Ncheb>(inlet_boundary_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u1(dx,c0);
                    });

    auto A11_side = create_matrix_operator(side_knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                        return psol_u1(dx,c0);
                    });

    

    auto A12_interior = create_zero_operator(interior_knots,knots);

    auto A12_outlet = create_matrix_operator(outlet_knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                        return psol_p2(dx,c0);
                    });

    auto A12_inlet_boundary = create_h2_operator<Ncheb>(inlet_boundary_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_u2(dx,c0);
                    
                    });

    auto A12_side = create_matrix_operator(side_knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                        return psol_u2(dx,c0);
                    });

    

    //-----B2: v = -flow_rate at inlet and v=0 at b, u = 0 at outlet
    //B2: v = -flow_rate at inlet and v=0 at b, dvdy = 0 at outlet
    auto A21_interior = create_zero_operator(interior_knots,knots);

    auto A21_outlet = create_matrix_operator(outlet_knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                        return psol_u1(dx,c0);
                    });

    auto A21_inlet_boundary = create_h2_operator<Ncheb>(inlet_boundary_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_v1(dx,c0);
                    });

    auto A21_side = create_matrix_operator(side_knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                        return psol_dv1dx(dx,c0);
                    });

    auto A22_interior = create_h2_operator<Ncheb>(interior_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return kernel_mq(dx,c0);
                    });
    

   auto A22_outlet = create_matrix_operator(outlet_knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                        return psol_u2(dx,c0);
                    });

    auto A22_inlet_boundary = create_h2_operator<Ncheb>(inlet_boundary_knots,knots,
            [&](const_position_reference dx,
               const_position_reference a,
               const_position_reference b) {
                        return psol_v2(dx,c0);
                    });

    auto A22_side = create_matrix_operator(side_knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                        return psol_dv2dx(dx,c0);
                    });


    auto A = create_block_operator<8,2>(
            A11_interior,A12_interior,
            A11_inlet_boundary,A12_inlet_boundary,
            A11_outlet,A12_outlet,
            A11_side,A12_side,
            A21_interior,A22_interior,
            A21_inlet_boundary,A22_inlet_boundary,
            A21_outlet,A22_outlet,
            A21_side,A22_side
            );


    /*
    //----B1: u = 0 at inlet and b, p = 0 at outlet, dudx = 0 at sides
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
                    } else if (get<boundary>(a)) { // b
                        return psol_u1(dx,get<kernel_constant>(b));
                    } else {
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
                    } else if (get<boundary>(a)) { // b
                        return psol_u2(dx,get<kernel_constant>(b));
                    } else {
                        return psol_u2(dx,get<kernel_constant>(b));
                    }
                    });
    auto A21 = create_dense_operator(knots,knots,
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    if (get<interior>(a)) {
                        return 0.0;
                    } else if (get<outlet>(a)) {
                        return psol_u1(dx,get<kernel_constant>(b));
                        //return psol_dv1dy(dx,get<kernel_constant>(b));
                    } else if (get<inlet>(a) || get<boundary>(a)) { // inlet or b
                        return psol_v1(dx,get<kernel_constant>(b));
                    } else {
                        return psol_dv1dx(dx,get<kernel_constant>(b));
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
                    } else if (get<inlet>(a) || get<boundary>(a)) { // inlet or b
                        return psol_v2(dx,get<kernel_constant>(b));
                    } else {
                        return psol_dv2dx(dx,get<kernel_constant>(b));
                    }
                    });
    






    for (int i = 0; i < N; ++i) {
        const size_t indexi = other_ids(get<id>(knots)[i]);
        for (int j = 0; j < N; ++j) {
            if (std::abs(A11.coeff(i,j)-A.coeff(indexi,j)) > 0.000001) {
                std::cout << "A11 = "<<A11.coeff(i,j)<<" A = "<<A.coeff(indexi,j) << std::endl; 
            }   
        }
    }
    for (int i = 0; i < N; ++i) {
        const size_t indexi = other_ids(get<id>(knots)[i]);
        for (int j = 0; j < N; ++j) {
            if (std::abs(A12.coeff(i,j)-A.coeff(indexi,j+N)) > 0.000001) {
                std::cout << "A12 = "<<A12.coeff(i,j)<<" A = "<<A.coeff(indexi,j+N) << std::endl; 
            }   
        }
    }
    for (int i = 0; i < N; ++i) {
        const size_t indexi = other_ids(get<id>(knots)[i]);
        for (int j = 0; j < N; ++j) {
            if (std::abs(A21.coeff(i,j)-A.coeff(indexi+N,j)) > 0.000001) {
                std::cout << "A21 = "<<A21.coeff(i,j)<<" A = "<<A.coeff(indexi+N,j) << std::endl; 
            }   
        }
    }
    for (int i = 0; i < N; ++i) {
        const size_t indexi = other_ids(get<id>(knots)[i]);
        for (int j = 0; j < N; ++j) {
            if (std::abs(A22.coeff(i,j)-A.coeff(indexi+N,j+N)) > 0.000001) {
                std::cout << "A22 = "<<A22.coeff(i,j)<<" A = "<<A.coeff(indexi+N,j+N) << std::endl; 
            }   
        }
    }
    */
    
    std::cout << "A size = ("<<A.rows()<<","<<A.cols()<<")"<<std::endl;
    std::cout << "setup equations..."<<std::endl;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Map<vector_type> map_type;


    vector_type source(2*N);
    vector_type alphas(2*N);
    //matrix_type A_eigen(2*N,2*N);
    //A.assemble(A_eigen);
    //std::cout << "A_eigen size = ("<<A_eigen.rows()<<","<<A_eigen.cols()<<")"<<std::endl;


    const int inlet_start = interior_knots.size();
    const int inlet_end = inlet_start + inlet_boundary_knots.size();
    for (int ii=0; ii<N; ++ii) {
        source(ii) = 0.0;
        if (ii >= inlet_start && ii < inlet_end && get<inlet>(inlet_boundary_knots)[ii-inlet_start]) {
            source(N+ii) = -flow_rate;
        } else {
            source(N+ii) = 0.0;
        }
    }

    
    std::cout << "solve w fMAPS ..."<<std::endl;
    /*
    Eigen::MINRES<decltype(A_eigen)> gmres;
    */
    Eigen::GMRES<decltype(A),Eigen::IdentityPreconditioner> gmres;
    gmres.set_restart(2*N+1);
    gmres.setMaxIterations(2*N);
    gmres.compute(A);
    alphas = gmres.solve(source);
    //alphas = A_eigen.householderQr().solve(source);
    double relative_error = (A*alphas - source).norm() / source.norm();
    //double true_error = (A_eigen*alphas - source).norm() / source.norm();
    double true_error = 0;
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << " gmres error: "<<relative_error<<" actual error: "<<true_error<<std::endl;

    
    std::cout << "done solve..."<<std::endl;

    /*
    for (int i = 0; i < N; ++i) {
        const size_t index = other_ids(get<id>(knots)[i]);
        if (index < interior_knots.size()) {
            if (get<id>(interior_knots)[index] != get<id>(knots)[i]) {
                std::cout << "interior id = "<< get<id>(interior_knots)[index]
                          << "my id = "<<get<id>(knots)[i] << std::endl;
            }
        }
        get<alpha1>(knots)[i] = alphas(i); 
        get<alpha2>(knots)[i] = alphas(i+N); 
    }
    */
    
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
    auto psol_u1_op = create_dense_operator(knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_u1(dx,c0);
            });

      auto psol_v1_op = create_dense_operator(knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_v1(dx,c0);
            });

      auto psol_p1_op = create_dense_operator(knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_p1(dx,c0);
            });

      auto psol_u2_op = create_dense_operator(knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_u2(dx,c0);
            });

      auto psol_v2_op = create_dense_operator(knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_v2(dx,c0);
            });

      auto psol_p2_op = create_dense_operator(knots,knots,
            [&](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
            return psol_p2(dx,c0);
            });
      */


    


    map_type(get<velocity_u>(knots).data(),N) = 
        psol_u1_op*map_type(get<alpha1>(knots).data(),N) +  
        psol_u2_op*map_type(get<alpha2>(knots).data(),N);
    map_type(get<velocity_v>(knots).data(),N) = 
        psol_v1_op*map_type(get<alpha1>(knots).data(),N) +  
        psol_v2_op*map_type(get<alpha2>(knots).data(),N);
    map_type(get<pressure>(knots).data(),N) = 
        psol_p1_op*map_type(get<alpha1>(knots).data(),N) +  
        psol_p2_op*map_type(get<alpha2>(knots).data(),N);
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
