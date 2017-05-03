#include "solve_stokes_Compact.h"

void solve_stokes_Compact(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0) {
    std::cout << "solving stokes..."<<knots.size()<<std::endl;

    const double flow_rate = 1.0;

    typedef typename KnotsType::position position;
    typedef typename position::value_type const & const_position_reference;
    typedef typename KnotsType::const_reference const_reference;

    const double h = 15.0;
    const double invh = 1.0/h;
    const double mu = 1.0;

    const double search_radius = h;

    auto A11 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(i) || get<inlet>(i) || get<outlet>(i)) {
                        if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                            return -kernel_yy(dx,invh);
                        } else {
                            return mu*laplace_yy(dx,invh);
                        }
                    } else {
                        if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                            return mu*laplace_yy(dx,invh);
                        } else {
                            return -mu*mu*laplace2_yy(dx,invh) - kernel_xx(dx,invh);
                        }
                    }
               });

        auto A12 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(i) || get<inlet>(i) || get<outlet>(i)) {
                        if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                            return kernel_xy(dx,invh);
                        } else {
                            return -mu*laplace_xy(dx,invh);
                        }
                    } else {
                        if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                            return -mu*laplace_xy(dx,invh);
                        } else {
                            return mu*mu*laplace2_xy(dx,invh) - kernel_xy(dx,invh);
                        }
                    }
               });

        auto A22 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(i) || get<inlet>(i) || get<outlet>(i)) {
                        if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                            return -kernel_xx(dx,invh);
                        } else {
                            return mu*laplace_xx(dx,invh);
                        }
                    } else {
                        if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                            return mu*laplace_xx(dx,invh);
                        } else {
                            return -mu*mu*laplace2_xx(dx,invh) - kernel_yy(dx,invh);
                        }
                    }
               });


        auto A = create_block_operator<2,2>(A11, A12,
                                            A12, A22);

        auto B11 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                        return -kernel_yy(dx,invh);
                    } else {
                        return mu*laplace_yy(dx,invh);
                    }
               });

        auto B12 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                        return kernel_xy(dx,invh);
                    } else {
                        return -mu*laplace_xy(dx,invh);
                    }
               });

        auto B22 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                        return -kernel_xx(dx,invh);
                    } else {
                        return mu*laplace_xx(dx,invh);
                    }
               });

        auto B31 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                        return 0.0;
                    } else {
                        return -kernel_x(dx,invh);
                    }
               });

        auto B32 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j) || get<inlet>(j) || get<outlet>(j)) {
                        return 0.0;
                    } else {
                        return -kernel_y(dx,invh);
                    }
               });

        
    /*
    auto A11 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(i) || get<inlet>(i)) {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return -kernel_yy(dx, invh);
                    } else if (get<outlet>(j)) {
                        return -kernel_yy(dx, invh);
                    } else {
                        return mu*laplace_yy(dx, invh);
                    }
                } else if (get<outlet>(i)) {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return -kernel_yy(dx, invh);
                    } else if (get<outlet>(j)) {
                        return -kernel_yy(dx, invh);
                    } else {
                        return mu*laplace_yy(dx, invh);
                    }
                } else {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return mu*laplace_yy(dx, invh);
                    } else if (get<outlet>(j)) {
                        return mu*laplace_yy(dx, invh);
                    } else {
                        return -mu*mu*laplace2_yy(dx, invh) - kernel_xx(dx, invh);
                    }
                }
           });
    auto A12 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(i) || get<inlet>(i)) {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return kernel_xy(dx, invh);
                    } else if (get<outlet>(j)) {
                        return 0.0;
                    } else {
                        return -mu*laplace_xy(dx, invh);
                    }
                } else if (get<outlet>(i)) {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return kernel_xy(dx, invh);
                    } else if (get<outlet>(j)) {
                        return 0.0;
                    } else {
                        return -mu*laplace_xy(dx, invh);
                    }
                } else {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return -mu*laplace_xy(dx, invh);
                    } else if (get<outlet>(j)) {
                        return kernel_x(dx, invh);
                    } else {
                        return mu*mu*laplace2_xy(dx, invh) - kernel_xy(dx, invh);
                    }
                }
           });

    auto A21 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(i) || get<inlet>(i)) {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return kernel_xy(dx, invh);
                    } else if (get<outlet>(j)) {
                        return kernel_xy(dx, invh);
                    } else {
                        return -mu*laplace_xy(dx, invh);
                    }
                   
                } else if (get<outlet>(i)) {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return 0.0;
                    } else if (get<outlet>(j)) {
                        return 0.0;
                    } else {
                        return -kernel_x(dx, invh);
                    }
                } else {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return -mu*laplace_xy(dx, invh);
                    } else if (get<outlet>(j)) {
                        return -mu*laplace_xy(dx, invh);
                    } else {
                        return mu*mu*laplace2_xy(dx, invh) - kernel_xy(dx, invh);
                    }
                }
           });

    auto A22 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(i) || get<inlet>(i)) {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return -kernel_xx(dx, invh);
                    } else if (get<outlet>(j)) {
                        return 0.0;
                    } else {
                        return mu*laplace_xx(dx, invh);
                    }
                } else if (get<outlet>(i)) {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return 0.0;
                    } else if (get<outlet>(j)) {
                        return kernel(dx, invh);
                    } else {
                        return -kernel_y(dx, invh);
                    }
                    
                } else {
                    if (get<boundary>(j) || get<inlet>(j)) {
                        return mu*laplace_xx(dx, invh);
                    } else if (get<outlet>(j)) {
                        return kernel_y(dx, invh);
                    } else {
                        return -mu*mu*laplace2_xx(dx, invh) - kernel_yy(dx, invh);
                    }
                }
           });

    
    auto A = create_block_operator<2,2>(A11, A12,
                                        A21, A22);


    auto B11 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(j) || get<inlet>(j)) {
                    return -kernel_yy(dx, invh);
                } else if (get<outlet>(j)) {
                    return -kernel_yy(dx, invh);
                } else {
                    return mu*laplace_yy(dx, invh);
                }
           });

    auto B12 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(j) || get<inlet>(j)) {
                    return kernel_xy(dx, invh);
                } else if (get<outlet>(j)) {
                    return 0.0;
                } else {
                    return -mu*laplace_xy(dx, invh);
                }
           });

    auto B21 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(j) || get<inlet>(j)) {
                    return kernel_xy(dx, invh);
                } else if (get<outlet>(j)) {
                    return kernel_xy(dx, invh);
                } else {
                    return -mu*laplace_xy(dx, invh);
                }
           });

    auto B22 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(j) || get<inlet>(j)) {
                    return -kernel_xx(dx, invh);
                } else if (get<outlet>(j)) {
                    return 0.0;
                } else {
                    return mu*laplace_xx(dx, invh);
                }
           });

    auto B31 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(j) || get<inlet>(j)) {
                    return 0.0;
                } else if (get<outlet>(j)) {
                    return 0.0;
                } else {
                    return -kernel_x(dx, invh);
                }
           });

    auto B32 = create_sparse_operator(knots,knots,search_radius,
            [&](const double2& dx,
                const_reference i,
                const_reference j) {
                if (get<boundary>(j) || get<inlet>(j)) {
                    return 0.0;
                } else if (get<outlet>(j)) {
                    return kernel(dx, invh);
                } else {
                    return -kernel_y(dx, invh);
                }
           });
    */


    std::cout << "setup equations..."<<std::endl;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Map<vector_type> map_type;

    const size_t N = knots.size();

    vector_type source(2*N);
    vector_type alphas(2*N);
    matrix_type A_eigen(2*N,2*N);
    std::cout << "assemble matrix..."<<std::endl;
    A.assemble(A_eigen);

    for (int ii=0; ii<N; ++ii) {
        source(ii) = 0.0;
        if (get<inlet>(knots)[ii]||get<outlet>(knots)[ii]) {
            source(N+ii) = -flow_rate;
        } else {
            source(N+ii) = 0.0;
        }
    }

    std::cout << "solve w Compact ..."<<std::endl;
    alphas = A_eigen.colPivHouseholderQr().solve(source);
    double relative_error = (A_eigen*alphas - source).norm() / source.norm();
    cout << "The relative error is:\n" << relative_error << endl;
    
    //alphas.setZero(2*N);
    
   // solve(A_eigen,alphas,source,
   //             max_iter_linear,restart_linear,(linear_solver)solver_in);


    std::cout << "done solve..."<<std::endl;

    map_type u(get<velocity_u>(knots).data(),N);
    map_type v(get<velocity_v>(knots).data(),N);
    map_type pr(get<pressure>(knots).data(),N);

    u = B11*alphas.head(N) + B12*alphas.tail(N);
    v = B12*alphas.head(N) + B22*alphas.tail(N);
    pr = B31*alphas.head(N) + B32*alphas.tail(N);

    vtkWriteGrid("MAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
}
