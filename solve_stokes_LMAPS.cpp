#include "solve_stokes_LMAPS.h"

template <typename F, typename matrix_type, typename sparse_matrix_type>
void set_submatrix_f(matrix_type& submatrix,sparse_matrix_type& matrix, KnotsType& knots, size_t ii, const F& function) {
    int k = 0;
    for (typename sparse_matrix_type::InnerIterator it(matrix,ii); it; ++it,++k) {
        int j = 0;
        int kk = it.index();   // col index 
        for (typename sparse_matrix_type::InnerIterator it(matrix,ii); it; ++it,++j) {
            int jj = it.index();   // col index 
            const double mat_val = matrix.coeff(jj,kk);
            if (std::abs(mat_val) > std::numeric_limits<double>::epsilon()) {
                submatrix(k,j) = mat_val;
            } else {
                typedef KnotsType::position::value_type double_d;
                const double_d dx = get<position>(knots)[jj]-get<position>(knots)[kk];
                submatrix(k,j) = function(dx,knots[jj],knots[kk]);
            }
        }
    }
}

void solve_stokes_LMAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0) {
    std::cout << "solving stokes..."<<knots.size()<<std::endl;

    const double flow_rate = 1.0;
    const double region_factor = 2;

    typedef typename KnotsType::position position;
    typedef typename position::value_type const & const_position_reference;
    typedef typename KnotsType::const_reference const_particle_reference;
    typedef Eigen::SparseMatrix<double,Eigen::RowMajor> sparse_matrix_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Map<vector_type> map_type;

    const double support = 3.0;
    double cmake = 0;
    const double cmax = *get<kernel_constant>(std::max_element(knots.begin(),knots.end(),
        [](const_particle_reference i, const_particle_reference j) { 
            return get<kernel_constant>(i) < get<kernel_constant>(j); 
            }));
    auto support_f = [&](const double c) { return c*support/c0; };
    const double radius = support_f(cmax);

    std::cout << "cmax = "<<cmax<<" search rad = "<<radius<<std::endl;

    const size_t N = knots.size();
    // generate sparse versions of particular solutions
    sparse_matrix_type kernel_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return kernel_mq(dx,a,b);
                }
            ).assemble(kernel_matrix);
    std::cout << "nonZeros = "<<kernel_matrix.nonZeros()<<std::endl;

    sparse_matrix_type psol_u1_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_u1(dx,a,b);
                }
            ).assemble(psol_u1_matrix);

    sparse_matrix_type psol_u2_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_u2(dx,a,b);
                }
            ).assemble(psol_u2_matrix);

    sparse_matrix_type psol_v1_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_v1(dx,a,b);
                }
            ).assemble(psol_v1_matrix);

    sparse_matrix_type psol_v2_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_v2(dx,a,b);
                }
            ).assemble(psol_v2_matrix);

    sparse_matrix_type psol_p1_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_p1(dx,a,b);
                }
            ).assemble(psol_p1_matrix);

    sparse_matrix_type psol_p2_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_p2(dx,a,b);
                }
            ).assemble(psol_p2_matrix);
    sparse_matrix_type psol_dudx1_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_du1dx(dx,a,b);
                }
            ).assemble(psol_dudx1_matrix);
    sparse_matrix_type psol_dudx2_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_du2dx(dx,a,b);
                }
            ).assemble(psol_dudx2_matrix);
    sparse_matrix_type psol_dvdy1_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_dv1dy(dx,a,b);
                }
            ).assemble(psol_dvdy1_matrix);
    sparse_matrix_type psol_dvdy2_matrix(N,N);
    create_sparse_operator(knots,knots,
            [&](const_particle_reference a) {
                    return support_f(get<kernel_constant>(a));
                },
            [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_dv2dy(dx,a,b);
                }
            ).assemble(psol_dvdy2_matrix);



    sparse_matrix_type solution_matrix(2*N,2*N);
    vector_type solution_source(2*N);

    typedef Eigen::Triplet<double> triplet_type;
    std::vector<triplet_type> triplets;
    triplets.reserve(5*kernel_matrix.nonZeros());
    for (int ii = 0; ii < N; ++ii) {
        const size_t n = kernel_matrix.innerVector(ii).nonZeros();
        std::cout << "n = "<<n<<std::endl;

        vector_type kernel_row_vector(n);
        vector_type psol_dudx1_row_vector(n);
        vector_type psol_dudx2_row_vector(n);
        vector_type psol_dvdy1_row_vector(n);
        vector_type psol_dvdy2_row_vector(n);
        matrix_type psol_u1_submatrix(n,n);
        matrix_type psol_u2_submatrix(n,n);
        matrix_type psol_v1_submatrix(n,n);
        matrix_type psol_v2_submatrix(n,n);
        vector_type psol_p1_row_vector(n);
        vector_type psol_p2_row_vector(n);
        auto set_row_vector_f = [&](vector_type& row_vector,sparse_matrix_type& matrix) {
            int k = 0;
            for (sparse_matrix_type::InnerIterator it(matrix,ii); it; ++it,++k) {
                row_vector(k) = it.value();
            }
        };
        /*
        auto set_submatrix_f = [&](matrix_type& submatrix,sparse_matrix_type& matrix) {
            int k = 0;
            for (sparse_matrix_type::InnerIterator it(matrix,ii); it; ++it,++k) {
                int j = 0;
                int kk = it.index();   // col index 
                for (sparse_matrix_type::InnerIterator it(matrix,ii); it; ++it,++j) {
                    int jj = it.index();   // col index 
                    submatrix(k,j) = matrix.coeff(kk,jj);
                }
            }
        };
        */
        set_row_vector_f(kernel_row_vector,kernel_matrix);
        /*
        set_row_vector_f(psol_dudx1_row_vector,psol_dudx1_matrix);
        set_row_vector_f(psol_dudx2_row_vector,psol_dudx2_matrix);
        set_row_vector_f(psol_dvdy1_row_vector,psol_dvdy1_matrix);
        set_row_vector_f(psol_dvdy2_row_vector,psol_dvdy2_matrix);
        */
        set_submatrix_f(psol_u1_submatrix,psol_u1_matrix,knots,ii,
                [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_u1(dx,a,b);
                });
        set_submatrix_f(psol_u2_submatrix,psol_u2_matrix,knots,ii,
                [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_u2(dx,a,b);
                });
        set_submatrix_f(psol_v1_submatrix,psol_v1_matrix,knots,ii,
                [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_v1(dx,a,b);
                });
        set_submatrix_f(psol_v2_submatrix,psol_v2_matrix,knots,ii,
                [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_v2(dx,a,b);
                });
        set_row_vector_f(psol_p1_row_vector,psol_p1_matrix);
        set_row_vector_f(psol_p2_row_vector,psol_p2_matrix);

        if (0) {
            int k = 0;
            for (sparse_matrix_type::InnerIterator it(psol_u1_matrix,ii); it; ++it,++k) {
                int j = 0;
                int kk = it.index();   // col index 
                for (sparse_matrix_type::InnerIterator it(psol_u1_matrix,ii); it; ++it,++j) {
                    int jj = it.index();   // col index 
                    typedef KnotsType::position::value_type double_d;
                    const double_d dx = get<position>(knots)[jj]-get<position>(knots)[kk];
                    std::cout << psol_u1_submatrix(k,j) << " " << psol_u1(dx,knots[kk],knots[jj]) << std::endl;
                }
            }
        }

        matrix_type matrix(2*n,2*n);
        vector_type source(2*n);
        //vector_type betas(2*n);
        vector_type alphas(2*n);
        matrix.topLeftCorner(n,n) = psol_u1_submatrix; 

        matrix.topRightCorner(n,n)= psol_v1_submatrix; 

        matrix.bottomLeftCorner(n,n) = psol_u2_submatrix;

        matrix.bottomRightCorner(n,n)= psol_v2_submatrix; 

        auto solver = matrix.colPivHouseholderQr();

        auto inner_solve_f = [&](int which) {
                alphas = solver.solve(source);
                if (0) {
                    double relative_error = (matrix*alphas - source).norm() / source.norm();
                    cout << "The relative error is:\n" << relative_error << endl;
                }

                // add to solution_matrix
                int k = 0;
                for (sparse_matrix_type::InnerIterator it(kernel_matrix,ii); it; ++it,++k) {
                    const int jj = it.index();   // col index 
                    triplets.push_back(triplet_type(which*N+ii,jj,alphas(k)));
                    triplets.push_back(triplet_type(which*N+ii,N+jj,alphas(n+k)));
                }
                std::cout << std::endl;
            };
        
        if (get<interior>(knots)[ii]) {
            // solve for first alpha row
            source.head(n) = kernel_row_vector;
            source.tail(n).setZero();
            inner_solve_f(0);
            solution_source[ii] = 0;
            // solve for second alpha row
            source.head(n).setZero();
            source.tail(n) = kernel_row_vector;
            inner_solve_f(1);
            solution_source[N+ii] = 0;
        } else if (get<boundary>(knots)[ii]) {
            triplets.push_back(triplet_type(ii,ii,1.0)); //u=0
            solution_source[ii] = 0;
            triplets.push_back(triplet_type(N+ii,N+ii,1.0)); //v=0
            solution_source[N+ii] = 0;
        } else if (get<inlet>(knots)[ii]) {
            triplets.push_back(triplet_type(ii,ii,1.0)); //u=0
            solution_source[ii] = 0;
            triplets.push_back(triplet_type(N+ii,N+ii,1.0)); //v=-flow_rate
            solution_source[N+ii] = -flow_rate;
        } else if (get<outlet>(knots)[ii]) {
            triplets.push_back(triplet_type(ii,ii,1.0)); //u=0
            solution_source[ii] = 0;
            // solve for pressure
            source.head(n) = psol_p1_row_vector;
            source.tail(n) = psol_p2_row_vector;
            inner_solve_f(1);
            solution_source[N+ii] = 0; //p=0
        } else  {
            std::cout << "this doesn't look good!!!!!!"<<std::endl;
            assert(false);
        }
    }


    solution_matrix.setFromTriplets(triplets.begin(),triplets.end());
    solution_matrix.makeCompressed();
    Eigen::GMRES<sparse_matrix_type,
                    Eigen::IdentityPreconditioner> gmres;
    gmres.set_restart(2*N+1);
    gmres.setMaxIterations(2*N);
    gmres.compute(solution_matrix);
    vector_type solution = gmres.solve(solution_source);
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << " actual error: "<<(solution_matrix*solution - solution_source).norm() / solution_source.norm()<<std::endl;

    /*
    Eigen::SparseQR<sparse_matrix_type,Eigen::COLAMDOrdering<int>> solver;
    solver.compute(solution_matrix);
    if (solver.info()!=Eigen::Success) {
        if (solver.info()==Eigen::NumericalIssue) {
            std::cout << "numerical issues" << std::endl;
        } else if (solver.info()==Eigen::InvalidInput) {
            std::cout << "invalid input" << std::endl;
        }
        std::cout << "decomposition failed" << std::endl;
        return;
    }
    std::cout << "rank = "<<solver.rank()<<std::endl;
    vector_type solution = solver.solve(solution_source);
    if (solver.info()!=Eigen::Success) {
        std::cout << "solving failed" << std::endl;
        return;
    }
    double relative_error = (solution_matrix*solution- solution_source).norm() / solution_source.norm();
    cout << "The relative error is:\n" << relative_error << endl;
    */

    map_type(get<velocity_u>(knots).data(),N) = solution.head(N);
    map_type(get<velocity_v>(knots).data(),N) = solution.tail(N);

    vtkWriteGrid("LMAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
}
