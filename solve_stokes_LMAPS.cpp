#include "solve_stokes_LMAPS.h"

void solve_stokes_LMAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0) {
    std::cout << "solving stokes..."<<knots.size()<<std::endl;

    const double flow_rate = 1.0;
    const double region_factor = 3;

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


    sparse_matrix_type solution_matrix(3*N,3*N);
    vector_type solution_source(3*N);

    typedef Eigen::Triplet<double> triplet_type;
    std::vector<triplet_type> triplets;
    triplets.reserve(10*kernel_matrix.nonZeros());
    for (int ii = 0; ii < N; ++ii) {
        if (get<boundary>(knots)[ii]) {
            triplets.push_back(triplet_type(ii,ii,1.0)); //u=0
            solution_source[ii] = 0;
            triplets.push_back(triplet_type(N+ii,N+ii,1.0)); //v=0
            solution_source[2*ii] = 0;
            //triplets.push_back(triplet_type(2*N+ii,2*N+ii,1.0)); //p=?
            solution_source[3*ii] = 0;
        } else if (get<inlet>(knots)[ii]) {
            triplets.push_back(triplet_type(ii,ii,1.0)); //u=0
            solution_source[ii] = 0;
            triplets.push_back(triplet_type(N+ii,N+ii,1.0)); //v=-flow_rate
            solution_source[2*ii] = -flow_rate;
            //triplets.push_back(triplet_type(2*N+ii,2*N+ii,1.0)); //p=?
            solution_source[3*ii] = 0;
        } else if (get<outlet>(knots)[ii]) {
            triplets.push_back(triplet_type(ii,ii,1.0)); //u=0
            solution_source[ii] = 0;
            //triplets.push_back(triplet_type(N+ii,N+ii,1.0)); //v=?
            solution_source[2*ii] = 0;
            triplets.push_back(triplet_type(2*N+ii,2*N+ii,1.0)); //p=0
            solution_source[3*ii] = 0;
        } else  {
            const size_t n = kernel_matrix.innerVector(ii).nonZeros();
            std::cout << "n = "<<n<<std::endl;
            
            vector_type kernel_row_vector(n);
            vector_type psol_u1_row_vector(n);
            vector_type psol_u2_row_vector(n);
            vector_type psol_v1_row_vector(n);
            vector_type psol_v2_row_vector(n);
            vector_type psol_p1_row_vector(n);
            vector_type psol_p2_row_vector(n);
            auto set_row_vector_f = [&](vector_type& row_vector,sparse_matrix_type& matrix) {
                int k = 0;
                for (sparse_matrix_type::InnerIterator it(kernel_matrix,ii); it; ++it,++k) {
                    row_vector(k) = it.value();
                }
            };
            set_row_vector_f(kernel_row_vector,kernel_matrix);
            set_row_vector_f(psol_u1_row_vector,psol_u1_matrix);
            set_row_vector_f(psol_u2_row_vector,psol_u2_matrix);
            set_row_vector_f(psol_v1_row_vector,psol_v1_matrix);
            set_row_vector_f(psol_v2_row_vector,psol_v2_matrix);
            set_row_vector_f(psol_p1_row_vector,psol_p1_matrix);
            set_row_vector_f(psol_p2_row_vector,psol_p2_matrix);

            matrix_type matrix(2*n,2*n);
            vector_type source(2*n);
            vector_type alphas(2*n);
            matrix.topLeftCorner(n,n) = psol_u1_row_vector*psol_u1_row_vector + 
                                          psol_v1_row_vector*psol_v1_row_vector + 
                                          psol_p1_row_vector*psol_p1_row_vector; 

            matrix.topRightCorner(n,n)= psol_u1_row_vector*psol_u2_row_vector + 
                                          psol_v1_row_vector*psol_v2_row_vector + 
                                          psol_p1_row_vector*psol_p2_row_vector; 

            matrix.bottomLeftCorner(n,n) = matrix.topRightCorner(n,n);

            matrix.bottomRightCorner(n,n)= psol_u2_row_vector*psol_u2_row_vector + 
                                          psol_v2_row_vector*psol_v2_row_vector + 
                                          psol_p2_row_vector*psol_p2_row_vector; 

            auto solver = matrix.colPivHouseholderQr();

            auto inner_solve_f = [&](decltype(psol_u1_row_vector)& vector1,
                                     decltype(psol_u1_row_vector)& vector2, int which) {
                source.head(n) = kernel_row_vector*vector1;
                source.tail(n) = kernel_row_vector*vector2;
                alphas = solver.solve(source);
                if (1) {
                    double relative_error = (matrix*alphas - source).norm() / source.norm();
                    cout << "The relative error is:\n" << relative_error << endl;
                }

                // add to solution_matrix
                int k = 0;
                for (sparse_matrix_type::InnerIterator it(kernel_matrix,ii); it; ++it,++k) {
                    const int jj = it.index();   // col index 
                    triplets.push_back(triplet_type(N+ii,which*N+jj,alphas(k)));
                    triplets.push_back(triplet_type(2*N+ii,which*N+jj,alphas(n+k)));
                }
            };

            // solve for u
            inner_solve_f(psol_u1_row_vector,psol_u2_row_vector,0);
            // solve for v
            inner_solve_f(psol_v1_row_vector,psol_v2_row_vector,1);
            // solve for p
            inner_solve_f(psol_p1_row_vector,psol_p2_row_vector,2);

        }

    }

    solution_matrix.setFromTriplets(triplets.begin(),triplets.end());
    Eigen::SparseLU<sparse_matrix_type> solver;
    solver.compute(solution_matrix);
    if (solver.info()!=Eigen::Success) {
        std::cout << "decomposition failed" << std::endl;
        return;
    }
    vector_type solution = solver.solve(solution_source);
    if (solver.info()!=Eigen::Success) {
        std::cout << "solving failed" << std::endl;
        return;
    }
    double relative_error = (solution_matrix*solution- solution_source).norm() / solution_source.norm();
    cout << "The relative error is:\n" << relative_error << endl;

    map_type(get<velocity_u>(knots).data(),N) = solution.head(N);
    map_type(get<velocity_v>(knots).data(),N) = solution.segment(N,N);
    map_type(get<pressure>(knots).data(),N) = solution.tail(N);

    vtkWriteGrid("LMAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
}
