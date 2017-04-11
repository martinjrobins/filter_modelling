#include "solve_stokes_LMAPS.h"

template <typename F, typename matrix_type, typename sparse_matrix_type>
void set_submatrix_f(matrix_type& submatrix,sparse_matrix_type& matrix, KnotsType& knots, size_t ii, const F& function, std::vector<std::pair<size_t,double>>& indicies) {
    int k = 0;

    for (std::pair<size_t,double>& kk: indicies) {
    //for (typename sparse_matrix_type::InnerIterator itk(matrix,ii); itk; ++itk,++k) {
        int j = 0;
        //int kk = itk.index();   // col index 
        for (std::pair<size_t,double>& jj: indicies) {
        //for (typename sparse_matrix_type::InnerIterator itj(matrix,ii); itj; ++itj,++j) {
            //int jj = itj.index();   // col index 
            const double mat_val = matrix.coeff(jj.first,kk.first);
            if (std::abs(mat_val) > std::numeric_limits<double>::epsilon()) {
                submatrix(k,j) = mat_val;
            } else {
                typedef KnotsType::position::value_type double_d;
                const double_d dx = get<position>(knots)[kk.first]-get<position>(knots)[jj.first];
                submatrix(k,j) = function(dx,knots[jj.first],knots[kk.first]);
            }
            ++j;
        }
        ++k;
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

    const double support = 6.0;
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
                    return kernel_mq(dx,get<kernel_constant>(b));
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
                    return psol_u1(dx,get<kernel_constant>(b));
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
                    return psol_u2(dx,get<kernel_constant>(b));
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
                    return psol_v1(dx,get<kernel_constant>(b));
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
                    return psol_v2(dx,get<kernel_constant>(b));
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
                    return psol_p1(dx,get<kernel_constant>(b));
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
                    return psol_p2(dx,get<kernel_constant>(b));
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
                    return psol_du1dx(dx,get<kernel_constant>(b));
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
                    return psol_du2dx(dx,get<kernel_constant>(b));
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
                    return psol_dv1dy(dx,get<kernel_constant>(b));
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
                    return psol_dv2dy(dx,get<kernel_constant>(b));
                }
            ).assemble(psol_dvdy2_matrix);



    sparse_matrix_type solution_matrix(2*N,2*N);
    sparse_matrix_type other_matrix(2*N,2*N);
    vector_type solution_source(2*N);

    typedef Eigen::Triplet<double> triplet_type;
    std::vector<triplet_type> triplets;
    std::vector<triplet_type> other_triplets;
    triplets.reserve(5*kernel_matrix.nonZeros());
    other_triplets.reserve(5*kernel_matrix.nonZeros());
    for (int ii = 0; ii < N; ++ii) {
        const size_t n = 30;

        std::cout << "n = "<<n<<" from "<<kernel_matrix.innerVector(ii).nonZeros() << std::endl;
        CHECK(kernel_matrix.innerVector(ii).nonZeros() >= n, "not enough neighbours");
        std::vector<std::pair<size_t,double>> indicies;
        indicies.reserve(3*n);
        for (sparse_matrix_type::InnerIterator it(kernel_matrix,ii); it; ++it) {
            indicies.push_back(std::make_pair(it.index(),(get<position>(knots)[it.index()]-get<position>(knots)[ii]).squaredNorm()));
        }
        std::sort(indicies.begin(),indicies.end(),[&](std::pair<size_t,double>& a, std::pair<size_t,double>& b) { return a.second < b.second; });

        
        indicies.resize(n);
        
        //const size_t n = kernel_matrix.innerVector(ii).nonZeros();

        vector_type kernel_row_vector(n);
        vector_type psol_dudx1_row_vector(n);
        vector_type psol_dudx2_row_vector(n);
        vector_type psol_dvdy1_row_vector(n);
        vector_type psol_dvdy2_row_vector(n);
        vector_type psol_u1_row_vector(n);
        vector_type psol_u2_row_vector(n);
        vector_type psol_v1_row_vector(n);
        vector_type psol_v2_row_vector(n);
        matrix_type psol_u1_submatrix(n,n);
        matrix_type psol_u2_submatrix(n,n);
        matrix_type psol_v1_submatrix(n,n);
        matrix_type psol_v2_submatrix(n,n);
        vector_type psol_p1_row_vector(n);
        vector_type psol_p2_row_vector(n);
        auto set_row_vector_f = [&](vector_type& row_vector,sparse_matrix_type& matrix) {
            int k = 0;
            double sum = 0;
            for (std::pair<size_t,double>& kk: indicies) {
                const double mat_val = matrix.coeff(ii,kk.first);
                std::cout << "ii = "<<ii<<" kk = "<<kk.first<<" value = "<<mat_val<<std::endl;
                //CHECK(std::abs(mat_val) > std::numeric_limits<double>::epsilon(),"asdf");
                row_vector(k) = mat_val;
                sum += mat_val;
                ++k;

            }
            return sum;
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
        //set_row_vector_f(psol_dudx1_row_vector,psol_dudx1_matrix);
        //set_row_vector_f(psol_dudx2_row_vector,psol_dudx2_matrix);
        /*
        set_row_vector_f(psol_dvdy1_row_vector,psol_dvdy1_matrix);
        set_row_vector_f(psol_dvdy2_row_vector,psol_dvdy2_matrix);
        */
        set_submatrix_f(psol_u1_submatrix,psol_u1_matrix,knots,ii,
                [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_u1(dx,get<kernel_constant>(b));
                },indicies);
        set_submatrix_f(psol_u2_submatrix,psol_u2_matrix,knots,ii,
                [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_u2(dx,get<kernel_constant>(b));
                },indicies);
        set_submatrix_f(psol_v1_submatrix,psol_v1_matrix,knots,ii,
                [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_v1(dx,get<kernel_constant>(b));
                },indicies);
        set_submatrix_f(psol_v2_submatrix,psol_v2_matrix,knots,ii,
                [](const_position_reference dx,
               const_particle_reference a,
               const_particle_reference b) {
                    return psol_v2(dx,get<kernel_constant>(b));
                },indicies);
        set_row_vector_f(psol_p1_row_vector,psol_p1_matrix);
        set_row_vector_f(psol_p2_row_vector,psol_p2_matrix);

        double scale = set_row_vector_f(psol_u1_row_vector,psol_u1_matrix);
        set_row_vector_f(psol_u2_row_vector,psol_u2_matrix);
        set_row_vector_f(psol_v1_row_vector,psol_v1_matrix);
        set_row_vector_f(psol_v2_row_vector,psol_v2_matrix);
        if (0) {
            int k = 0;
            for (sparse_matrix_type::InnerIterator it(psol_u1_matrix,ii); it; ++it,++k) {
                int j = 0;
                int kk = it.index();   // col index 
                for (sparse_matrix_type::InnerIterator it(psol_u1_matrix,ii); it; ++it,++j) {
                    int jj = it.index();   // col index 
                    typedef KnotsType::position::value_type double_d;
                    const double_d dx = get<position>(knots)[jj]-get<position>(knots)[kk];
                    //std::cout << psol_u1_submatrix(k,j) << " " << psol_u1(dx,knots[kk],knots[jj]) << std::endl;
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

        matrix = matrix;

        auto solver = matrix.colPivHouseholderQr();

        auto inner_solve_f = [&](int which, std::vector<triplet_type>& triplets_in) {
                alphas = solver.solve(source);
                if (1) {
                    double relative_error = (matrix*alphas - source).norm() / source.norm();
                    cout << "The relative error is:\n" << relative_error << endl;
                }

                // add to solution_matrix
                int k = 0;
                //for (sparse_matrix_type::InnerIterator it(kernel_matrix,ii); it; ++it,++k) {
                for(std::pair<size_t,double>& jj: indicies) {
                    std::cout << " solve ii = "<<ii<<" jj = "<<jj.first<<" value1 = "<<alphas(k)<<" value2 = "<<alphas(n+k)<<std::endl;
                    //const int jj = it.index();   // col index 
                    triplets_in.push_back(triplet_type(which*N+ii,jj.first,alphas(k)));
                    triplets_in.push_back(triplet_type(which*N+ii,N+jj.first,alphas(n+k)));
                    ++k;
                }
                std::cout << std::endl;
            };

        source.head(n) = psol_dudx1_row_vector;
        source.tail(n) = psol_dudx2_row_vector;
        std::cout << "solving for dudx"<<std::endl;
        inner_solve_f(0,other_triplets);
        source.head(n) = psol_p1_row_vector;
        source.tail(n) = psol_p2_row_vector;
        std::cout << "solving for pressure"<<std::endl;
        inner_solve_f(1,other_triplets);
        std::cout << "solving for v to check"<<std::endl;
        source.head(n) = psol_v1_row_vector;
        source.tail(n) = psol_v2_row_vector;
        alphas = solver.solve(source);
        int k = 0;
        for(std::pair<size_t,double>& jj: indicies) {
            std::cout << " solve ii = "<<ii<<" jj = "<<jj.first<<" value1 = "<<alphas(k)<<" value2 = "<<alphas(n+k)<<std::endl;
            ++k;
        }

        
        if (get<interior>(knots)[ii]) {
            // solve for first alpha row
            std::cout << "solving for stokes1"<<std::endl;
            source.head(n) = kernel_row_vector;
            source.tail(n).setZero();
            inner_solve_f(0,triplets);
            std::cout << "solving for stokes2"<<std::endl;
            solution_source[ii] = 0;
            // solve for second alpha row
            source.head(n).setZero();
            source.tail(n) = kernel_row_vector;
            inner_solve_f(1,triplets);
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
            inner_solve_f(1,triplets);
            solution_source[N+ii] = 0; //p=0
        } else  {
            std::cout << "this doesn't look good!!!!!!"<<std::endl;
            assert(false);
        }
    }


    solution_matrix.setFromTriplets(triplets.begin(),triplets.end());
    std::cout << "solution_matrix nonZeros = "<<solution_matrix.nonZeros()<<std::endl;
    solution_matrix.makeCompressed();
    for (sparse_matrix_type::InnerIterator it(solution_matrix,N-1); it; ++it) {
        std::cout << "ii = "<<N-1<<" jj = "<<it.index()<<" value = "<<it.value()<<std::endl;
    }
    /*
    Eigen::GMRES<sparse_matrix_type,
                    Eigen::IdentityPreconditioner> gmres;
    gmres.set_restart(2*N+1);
    gmres.setMaxIterations(2*N);
    gmres.compute(solution_matrix);
    vector_type solution(2*N);
    solution = gmres.solve(solution_source);
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << " actual error: "<<(solution_matrix*solution - solution_source).norm() / solution_source.norm()<<std::endl;
    */

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


    other_matrix.setFromTriplets(other_triplets.begin(),other_triplets.end());
    vector_type other_vector(2*N); 
    other_vector = other_matrix*solution;

    for (int i = 0; i < N; ++i) {
        get<velocity_u>(knots)[i] = solution(i);
        get<velocity_v>(knots)[i] = solution(N+i);
        get<velocity_dudx>(knots)[i] = other_vector(i);
        get<pressure>(knots)[i] = other_vector(N+i);
    }
    /*
    map_type(get<velocity_u>(knots).data(),N) = solution.head(N);
    map_type(get<velocity_v>(knots).data(),N) = solution.tail(N);
    map_type(get<velocity_dudx>(knots).data(),N) = other_vector.head(N);
    map_type(get<pressure>(knots).data(),N) = other_vector.tail(N);
    */

    vtkWriteGrid("LMAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
}
