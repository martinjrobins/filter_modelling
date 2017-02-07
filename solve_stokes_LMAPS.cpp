#include "solve_stokes_LMAPS.h"

void solve_stokes_LMAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0) {
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
    Accumulate<Aboria::max<double> > max;
    max.set_init(0);

    auto kernel_mq = gen_kernel_mq(i,j,c);
    auto psol_p1 = gen_psol_p1(i,j,c);
    auto psol_p2 = gen_psol_p2(i,j,c);
    auto psol_u1 = gen_psol_u1(i,j,c);
    auto psol_u2 = gen_psol_u2(i,j,c);
    auto psol_v1 = gen_psol_v1(i,j,c);
    auto psol_v2 = gen_psol_v2(i,j,c);

    //B1: u = 0 at inlet and b, p = 0 at outlet
    const double support = 1.0;
    const double cmax = eval(max(i,true,c[i]));
    std::cout << "cmax = "<<cmax<<" search rad = "<<cmax*support/c0<<std::endl;
    knots.reset_neighbour_search(1.01*cmax*support/c0);
    /*

    // generate sparse versions of particular solutions
    sparse_matrix_type kernel_matrix(knots.size(),knots.size());
    create_eigen_operator(i,j,
               kernel_mq
            , norm(dx) < (c[i]+c[j])*support/c0
            ).assemble(kernel_matrix);

    sparse_matrix_type psol_u1_matrix(knots.size(),knots.size());
    create_eigen_operator(i,j,
               psol_u1
            , norm(dx) < (c[i]+c[j])*support/c0
            ).assemble(psol_u1_matrix);

    sparse_matrix_type psol_u2_matrix(knots.size(),knots.size());
    create_eigen_operator(i,j,
               psol_u2
            , norm(dx) < (c[i]+c[j])*support/c0
            ).assemble(psol_u2_matrix);

    sparse_matrix_type psol_v1_matrix(knots.size(),knots.size());
    create_eigen_operator(i,j,
               psol_v1
            , norm(dx) < (c[i]+c[j])*support/c0
            ).assemble(psol_v1_matrix);

    sparse_matrix_type psol_v2_matrix(knots.size(),knots.size());
    create_eigen_operator(i,j,
               psol_v2
            , norm(dx) < (c[i]+c[j])*support/c0
            ).assemble(psol_v2_matrix);

    sparse_matrix_type psol_p1_matrix(knots.size(),knots.size());
    create_eigen_operator(i,j,
               psol_p1
            , norm(dx) < (c[i]+c[j])*support/c0
            ).assemble(psol_p1_matrix);

    sparse_matrix_type psol_p2_matrix(knots.size(),knots.size());
    create_eigen_operator(i,j,
               psol_p2
            , norm(dx) < (c[i]+c[j])*support/c0
            ).assemble(psol_p2_matrix);


    sparse_matrix_type alpha(psol_p1_matrix);

    typedef Eigen::Triplet<Scalar> triplet_type;
    std::vector<triplet_type> triplets_u1,triplets_u2,triplets_v1,triplets_v2,triplets_p1,triplets_p2;
    triplets_u1.reserve(support*support);
    triplets_u2.reserve(support*support);
    triplets_v1.reserve(support*support);
    triplets_v2.reserve(support*support);
    triplets_p1.reserve(support*support);
    triplets_p2.reserve(support*support);
    for (int ii = 0; ii < knots.size(); ++ii) {
        if (get<interior>(knots)[ii]) {
            auto range = box_search(knots,get<position>(knots)[ii]);
            typedef decltype(range)::iterator iterator;
            typedef iterator::reference reference;
            const size_t np = range.end()-range.begin();

            // find global indexes
            std::vector<size_t> indexs(np);
            std::transform(range.begin(),range.end(),indexs.begin(),[](reference i)
                {
                    return get<position>(i)-get<position>(knots.begin());
                });

            //
            // solve and build up weights matrix
            //
            matrix_type particular_solutions(np,np);
            vector_type kernel(np);
            vector_type alphas(np);
            // copy kernel from global to local
            for (int jj = 0; jj < np; ++jj) {
                const int gjj = indexs[jj];
                kernel(jj) = kernel_matrix(ii,jj);
            }

            auto solve_local = [&](sparse_matrix& gmatrix, triplets_type& triplets) {
                // copy psols from global to local
                for (int jj = 0; jj < np; ++jj) {
                    for (int kk = 0; kk < np; ++kk) {
                        const int gjj = indexs[jj];
                        const int gkk = indexs[kk];
                        particular_solutions(jj,kk) = gmatrix(gjj,gkk);
                    }
                }
                // and solve for alpha
                alphas = particular_solutions.colPivHouseholderQr().solve(kernel);

                // copy result to global weights matrix
                for (int jj = 0; jj < np; ++jj) {
                    triplets.push_back(triplet_type(ii,jj,alphas(jj)));
                }
            };
            solve_local(psol_u1_matrix,triplets_u1);
            solve_local(psol_u2_matrix,triplets_u2);
            solve_local(psol_v1_matrix,triplets_v1);
            solve_local(psol_v2_matrix,triplets_v2);
            solve_local(psol_p1_matrix,triplets_p1);
            solve_local(psol_p2_matrix,triplets_p2);
        }

    }


    matrix.setFromTriplets(tripletList.begin(),tripletList.end());
    */
    auto A11 = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,kernel_mq
               ,if_else(is_out[i]
                   ,psol_p1
                   ,psol_u1
                   )
               )
            , norm(dx) < c[i]*support/c0
            );

    auto A12 = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,0.0
               ,if_else(is_out[i]
                   ,psol_p2
                   ,psol_u2
                   )
               )
            , norm(dx) < c[i]*support/c0
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
            , norm(dx) < c[i]*support/c0
            );

    auto A22 = create_eigen_operator(i,j,
            if_else(is_i[i]
               ,kernel_mq
               ,if_else(is_out[i]
                   ,psol_u2
                   ,psol_v2
                   )
               )
            , norm(dx) < c[i]*support/c0
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

    vector_type source(2*knots.size());
    vector_type alphas(2*knots.size());
    sparse_matrix_type A_eigen(2*knots.size(),2*knots.size());
    matrix_type A_eigen_test(2*knots.size(),2*knots.size());
    std::cout << "assemble matrix..."<<std::endl;
    A.assemble(A_eigen);
    /*
    A.assemble(A_eigen_test);
    for (int ii=0; ii<2*knots.size(); ++ii) {
        for (int jj=0; jj<2*knots.size(); ++jj) {
            std::cout << "(ii,jj) = "<<ii<<" "<<jj<<std::endl;
            std::cout << "diff = "<<A_eigen.coeff(ii,jj) - A_eigen_test(ii,jj)<<std::endl;
            if (ii < knots.size() && jj < knots.size()) {
                std::cout << "dx = "<<(get<position>(knots)[ii]-get<position>(knots)[jj]).norm()<<std::endl;
            }
            std::cout << "A_eigen = "<<A_eigen.coeff(ii,jj)<<" A_eigen_test = "<< A_eigen_test(ii,jj)<<std::endl;
        }
    }
    */
    std::cout << "A_eigen has "<<A_eigen.nonZeros() << " non zeros"<<std::endl;

    for (int ii=0; ii<knots.size(); ++ii) {
        source(ii) = 0.0;
        if (get<inlet>(knots)[ii]) {
        //if ((get<inlet>(knots)[ii])||(get<outlet>(knots)[ii])) {
            source(knots.size()+ii) = -flow_rate;
        } else {
            source(knots.size()+ii) = 0.0;
        }

        alphas[ii] = 0.0;
        alphas[knots.size()+ii] = 0.0;
    }

    std::cout << "solve w LMAPS ..."<<std::endl;
    Eigen::SparseLU<sparse_matrix_type> solver;
    solver.compute(A_eigen);
    if(solver.info()!=Eigen::Success) {
        std::cout << "solve decomp failed  ..."<<std::endl;
        // decomposition failed
        return;
    }
    alphas = solver.solve(source);
    if(solver.info()!=Eigen::Success) {
        std::cout << "solve failed  ..."<<std::endl;
        // solving failed
        return;
    }

    double relative_error = (A_eigen*alphas - source).norm() / source.norm();
    cout << "The relative error is:\n" << relative_error << endl;
    //solve(A_eigen,alphas,source,
                //max_iter_linear,restart_linear,(linear_solver)solver_in);


    std::cout << "done solve..."<<std::endl;

    for (int ii=0; ii<knots.size(); ++ii) {
        get<alpha>(knots)[ii][0] = alphas[ii];
        get<alpha>(knots)[ii][1] = alphas[ii+knots.size()];
    }

    vu[i] = sum(j,norm(dx)<c[i]*support/c0,psol_u1*al[j][0] + psol_u2*al[j][1]);
    vv[i] = sum(j,norm(dx)<c[i]*support/c0,psol_v1*al[j][0] + psol_v2*al[j][1]);
    pr[i] = sum(j,norm(dx)<c[i]*support/c0,psol_p1*al[j][0] + psol_p2*al[j][1]);

    /*
    vu[i] = sum(j,true,psol_u1*al[j][0] + psol_u2*al[j][1]);
    vv[i] = sum(j,true,psol_v1*al[j][0] + psol_v2*al[j][1]);
    pr[i] = sum(j,true,psol_p1*al[j][0] + psol_p2*al[j][1]);
    */

    vtkWriteGrid("MAPS",0,knots.get_grid(true));

    std::cout << "done solving stokes"<<std::endl;
}
