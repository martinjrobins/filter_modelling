#include "solve_stokes_MAPS.h"

void solve_stokes_MAPS(KnotsType &knots, unsigned int max_iter_linear,unsigned int restart_linear,unsigned int solver_in, double c0) {
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
    typedef Eigen::Map<vector_type> map_type;


    vector_type source(2*knots.size());
    vector_type alphas(2*knots.size());
    vector_type alphas_test(2*knots.size());
    matrix_type A_eigen(2*knots.size(),2*knots.size());
    std::cout << "assemble matrix..."<<std::endl;
    //c[i] = c0;
    A.assemble(A_eigen);

    for (int ii=0; ii<knots.size(); ++ii) {
        source(ii) = 0.0;
        //if ((get<inlet>(knots)[ii])||(get<outlet>(knots)[ii])) {
        if (get<inlet>(knots)[ii]) {
            source(knots.size()+ii) = -flow_rate;
        } else {
            source(knots.size()+ii) = 0.0;
        }

        alphas[ii] = 0.0;
        alphas[knots.size()+ii] = 0.0;
    }

    std::cout << "solve w MAPS ..."<<std::endl;
    alphas = A_eigen.colPivHouseholderQr().solve(source);
    double relative_error = (A_eigen*alphas - source).norm() / source.norm();
    cout << "The relative error is:\n" << relative_error << endl;
   // solve(A_eigen,alphas,source,
   //             max_iter_linear,restart_linear,(linear_solver)solver_in);


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
