#include "solve_laplace_BEM.h"
#include <algorithm>


double solve_laplace_BEM(KnotsType &knots, ElementsType& elements, const double fibre_charge) {
    std::cout << "solving stokes with BEM..."<<elements.size()<<std::endl;

    const double pot = fibre_charge;
    const double mu = 1.0;

    typedef typename KnotsType::position position;
    typedef typename position::value_type const & const_position_reference;
    typedef typename KnotsType::const_reference const_particle_reference;

    const vdouble2 min = elements.get_min();
    const vdouble2 max = elements.get_max();
    const double h = (get<point_b>(elements)[0] - get<point_a>(elements)[0]).norm();


    auto Aslp = create_dense_operator(elements,elements,
            make_laplace_SLP_2d1p(min,max,true));

    std::cout << "setup equations..."<<std::endl;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Map<vector_type> map_type;

    const size_t N = elements.size();

    vector_type potential(N);
    vector_type source(N);
    vector_type alphas(N);
    matrix_type A_eigen(N,N);

    std::cout << "assemble matrix..."<<std::endl;
    Aslp.assemble(A_eigen);

    potential = vector_type::Ones(N)*pot;
    source = -2*PI*potential;

    std::cout << "solve w BEM ..."<<std::endl;
    alphas = A_eigen.householderQr().solve(source);
    //alphas = A_eigen.colPivHouseholderQr().solve(source);
    double relative_error = (A_eigen*alphas - source).norm() / source.norm();
    std::cout << "The relative error is:\n" << relative_error << std::endl;
    
    std::cout << "done solve..."<<std::endl;

    for (int ii=0; ii<N; ++ii) {
        get<gradP>(elements)[ii] = alphas[ii];
    }

    vtkWriteGrid("BEMelements",0,elements.get_grid(true));

    std::cout << "assemble knot matrix..."<<std::endl;
    auto AknotsSLP = create_dense_operator(knots,elements,
            make_laplace_SLP_2d1p(min,max,false));
    auto AknotsDLP = create_dense_operator(knots,elements,
            make_laplace_gradSLP_2d1p(min,max,false));

    AknotsSLP.get_first_kernel().evaluate(get<knotpot>(knots),get<gradP>(elements));
    AknotsDLP.get_first_kernel().evaluate(get<gradknotpot>(knots),get<gradP>(elements));
    std::transform(std::begin(get<gradknotpot>(knots)),
                   std::end(get<gradknotpot>(knots)),
                    std::begin(get<gradknotpot>(knots)),
                   [](auto& i) { return i/(-2*PI); });
    std::transform(std::begin(get<knotpot>(knots)),
                   std::end(get<knotpot>(knots)),
                   std::begin(get<knotpot>(knots)),
                   [](auto& i) { return i/(-2*PI); });
    
    vtkWriteGrid("BEMknots",0,knots.get_grid(true));

    std::cout << "done solving laplace"<<std::endl;
    return relative_error;
}
