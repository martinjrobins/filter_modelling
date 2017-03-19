#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_LIST_SIZE 30

#include "Aboria.h"

using namespace Aboria;
ABORIA_VARIABLE(kernel_constant, double, "kernel_constant")

#include "particular_MAPS.h"

int main() {

    const size_t N = 100;
    ABORIA_VARIABLE(u1, double, "u1");
    ABORIA_VARIABLE(u2, double, "u2");
    ABORIA_VARIABLE(v1, double, "v1");
    ABORIA_VARIABLE(v2, double, "v2");
    ABORIA_VARIABLE(p1, double, "p1");
    ABORIA_VARIABLE(p2, double, "p2");
    ABORIA_VARIABLE(dudx1, double, "dudx1");
    ABORIA_VARIABLE(dudx2, double, "dudx2");
    ABORIA_VARIABLE(dudy1, double, "dudy1");
    ABORIA_VARIABLE(dudy2, double, "dudy2");
    ABORIA_VARIABLE(dvdx1, double, "dvdx1");
    ABORIA_VARIABLE(dvdx2, double, "dvdx2");
    /*
    ABORIA_VARIABLE(dvdy1, double, "dvdy1");
    ABORIA_VARIABLE(dvdy2, double, "dvdy2");
    */
    ABORIA_VARIABLE(stokes1, double, "stokes1");
    ABORIA_VARIABLE(stokes2, double, "stokes2");
    ABORIA_VARIABLE(mq, double, "mq");
    typedef Particles<std::tuple<u1,u2,v1,v2,p1,p2,dudx1,dudx2,dudy1,dudy2,dvdx1,dvdx2,mq,stokes1,stokes2>, 2> KnotsType; 
    typedef Particles<std::tuple<kernel_constant>, 2> BasisType; 
    KnotsType knots(N*N);
    BasisType basis(1);
    typedef KnotsType::reference reference;
    typedef KnotsType::const_reference const_reference;
    typedef KnotsType::position position;

    const double dx = 1.0/N;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const size_t index = i*N + j;
            get<position>(knots)[index] = double2(i*dx, j*dx);
            if ((i==N/2) && (j==N/2)) {
                get<position>(basis)[0] = double2(i*dx,j*dx);
                get<kernel_constant>(basis)[0] = 1.0;
            }

        }
    }

    for (reference knot: knots) {
        const double2 dx = get<position>(basis)[0] - get<position>(knot);
        get<u1>(knot) = psol_u1(dx,knot,basis[0]);
        get<u2>(knot) = psol_u2(dx,knot,basis[0]);
        get<v1>(knot) = psol_v1(dx,knot,basis[0]);
        get<v2>(knot) = psol_v2(dx,knot,basis[0]);
        get<p1>(knot) = psol_p1(dx,knot,basis[0]);
        get<p2>(knot) = psol_p2(dx,knot,basis[0]);
        get<dudx1>(knot) = psol_du1dx(dx,knot,basis[0]);
        get<dudx2>(knot) = psol_du2dx(dx,knot,basis[0]);
        get<dudy1>(knot) = psol_du1dy(dx,knot,basis[0]);
        get<dudy2>(knot) = psol_du2dy(dx,knot,basis[0]);
        get<dvdx1>(knot) = psol_dv1dx(dx,knot,basis[0]);
        get<dvdx2>(knot) = psol_dv2dx(dx,knot,basis[0]);
        /*
        get<dvdy1>(knot) = psol_dv1dy(dx,knot,basis[0]);
        get<dvdy2>(knot) = psol_dv2dy(dx,knot,basis[0]);
        */
        get<mq>(knot) = kernel_mq(dx,knot,basis[0]);

            }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const size_t index = i*N + j;
            // check that dpdx = dudxx + dudyy and 
            //            dpdy = dvdxx + dvdyy with FD
            if ((i >= 1) && (j >= 1) && (i < N-1) && (j < N-1)) {
                const size_t up = i*N+j+1;
                const size_t down = i*N+j-1;
                const size_t left = (i-1)*N+j;
                const size_t right = (i+1)*N+j;
                const double dpdx = (get<p1>(knots)[right]+get<p2>(knots)[right] - 
                                    get<p1>(knots)[left]-get<p2>(knots)[left])/(2*dx);
                const double dpdy = (get<p1>(knots)[up]+get<p2>(knots)[up] - 
                                    get<p1>(knots)[down]-get<p2>(knots)[down])/(2*dx);
                const double dudxx =    ((get<u1>(knots)[left]+get<u2>(knots)[left])
                                    - 2*(get<u1>(knots)[index]+get<u2>(knots)[index])
                                    +   (get<u1>(knots)[right]+get<u2>(knots)[right]))
                                            /(dx*dx);
                const double dudyy =    ((get<u1>(knots)[up]+get<u2>(knots)[up])
                                    - 2*(get<u1>(knots)[index]+get<u2>(knots)[index])
                                    +   (get<u1>(knots)[down]+get<u2>(knots)[down]))
                                            /(dx*dx);
                const double dvdxx =    ((get<v1>(knots)[left]+get<v2>(knots)[left])
                                    - 2*(get<v1>(knots)[index]+get<v2>(knots)[index])
                                    +   (get<v1>(knots)[right]+get<v2>(knots)[right]))
                                            /(dx*dx);
                const double dvdyy =    ((get<v1>(knots)[up]+get<v2>(knots)[up])
                                    - 2*(get<v1>(knots)[index]+get<v2>(knots)[index])
                                    +   (get<v1>(knots)[down]+get<v2>(knots)[down]))
                                            /(dx*dx);
                get<stokes1>(knots)[index] = -dpdx + (dudxx+dudyy) - get<mq>(knots)[index];
                get<stokes2>(knots)[index] = -dpdy + (dvdxx+dvdyy) - get<mq>(knots)[index];
            }
        }
    }



    vtkWriteGrid("check_solutions",0,knots.get_grid(true));

}
