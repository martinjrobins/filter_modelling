#include "setup_knots.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <iostream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

void setup_knots(KnotsType &knots, ParticlesType &fibres, const double fibre_radius, const double fibre_resolution_factor, const double nx, double2 domain_min, double2 domain_max, const double c0, const double k, const double epsilon_strength, const double epsilon_falloff, const bool periodic) {

    std::cout << "setup knots..." << std::endl;

    typedef typename KnotsType::position position;
    const double L = domain_max[0] - domain_min[0];
    const double delta = L/nx;
    const double fibre_resolution = fibre_resolution_factor;
    const double s = 1.50*delta;

    const double boundary_layer = delta/5;
    const double volume = ((domain_max-domain_min).prod()
                            - fibres.size()*PI*std::pow(fibre_radius,2))
                            /(domain_max-domain_min).prod();
    const int N = nx*nx*(domain_max[1]-domain_min[1])*volume/(domain_max[0]-domain_min[0]);
    double2 ns_buffer(L/2,L/2);

    if (periodic) {
        ns_buffer[0] = 0;
    }

    const int layers = 1;


    const double dt_adapt = (1.0/100.0)*PI/sqrt(2*k);
    //const double dt_adapt = 0.5*delta;

    typename KnotsType::value_type p;
    
    CDT cdt;

    // fibre boundaries
    for (int ii=0; ii<fibres.size(); ++ii) {
        Vertex_handle v1,va,vb;
        const double fdelta = delta*fibre_resolution;
        const double2 origin = get<position>(fibres)[ii];
        for (int jj=0; jj<layers; ++jj) {
            const double radius = fibre_radius-jj*fdelta;
            const double dtheta_aim = fdelta/radius;
            const int n = std::ceil(2*PI/dtheta_aim);
            const double dtheta = 2*PI/n;
            for (int kk=0; kk<n; ++kk) {
                get<position>(p) = origin + radius*double2(cos(kk*dtheta),sin(kk*dtheta));
                if (jj==0) {
                    if (kk==0) {
                        v1 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        va = v1;
                    } else {
                        vb = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        cdt.insert_constraint(va,vb);
                        va = vb;
                    } 
                }
                if (jj==0) {
                    get<target>(p) = true;
                } else {
                    get<target>(p) = false;
                }
                get<boundary>(p) = true;
                get<inlet>(p) = false;
                get<outlet>(p) = false;
                get<interior>(p) = false;
                get<kernel_constant>(p) = c0;
                knots.push_back(p);
            }
            if (jj==0) {
                cdt.insert_constraint(va,v1);
            }

        }
    }

    Vertex_handle vstart_left,vend_left,vstart_right,vend_right;
    Vertex_handle vstart_top,vend_top,vstart_bottom,vend_bottom;
    Vertex_handle va1,vb1,va2,vb2;
    if (!periodic) {
        const int ny = (domain_max[1]-domain_min[1])/(delta*fibre_resolution);
        const double deltay = (domain_max[1]-domain_min[1])/ny;
        for (int ii=1-layers; ii<ny+layers; ++ii) {
            for (int jj=0; jj<layers; ++jj) {
                // boundary - left
                get<position>(p) = double2(domain_min[0]-jj*deltay,
                        domain_min[1]+ii*deltay);
                if (jj==0) {
                    if (ii==1) {
                        vstart_left = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        va1 = vstart_left;
                    } else if (ii==ny) {
                        vend_left = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        cdt.insert_constraint(va1,vend_left);
                    } else {
                        vb1 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        cdt.insert_constraint(va1,vb1);
                        va1 = vb1;
                    }
                }


                if ((jj==0) && (ii >= 0) && (ii <= ny)) {
                    get<target>(p) = true;
                } else {
                    get<target>(p) = false;
                }
                get<boundary>(p) = true;
                get<inlet>(p) = false;
                get<interior>(p) = false;
                get<outlet>(p) = false;
                get<kernel_constant>(p) = c0;
                knots.push_back(p);

                // boundary - right
                get<position>(p) = double2(domain_max[0]+jj*deltay,
                        domain_min[1]+ii*deltay);
                if (jj==0) {
                    if (ii==1) {
                        vstart_right = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        va2 = vstart_right;
                    } else if (ii==ny) {
                        vend_right = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        cdt.insert_constraint(va1,vend_right);
                    } else {
                        vb2 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        cdt.insert_constraint(va2,vb2);
                        va2 = vb2;
                    }
                }
                if ((jj==0) && (ii >= 0) && (ii <= ny)) {
                    get<target>(p) = true;
                } else {
                    get<target>(p) = false;
                }
                get<boundary>(p) = true;
                get<inlet>(p) = false;
                get<interior>(p) = false;
                get<outlet>(p) = false;
                get<kernel_constant>(p) = c0;
                knots.push_back(p);
            }
        }
    }

    const int nxb = (domain_max[0]-domain_min[0])/(delta*fibre_resolution);
    const double deltax = (domain_max[0]-domain_min[0])/nxb;
    for (int ii=1; ii<nxb; ++ii) {
        for (int jj=0; jj<layers; ++jj) {
            // inlet
            get<position>(p) = double2(domain_min[0] + ii*deltax,domain_max[1]+jj*deltax);
            if (jj==0) {
                if (ii==1) {
                    vstart_top = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    va1 = vstart_top;
                } else if (ii==nxb) {
                    vend_top = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    cdt.insert_constraint(va1,vend_top);
                } else {
                    vb1 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    cdt.insert_constraint(va1,vb1);
                    va1 = vb1;
                }
            }
            get<target>(p) = true;
            get<boundary>(p) = false;
            get<inlet>(p) = true;
            get<outlet>(p) = false;
            get<interior>(p) = false;
            get<kernel_constant>(p) = c0;
            knots.push_back(p);
        }
    }

    for (int ii=1; ii<nxb; ++ii) {
        for (int jj=0; jj<layers; ++jj) {
            // outlet
            get<position>(p) = double2(domain_min[0] + ii*deltax,domain_min[1]-jj*deltax);
            if (jj==0) {
                if (ii==1) {
                    vstart_bottom = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    va1 = vstart_top;
                } else if (ii==nxb) {
                    vend_bottom = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    cdt.insert_constraint(va1,vend_bottom);
                } else {
                    vb1 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    cdt.insert_constraint(va1,vb1);
                    va1 = vb1;
                }
            }
            get<boundary>(p) = false;
            get<target>(p) = true;
            get<inlet>(p) = false;
            get<outlet>(p) = true;
            get<interior>(p) = false;
            get<kernel_constant>(p) = c0;
            knots.push_back(p);
        }
    }

    cdt.insert_constraint(vend_left,vstart_top);
    cdt.insert_constraint(vend_top,vend_right);
    cdt.insert_constraint(vstart_right,vend_bottom);
    cdt.insert_constraint(vstart_bottom,vstart_left);
    
    std::list<Point> list_of_seeds;
    list_of_seeds.push_back(Point(deltax, deltax));
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Meshing the domain..." << std::endl;
    CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
            Criteria());
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Number of finite faces: " << cdt.number_of_faces() << std::endl;
    int mesh_faces_counter = 0;
    while (knots.size() < N) {
        get<position>(p) = double2(uniformx(generator),uniformy(generator));
        bool outside_fibres = true;
        for (auto tpl: euclidean_search(fibres.get_query(),get<position>(p),fibre_radius+fibre_resolution*delta)) {
            outside_fibres = false;
        }
        if (outside_fibres) {
            get<boundary>(p) = false;
            get<inlet>(p) = false;
            get<outlet>(p) = false;
            get<target>(p) = true;
            get<interior>(p) = true;
            get<kernel_constant>(p) = c0;
            knots.push_back(p);
        }

    }
    std::cout << "added "<<knots.size()<<" interior knots"<<std::endl;


    for(CDT::finite_vertices_iterator v = cdt.finite_vertices_begin();
            v != cdt.finite_vertices_end(); ++v) {
        const double x = v->point().x;
        const double y = v->point().y;
        get<position>(p) = double2(x,y);
        get<boundary>(p) = false;
        get<inlet>(p) = false;
        get<outlet>(p) = false;
        get<target>(p) = true;
        get<interior>(p) = true;
        get<kernel_constant>(p) = c0;
        knots.push_back(p);
    }

    knots.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,bool2(periodic,false));
    
    std::cout << "added "<<knots.size()<<" knots with c0 = " <<c0<< std::endl;


    Symbol<boundary> is_b;
    Symbol<target> is_t;
    Symbol<inlet> is_in;
    Symbol<outlet> is_out;
    Symbol<interior> is_i;
    Symbol<position> r;
    Symbol<pressure> pr;
    Symbol<velocity_u> refd;
    Symbol<alive> alive_;
    Symbol<id> id_;
    Symbol<kernel_constant> c;

    Label<0,KnotsType> i(knots);
    Label<1,KnotsType> j(knots);
    Label<1,ParticlesType> bf(fibres);
    auto dx = create_dx(i,j);
    auto dkf = create_dx(i,bf);
    Accumulate<std::plus<double2> > sumv;
    AccumulateWithinDistance<std::plus<double2> > sumv_sparse(s);
    Accumulate<std::plus<double> > sum;
    VectorSymbolic<double,2> vector;
    Accumulate<Aboria::min<double> > min;
    min.set_init(100);

    //c[i] = min(bf,norm(dkf));
    vtkWriteGrid("init_knots",0,knots.get_grid(true));

    //auto spring_force_kk = gen_spring(i,j,k,s);
    auto spring_force_kk = deep_copy(
                if_else(dot(dx,dx)==0
                    ,vector(0,0)
                    ,(-k*(s-norm(dx))/norm(dx))*dx
                  )
            );

    auto spring_force_kf = deep_copy(
        if_else(dot(dkf,dkf)==0
            ,vector(0,0)
            ,if_else(norm(dkf)>fibre_radius+boundary_layer
                //,0.0
                ,epsilon_strength*delta*exp(1.0/epsilon_falloff*(fibre_radius+boundary_layer-norm(dkf)))
                ,-10*k*(fibre_radius+boundary_layer-norm(dkf))
                )*dkf/norm(dkf)
            )
        );


    auto spring_force_kb = deep_copy(
        vector(
                if_else(r[i][0] < domain_min[0]+boundary_layer
                ,10*k*(domain_min[0]+boundary_layer-r[i][0])
                ,if_else(r[i][0] > domain_max[0]-boundary_layer
                    ,10*k*(domain_max[0]-boundary_layer-r[i][0])
                    //,0.0
                    ,-epsilon_strength*delta*(
                            exp(1.0/epsilon_falloff*(domain_min[0]+boundary_layer-r[i][0]))
                           - exp(1.0/epsilon_falloff*(r[i][0]-domain_max[0]+boundary_layer))
                           )
                    )
                )
                ,if_else(r[i][1] < domain_min[1]+boundary_layer
                ,10*k*(domain_min[1]+boundary_layer-r[i][1])
                ,if_else(r[i][1] > domain_max[1]-boundary_layer
                    ,10*k*(domain_max[1]-boundary_layer-r[i][1])
                    //,0.0
                    ,-epsilon_strength*delta*(
                            exp(1.0/epsilon_falloff*(domain_min[1]+boundary_layer-r[i][1]))
                           - exp(1.0/epsilon_falloff*(r[i][1]-domain_max[1]+boundary_layer))
                           )
                    )
                )
            )
        );

    auto spring_force_periodic_kb = deep_copy(
        vector(
                if_else(r[i][0] < domain_min[0]+boundary_layer
                ,10*k*(domain_min[0]+boundary_layer-r[i][0])
                ,if_else(r[i][0] > domain_max[0]-boundary_layer
                    ,10*k*(domain_max[0]-boundary_layer-r[i][0])
                    ,0.0
                    //,-epsilon_strength*delta*(
                    //        exp(1.0/epsilon_falloff*(domain_min[0]+boundary_layer-r[i][0]))
                    //       - exp(1.0/epsilon_falloff*(r[i][0]-domain_max[0]+boundary_layer))
                    //       )
                    )
                )
                ,if_else(r[i][1] < domain_min[1]+boundary_layer
                ,10*k*(domain_min[1]+boundary_layer-r[i][1])
                ,if_else(r[i][1] > domain_max[1]-boundary_layer
                    ,10*k*(domain_max[1]-boundary_layer-r[i][1])
                    //,0.0
                    ,-epsilon_strength*delta*(
                            exp(1.0/epsilon_falloff*(domain_min[1]+boundary_layer-r[i][1]))
                           - exp(1.0/epsilon_falloff*(r[i][1]-domain_max[1]+boundary_layer))
                           )
                    )
                )
            )
        );


    /*
    const double mult = 2.0;
    const double scale = 7.0/(64.0*PI);
    auto kernel_wendland = deep_copy(
        scale*pow(2.0-norm(dx)/c[i],4)*(1.0+2*norm(dx)/c[i])/(c[i]*c[i])
        );

    auto kernel_grad_wendland_a = deep_copy(
        scale*(
            -4*pow(2-norm(dx)/c[i],3)*(1+2*norm(dx)/c[i])
            + 2*pow(2-norm(dx)/c[i],4)
            )/(norm(dx)*c[i]*c[i]*c[i])
            );

    auto kernel_grad_wendland_b = deep_copy(
        scale*(
            -4*pow(2-norm(dx)/c[j],3)*(1+2*norm(dx)/c[j])
            + 2*pow(2-norm(dx)/c[j],4)
            )/(norm(dx)*c[j]*c[j]*c[j])
            );



    c[i] = mult*delta;
    //const double refd = 1.0/std::pow(delta,2);

    for (int ii=0; ii<1000; ii++) {
        pr[i] = sum(j,norm(dx)<2*c[i],kernel_wendland);
        c[i] = mult*sqrt(1.0/pr[i]);

        refd[i] = min(bf,true,abs(norm(dkf)-fibre_radius));
        for(int k=0; k<knots.size(); ++k) {
            double min_dist = std::min(
                    std::abs(L-get<position>(knots)[k][0]),
                    std::abs(get<position>(knots)[k][0])
                    );
            min_dist = std::min(min_dist,get<velocity_u>(knots)[k]);
            get<velocity_u>(knots)[k] = (1.0+std::sqrt(((1-fibre_resolution)/fibre_resolution))*std::exp(-min_dist/epsilon_falloff))/std::pow(delta,2);
        }

        pr[i] = (pr[i]/refd[i] - 1.0)/pow(pr[i],2);
        //pr[i] = (pr[i]/refd[i] - 1.0);
        r[i] += if_else(is_i[i]
                ,dt_adapt*(
                    sumv(j,norm(dx)<2*c[i] && id_[i] != id_[j],
                        (pr[i]*kernel_grad_wendland_a
                        +pr[j]*kernel_grad_wendland_b)*dx)
                  //+ sumv(bf,true,spring_force_kf)
                  //+ spring_force_kb
                  )
                ,vector(0,0));
    }


    c[i] = c0*c[i]/mult;
    */

    // adapt knot locations
    for (int ii=0; ii<1000; ii++) {
        if (periodic) {
            r[i] += dt_adapt*if_else(is_i[i]
                    ,sumv_sparse(j,spring_force_kk)
                        + sumv(bf,spring_force_kf)
                        + spring_force_periodic_kb
                    ,vector(0,0)
                );
        } else {
             r[i] += dt_adapt*if_else(is_i[i]
                    ,sumv_sparse(j,spring_force_kk)
                        + sumv(bf,spring_force_kf)
                        + spring_force_kb
                    ,vector(0,0)
                );

        }

    }

    //kill anything outside domain
    alive_[i] = r[i][0] >= domain_min[0]
                    && r[i][0] <= domain_max[0]
                    && r[i][1] >= domain_min[1]
                    && r[i][1] <= domain_max[1];

    std::cout << "finshed adapting. Writing to vtk file init_knots"<<std::endl;

    vtkWriteGrid("init_knots",1,knots.get_grid(true));

    std::cout << "finshed writing to file."<<std::endl;
}

void calculate_c(KnotsType &knots, double c0, const double nx, double2 domain_min, double2 domain_max) {
    std::cout << "calculate c..."<<knots.size()<<std::endl;

    typedef typename KnotsType::position position;

    const double L = domain_max[0] - domain_min[0];
    const double delta = L/nx;
    Symbol<kernel_constant> c;

    Label<0,KnotsType> i(knots);
    Label<1,KnotsType> j(knots);
    auto dx = create_dx(i,j);
    AccumulateWithinDistance<std::plus<double> > sum;
    Accumulate<Aboria::min<double> > max;
    max.set_init(0);

    const double mult = 2.0;
    c[i] = mult*delta;
    const double scale = 7.0/(64.0*PI);
    auto kernel_wendland = deep_copy(
        scale*pow(2.0-norm(dx)/c[i],4)*(1.0+2*norm(dx)/c[i])/(c[i]*c[i])
        );

    for (int ii=0; ii<20; ++ii) {
        const double cmax = eval(max(i,c[i]));
        sum.set_max_distance(2*cmax);
        c[i] = mult*sqrt(1.0/sum(j,if_else(norm(dx)<2*c[i],kernel_wendland,0)));
    }

    c[i] = c0*c[i]/mult;
    //c[i] = c0/mult/2;
    //c[i] = delta;

    vtkWriteGrid("init_knots",2,knots.get_grid(true));
    std::cout << "done calculate c."<<std::endl;
}
