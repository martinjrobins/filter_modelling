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
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

template <class CDT>
class My_delaunay_mesh_size_criteria_2 : 
    public virtual CGAL::Delaunay_mesh_size_criteria_2<CDT>
{
protected:
  typedef typename CDT::Geom_traits Geom_traits;

public:
  typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Base;

  My_delaunay_mesh_size_criteria_2(const double aspect_bound = 0.125, 
                                const double size_bound = 0,
                                const Geom_traits& traits = Geom_traits())
    : Base(aspect_bound,size_bound,traits) {}

    typedef typename Base::Quality Quality; 

    class Is_bad: public Base::Is_bad
    {
    protected:
        const double squared_size_bound; // squared size bound on edge length
    public:
        typedef typename Base::Is_bad::Point_2 Point_2;

        Is_bad(const double aspect_bound,
                const double size_bound,
                const Geom_traits& traits)
            : Base::Is_bad(aspect_bound, traits),
            squared_size_bound(size_bound * size_bound) {}

        Mesh_2::Face_badness operator()(const Quality q) const
        {
            if( q.size() > 1 )
                return Mesh_2::IMPERATIVELY_BAD;
            if( q.sine() < this->B )
                return Mesh_2::BAD;
            else
                return Mesh_2::NOT_BAD;
        }

        Mesh_2::Face_badness operator()(const typename CDT::Face_handle& fh,
                Quality& q) const
        {
            typedef typename CDT::Geom_traits Geom_traits;
            typedef typename Geom_traits::Compute_area_2 Compute_area_2;
            typedef typename Geom_traits::Compute_squared_distance_2
                Compute_squared_distance_2;

            Geom_traits traits; /** @warning traits with data!! */

            Compute_squared_distance_2 squared_distance = 
                traits.compute_squared_distance_2_object();

            const Point_2& pa = fh->vertex(0)->point();
            const Point_2& pb = fh->vertex(1)->point();
            const Point_2& pc = fh->vertex(2)->point();

            double
                a = CGAL::to_double(squared_distance(pb, pc)),
                  b = CGAL::to_double(squared_distance(pc, pa)),
                  c = CGAL::to_double(squared_distance(pa, pb));

            double max_sq_length; // squared max edge length
            double second_max_sq_length;

            if(a<b)
            {
                if(b<c) {
                    max_sq_length = c;
                    second_max_sq_length = b;
                }
                else { // c<=b
                    max_sq_length = b;
                    second_max_sq_length = ( a < c ? c : a );
                }
            }
            else // b<=a
            {
                if(a<c) {
                    max_sq_length = c;
                    second_max_sq_length = a;
                }
                else { // c<=a
                    max_sq_length = a;
                    second_max_sq_length = ( b < c ? c : b );
                }
            }

            q.second = 0;
            if( squared_size_bound != 0 )
            {
                //	  std::cerr << squared_size_bound << std::endl;
                q.second = max_sq_length / squared_size_bound;
                // normalized by size bound to deal
                // with size field
                if( q.size() > 1 )
                {
                    q.first = 1; // (do not compute sine)
                    return Mesh_2::IMPERATIVELY_BAD;
                }
            }

            Compute_area_2 area_2 = traits.compute_area_2_object();

            double area = 2*CGAL::to_double(area_2(pa, pb, pc));

            q.first = (area * area) / (max_sq_length * second_max_sq_length); // (sine)

            if( q.sine() < this->B )
                return Mesh_2::BAD;
            else
                return Mesh_2::NOT_BAD;
        }
    };

    Is_bad is_bad_object() const
    { return Is_bad(this->bound(), size_bound(), 
            this->traits /* from the bad class */); }
};


typedef My_delaunay_mesh_size_criteria_2<CDT> Criteria;

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
    std::list<Point> list_of_seeds;
    // fibre boundaries
    for (int ii=0; ii<fibres.size(); ++ii) {
        Vertex_handle v1,va,vb;
        const double fdelta = delta*fibre_resolution;
        const double2 origin = get<position>(fibres)[ii];
        list_of_seeds.push_back(Point(origin[0], origin[1]));
        for (int jj=0; jj<layers; ++jj) {
            const double radius = fibre_radius-jj*fdelta;
            const double dtheta_aim = fdelta/radius;
            const int n = std::ceil(2*PI/dtheta_aim);
            const double dtheta = 2*PI/n;
            for (int kk=0; kk<n; ++kk) {
                get<position>(p) = origin + radius*double2(cos(kk*dtheta),sin(kk*dtheta));
                std::cout << "adding f at "<<get<position>(p)<<std::endl;
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
                std::cout << "adding bl at "<<get<position>(p)<<std::endl;
                if (jj==0) {
                    if (ii==0) {
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
                std::cout << "adding br at "<<get<position>(p)<<std::endl;
                if (jj==0) {
                    if (ii==0) {
                        vstart_right = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        va2 = vstart_right;
                    } else if (ii==ny) {
                        vend_right = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        cdt.insert_constraint(va2,vend_right);
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
            std::cout << "adding inlet at "<<get<position>(p)<<std::endl;
            if (jj==0) {
                if (ii==1) {
                    vstart_top = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    va1 = vstart_top;
                } else if (ii==nxb-1) {
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
            std::cout << "adding outlet at "<<get<position>(p)<<std::endl;
            if (jj==0) {
                if (ii==1) {
                    vstart_bottom = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    va1 = vstart_bottom;
                } else if (ii==nxb-1) {
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
    
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Meshing the domain..." << std::endl;
    CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
            Criteria(0.125,delta));
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Number of finite faces: " << cdt.number_of_faces() << std::endl;

    for(CDT::Finite_vertices_iterator v = cdt.finite_vertices_begin();
            v != cdt.finite_vertices_end(); ++v) {
        if (cdt.are_there_incident_constraints(v)) continue; 
        const double x = v->point().x();
        const double y = v->point().y();
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
