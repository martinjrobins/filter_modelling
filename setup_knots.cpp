#include "setup_knots.h"

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

template <class CDT, typename Fibres>
class My_delaunay_mesh_size_criteria_2 : 
    public CGAL::Delaunay_mesh_size_criteria_2<CDT>
{
protected:
  typedef typename CDT::Geom_traits Geom_traits;
  double transition_size;
  double transition_dist;
  const Fibres* fibres;
  const double fibre_radius;
  const vdouble2 min,max;

public:
  typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Base;

  My_delaunay_mesh_size_criteria_2(
                                const vdouble2 min,
                                const vdouble2 max,
                                const Fibres& fibres,
                                const double fibre_radius,
                                const double aspect_bound = 0.125, 
                                const double size_bound = 0,
                                const double transition_size = 0,
                                const double transition_dist = 0,
                                const Geom_traits& traits = Geom_traits())
    : Base(aspect_bound,size_bound,traits),
      transition_size(transition_size),transition_dist(transition_dist),
      fibres(&fibres),fibre_radius(fibre_radius),
      min(min),max(max)
    {}
            

    typedef typename Base::Quality Quality; 

    class Is_bad: public Base::Is_bad
    {
    protected:
        double transition_size;
        double transition_dist;
        const Fibres* fibres;
        const double fibre_radius;
        const vdouble2 min,max;
    public:
        typedef typename Base::Is_bad::Point_2 Point_2;

        Is_bad( const vdouble2 min,
                const vdouble2 max,
                const Fibres* fibres,
                const double fibre_radius,
                const double aspect_bound,
                const double size_bound,
                const double transition_size,
                const double transition_dist,
                const Geom_traits& traits)
            : Base::Is_bad(aspect_bound, size_bound, traits),
            transition_size(transition_size),transition_dist(transition_dist),
            fibres(fibres),fibre_radius(fibre_radius),
            min(min),max(max)
        {}


        using Base::Is_bad::operator();
        /*
        CGAL::Mesh_2::Face_badness operator()(const Quality q) const
        {
            return this->operator()(q);
        }
        */

        CGAL::Mesh_2::Face_badness operator()(const typename CDT::Face_handle& fh,
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
            const vdouble2 p_centre((pa[0]+pb[0]+pc[0])/3.0,(pa[1]+pb[1]+pc[1])/3.0);

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
            if( this->squared_size_bound != 0 )
            {
                //	  std::cerr << squared_size_bound << std::endl;
                q.second = max_sq_length / this->squared_size_bound;
                double min_dist = std::numeric_limits<double>::max();
                
                for (auto pref: {pa,pb,pc}) {
                    const vdouble2 p(pref[0],pref[1]);
                    for (typename Fibres::const_reference f: *fibres) {
                        const double dist = (get<typename Fibres::position>(f)-p).norm()-fibre_radius;
                        if (dist < min_dist) min_dist = dist;
                    }
                    for (int d = 0; d < 2; ++d) {
                        const double dist_min = std::abs(p[d]-min[d]);
                        const double dist_max = std::abs(max[d]-p[d]);
                        if (dist_min < min_dist) min_dist = dist_min;
                        if (dist_max < min_dist) min_dist = dist_max;
                    }
                }

                q.second /= std::pow((1+transition_size)/2
                            + ((1-transition_size)/2)*std::tanh((2.0/transition_dist)*(2*min_dist-transition_dist)),2);
                // normalized by size bound to deal
                // with size field
                if( q.size() > 1 )
                {
                    q.first = 1; // (do not compute sine)
                    return CGAL::Mesh_2::IMPERATIVELY_BAD;
                }
            }

            Compute_area_2 area_2 = traits.compute_area_2_object();

            double area = 2*CGAL::to_double(area_2(pa, pb, pc));

            q.first = (area * area) / (max_sq_length * second_max_sq_length); // (sine)

            if( q.sine() < this->B )
                return CGAL::Mesh_2::BAD;
            else
                return CGAL::Mesh_2::NOT_BAD;
        }
    };

    Is_bad is_bad_object() const
    { return Is_bad(min,max,fibres,fibre_radius,this->bound(), this->size_bound(), 
            transition_size, transition_dist,
            this->traits /* from the bad class */); }
};


typedef My_delaunay_mesh_size_criteria_2<CDT,ParticlesType> Criteria;

CDT setup_knots(KnotsType &knots, ParticlesType &fibres, const double fibre_radius, const double fibre_resolution_factor, const double nx, vdouble2 domain_min, vdouble2 domain_max, const double k, const bool periodic,const int nbucket) {
    CDT cdt;

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
    vdouble2 ns_buffer(L/2,L/2);

    if (periodic) {
        ns_buffer[0] = 0;
    }

    const int layers = 1;


    const double dt_adapt = (1.0/100.0)*PI/sqrt(2*k);
    //const double dt_adapt = 0.5*delta;

    typename KnotsType::value_type p;
    
    std::list<Point> list_of_seeds;
    // fibre boundaries
    std::vector<std::pair<Vertex_handle,double>> cuts_min,cuts_max;
    for (int ii=0; ii<fibres.size(); ++ii) {
        Vertex_handle v1,va,vb;
        const vdouble2 origin = get<position>(fibres)[ii];
        const double dtheta = 2*PI/nx;
        list_of_seeds.push_back(Point(origin[0], origin[1]));
        bool outside,started;
        started = false;
        for (int kk=0; kk<nx; ++kk) {
            get<position>(p) = origin + fibre_radius*vdouble2(std::cos(kk*dtheta),std::sin(kk*dtheta));
            if (get<position>(p)[0] < domain_min[0]) {
                if (!started) {
                    outside = true;
                } else if (!outside) {
                    outside = true;
                    cuts_min.push_back(std::make_pair(va,va->point()[1]));
                }
            } else if (get<position>(p)[0] > domain_max[0]) {
                if (!started) {
                    outside = true;
                } else if (!outside) {
                    outside = true;
                    cuts_max.push_back(std::make_pair(va,va->point()[1]));
                }
            } else if (started && outside) {
                outside = false;
                vb = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                va = vb;
                if (get<position>(p)[0] < 0.5*(domain_max[0]+domain_min[0])) {
                    cuts_min.push_back(std::make_pair(va,get<position>(p)[1]));
                } else {
                    cuts_max.push_back(std::make_pair(va,get<position>(p)[1]));
                }
            } else if (started && !outside) {
                outside = false;
                vb = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                cdt.insert_constraint(va,vb);
                va = vb;
            } else {
                if (kk==0) {
                    v1 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    va = v1;
                    outside = false;
                    started = true;
                } else {
                    v1 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                    va = v1;
                    outside = false;
                    started = true;
                    if (get<position>(p)[0] < 0.5*(domain_max[0]+domain_min[0])) {
                        cuts_min.push_back(std::make_pair(va,get<position>(p)[1]));
                    } else {
                        cuts_max.push_back(std::make_pair(va,get<position>(p)[1]));
                    }
                }

            }


            get<target>(p) = true;
            get<boundary>(p) = true;
            get<inlet>(p) = false;
            get<outlet>(p) = false;
            get<interior>(p) = false;
            if (!outside) knots.push_back(p);
        }
        get<position>(p) = origin + fibre_radius*vdouble2(1,0);
        if (get<position>(p)[0] > domain_max[0]) {
            if (!outside) {
                cuts_max.push_back(std::make_pair(va,va->point()[1]));
            }
        } else {
            cdt.insert_constraint(va,v1);
        }
    }

    std::sort(std::begin(cuts_min),std::end(cuts_min),
            [](std::pair<Vertex_handle,double> a,std::pair<Vertex_handle,double> b) {
                return a.second < b.second;
                });
    std::sort(std::begin(cuts_max),std::end(cuts_max),
            [](std::pair<Vertex_handle,double> a,std::pair<Vertex_handle,double> b) {
                return a.second < b.second;
                });

    for (auto &i: cuts_min) {
        std::cout << "min cut at "<<i.second<< std::endl;
    }
    for (auto &i: cuts_max) {
        std::cout << "max cut at "<<i.second<< std::endl;
    }

    Vertex_handle vstart_left,vend_left,vstart_right,vend_right;
    Vertex_handle vstart_top,vend_top,vstart_bottom,vend_bottom;
    Vertex_handle va1,vb1,va2,vb2;
    if (!periodic) {
        const int ny = (domain_max[1]-domain_min[1])/(delta*fibre_resolution);
        const double deltay = (domain_max[1]-domain_min[1])/ny;
        bool in_fibre_min = false;
        bool in_fibre_max = false;
        int fibre_index_min = 0;
        int fibre_index_max = 0;
        for (int ii=1-layers; ii<ny+layers; ++ii) {
            for (int jj=0; jj<layers; ++jj) {
                // boundary - left
                get<position>(p) = vdouble2(domain_min[0]-jj*deltay,
                        domain_min[1]+ii*deltay);
                if (jj==0) {
                    if (ii==0) {
                        vstart_left = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        va1 = vstart_left;
                    } else if (ii==ny) {
                        vend_left = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        cdt.insert_constraint(va1,vend_left);
                    } else {
                        if (fibre_index_min < cuts_min.size() 
                            && get<position>(p)[1] > cuts_min[fibre_index_min].second) {
                            if (!in_fibre_min) {
                                in_fibre_min = true;
                                cdt.insert_constraint(va1,cuts_min[fibre_index_min].first);
                                ++fibre_index_min;
                            } else {
                                in_fibre_min = false;
                                vb1 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                                cdt.insert_constraint(cuts_min[fibre_index_min].first,vb1);
                                va1 = vb1;
                                ++fibre_index_min;
                            }
                        } else if (!in_fibre_min) {
                            vb1 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                            cdt.insert_constraint(va1,vb1);
                            va1 = vb1;
                        }
                    }
                }


                if ((jj==0) && (ii >= 0) && (ii <= ny)) {
                    get<target>(p) = true;
                } else {
                    get<target>(p) = false;
                }
                get<boundary>(p) = false;
                get<inlet>(p) = false;
                get<interior>(p) = false;
                get<outlet>(p) = false;
                if (!in_fibre_min) knots.push_back(p);

                // boundary - right
                get<position>(p) = vdouble2(domain_max[0]+jj*deltay,
                        domain_min[1]+ii*deltay);
                if (jj==0) {
                    if (ii==0) {
                        vstart_right = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        va2 = vstart_right;
                    } else if (ii==ny) {
                        vend_right = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                        cdt.insert_constraint(va2,vend_right);
                    } else {

                        if (fibre_index_max < cuts_max.size() 
                            && get<position>(p)[1] > cuts_max[fibre_index_max].second) {
                            if (!in_fibre_max) {
                                in_fibre_max = true;
                                cdt.insert_constraint(cuts_max[fibre_index_max].first,va2);
                                ++fibre_index_max;
                            } else {
                                in_fibre_max = false;
                                vb2 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                                cdt.insert_constraint(vb2,cuts_max[fibre_index_max].first);
                                va2 = vb2;
                                ++fibre_index_max;
                            }
                        } else if (!in_fibre_max) {
                            vb2 = cdt.insert(Point(get<position>(p)[0],get<position>(p)[1]));
                            cdt.insert_constraint(va2,vb2);
                            va2 = vb2;
                        }
                    }
                }
                if ((jj==0) && (ii >= 0) && (ii <= ny)) {
                    get<target>(p) = true;
                } else {
                    get<target>(p) = false;
                }
                get<boundary>(p) = false;
                get<inlet>(p) = false;
                get<interior>(p) = false;
                get<outlet>(p) = false;
                if (!in_fibre_max) knots.push_back(p);
            }
        }
    }

    const int nxb = (domain_max[0]-domain_min[0])/(delta*fibre_resolution);
    const double deltax = (domain_max[0]-domain_min[0])/nxb;
    for (int ii=1; ii<nxb; ++ii) {
        for (int jj=0; jj<layers; ++jj) {
            // inlet
            get<position>(p) = vdouble2(domain_min[0] + ii*deltax,domain_max[1]+jj*deltax);
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
            knots.push_back(p);
        }
    }

    for (int ii=1; ii<nxb; ++ii) {
        for (int jj=0; jj<layers; ++jj) {
            // outlet
            get<position>(p) = vdouble2(domain_min[0] + ii*deltax,domain_min[1]-jj*deltax);
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
            Criteria(domain_min,domain_max,fibres,fibre_radius,0.125,delta,fibre_resolution,0.3));

    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Number of finite faces: " << cdt.number_of_faces() << std::endl;

    std::cout << "Run Lloyd optimization...";
    CGAL::lloyd_optimize_mesh_2(cdt,
            CGAL::parameters::max_iteration_number = 10);
    std::cout << " done." << std::endl;
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

    for(CDT::Finite_vertices_iterator v = cdt.finite_vertices_begin();
            v != cdt.finite_vertices_end(); ++v) {
        if (cdt.are_there_incident_constraints(v)) continue; 
        const double x = v->point().x();
        const double y = v->point().y();
        get<position>(p) = vdouble2(x,y);
        get<boundary>(p) = false;
        get<inlet>(p) = false;
        get<outlet>(p) = false;
        get<target>(p) = true;
        get<interior>(p) = true;
        knots.push_back(p);
    }

    knots.init_neighbour_search(domain_min-ns_buffer,domain_max+ns_buffer,vbool2(periodic,false),nbucket);
    
    std::cout << "added "<<knots.size()<<" knots"<< std::endl;


    std::cout << "finshed adapting. Writing to vtk file init_knots"<<std::endl;

    //vtkWriteGrid("init_knots",1,knots.get_grid(true));
    //vtkWriteGrid("init_knots_fibres",1,fibres.get_grid(true));

    std::cout << "finshed writing to file."<<std::endl;

    return cdt;
}

/*
void calculate_c(KnotsType &knots, double c0, const double nx, vdouble2 domain_min, vdouble2 domain_max) {
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
    c[i] = c0;

    //vtkWriteGrid("init_knots",2,knots.get_grid(true));
    std::cout << "done calculate c."<<std::endl;
}
*/
