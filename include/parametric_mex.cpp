//
// Solving the parametric max-flow problem for two parameters.
//
// Notes
// a: The code relies on the fact that E(0)=0 that is the zero solution has zero flow
//    1) It relies on the fact that (0,0,0) lies inside the solution polyhedron.
//    2) The bounding box is calculated with this assumption.
//
// b: The code only supports integer unary and pairwise terms.

// c: CGAL does not return the added intersection points when intersections are performed.
//    On intersection of polyhedrons CGAL recreates each face and vertex
//    that is why the solution and vertex queue are stored separate.
//
// Johannes Ulén

// MATLAB and Matrices
#define TIMING
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <memory>
#include <cmath>
#include <vector>
using namespace std;
#include "mex.h"
#include "utils/cppmatrix.h"
#include "utils/mexutils.h"
#include "utils/mextiming.h"
#include "maxflow-v3.02.src/graph.h"
#include <stdint.h>
#include <limits>
#include <list>
#include <stack>
#include <set>

// CGAL
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>

// Ability to test code with inexact kernel
#if 1
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
#else
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
#endif

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Vertex_iterator Vertex_iterator;
typedef Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;
typedef Nef_polyhedron::Halffacet Halffacet;
typedef Nef_polyhedron::Halffacet_const_iterator Halffacet_const_iterator;
typedef Nef_polyhedron::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;
typedef Nef_polyhedron::SFace_cycle_const_iterator SFace_cycle_const_iterator;
typedef Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;
typedef Nef_polyhedron::Vertex Vertex;
typedef Nef_polyhedron::Plane_3  Plane;
typedef Nef_polyhedron::Point_3  Point;
typedef Nef_polyhedron::Line_3 Line;
typedef Kernel::Direction_3  Direction;
typedef Nef_polyhedron::Vector_3 Vector;
typedef Nef_polyhedron::Object_handle Object_handle;
typedef Nef_polyhedron::Vertex_const_handle Vertex_const_handle;
typedef Nef_polyhedron::Volume_const_handle Volume_const_handle;
typedef Nef_polyhedron::Halffacet_const_iterator Halffacet_const_iterator;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

typedef Polyhedron::Halfedge_handle Halfedge_handle;

// Datatypes used in CGAL, field and ring.
typedef Kernel::FT  FT;
typedef Kernel::RT  RT;

typedef double tcaptype; // May be negative
typedef double captype; // Always positive
typedef double flowtype; //

typedef matrix<double> matrixtype;
typedef Graph<captype, tcaptype, flowtype> graphtype;
typedef graphtype* graph_ptr;

typedef double lambdatype;

void error_function(const char *msg) 
{
    throw runtime_error(msg);
}


struct Interval
{
    flowtype min_flow;
    flowtype max_flow;
    
    flowtype min_lambda;
    flowtype max_lambda;
    flowtype min_mu;
    flowtype max_mu;
};

class breakpoint
{
    private:
        mutable bool verified;

    public:
        Point p;

        breakpoint()
        {
            verified = false;
        }
        
        breakpoint(flowtype x, flowtype y, flowtype z)
        {
            p = Point(x,y,z);
            verified = false;
        }
        
        breakpoint(Point external)
        {
            p = external;
            verified = false;
        }

        bool operator< (const breakpoint &rhs) const
        {
            return (p < rhs.p);
        }

        bool operator== (const breakpoint &rhs) const
        {
            return (p == rhs.p);
        }

        void verify() const
        {
            verified = true;
        }

        bool is_verified() const
        {
            return verified;
        }
        
        Point get_point() const
        {
            return p;
        }
        
        FT x() const { return p.x(); }
        FT y() const { return p.y(); }
        FT z() const { return p.z(); }
};

class solution_facet
{
    public:
        Plane plane_solution;
        vector<bool> solution;
    
        
    solution_facet(Plane plane_external)
    {
        plane_solution = plane_external;
    }
    
    solution_facet(Plane plane_external, vector<bool> solution_external)
    {
        plane_solution = plane_external;
        solution = solution_external;
    }
    
    bool operator==(const solution_facet &rhs) const
    {
        return (plane_solution == rhs.plane_solution);
    }

    // Note ax+by+cd+d and f(a*x+by+cd+d) is the same
    // tanget plane. For this application this will not be a problem.
    bool operator< (const solution_facet &rhs) const
    {
        if (an() < rhs.an()) { return true; }
        if (an() > rhs.an()) { return false; }

        if (bn() < rhs.bn()) { return true; }
        if (bn() > rhs.bn()) { return false; }

        if (cn() < rhs.cn()) { return true; }
        if (cn() > rhs.cn()) { return false; }

        if (dn() < rhs.dn()) { return true; }
        if (dn() > rhs.d()) { return false; }

        return false;
    }

    RT a() const {return plane_solution.a();}
    RT b() const {return plane_solution.b();}
    RT c() const {return plane_solution.c();}
    RT d() const {return plane_solution.d();}
    
    RT an() const {return norm()*a();};
    RT bn() const {return norm()*b();};
    RT cn() const {return norm()*c();};
    RT dn() const {return norm()*d();};
    
    RT norm() const
    {
        RT normalizer;
       
        if (a() > 0)
        {
            normalizer = 1;
        } else
        {
            normalizer = -1;
        }

        return normalizer;
    }
    
    int solution_size() const
    {
        return solution.size();
    }
    
    vector<bool>::const_iterator begin() const
    {
        return solution.begin();
    }
    
    vector<bool>::const_iterator end() const
    {
        return solution.end();
    }
    
};

matrixtype *unary;
matrixtype *slope1;
matrixtype *slope2;
matrixtype *pairwise;

int n;

class Polytope
{
public:
    Nef_polyhedron Solution;
    Nef_polyhedron::Vertex_const_iterator  verticeIterator;
    std::set<breakpoint> breakpoints;
    std::set<solution_facet> solution_facets;
    Interval boundary;
    
    flowtype max_flow;
    flowtype min_flow;   
    
    // Search through all vertices
    bool inside_polytope(Point p)
    {
        bool inside = false;
         
        Vertex_const_iterator vi;
        CGAL_forall_vertices(vi, Solution) 
        {
           Point vp = vi->point();
           
           if ((p.x() == vp.x()) && (p.y() == vp.y()))
           {
             
               if (p.z() < vp.z())
               {
                   inside = true;
                   break;
               }
           }
        }
             
        return inside;
    }
    
    // Searches through list and validate all points along this ray
    void verify_point(const breakpoint * p)
    {
        std::set<breakpoint>::iterator it;
        for (it = breakpoints.begin();
             it != breakpoints.end();
             it++)
        {
            if ((it->x() == p->x()) && (it->y() == p->y()))
            {
                it->verify();
            }
        }
    }
        
    // Main loop: 
    // return false if there a re no more iterations to perform.
    bool iterate()
    {
        if (!any_unverified())
            return false;

        const breakpoint *bp_pointer = get_unverified();        
        Evaluate(bp_pointer);
        return true;
    }

    Polytope(Interval I)
    {
        boundary = I;
        
        // Over estimation
        max_flow = 0;
        min_flow = 0;
        
        for (int i = 0; i < n; i++)
        {
            min_flow += min(unary->get_value(i),0.0);
            max_flow += max(unary->get_value(i),0.0);
            
            if (slope1->get_value(i) < 0)
            {
                min_flow += min(I.max_lambda * slope1->get_value(i),0.0);
                max_flow += max(I.min_lambda * slope1->get_value(i),0.0);
            } else
            {
                min_flow += min(I.min_lambda * slope1->get_value(i),0.0);
                max_flow += max(I.max_lambda * slope1->get_value(i),0.0);;
            }
            
            if (slope2->get_value(i) < 0)
            {
                min_flow += min(I.max_lambda * slope2->get_value(i),0.0);
                max_flow += max(I.min_lambda * slope2->get_value(i),0.0);
            } else
            {
                min_flow += min(I.min_lambda * slope2->get_value(i),0.0);
                max_flow += max(I.max_lambda * slope2->get_value(i),0.0);;
            }
        }

        I.max_flow = max_flow;
        I.min_flow = min_flow;

        
        Point P000(I.min_lambda,I.min_mu,I.min_flow);
        Point P100(I.max_lambda,I.min_mu,I.min_flow);
        Point P110(I.max_lambda,I.max_mu,I.min_flow);
        Point P111(I.max_lambda,I.max_mu,I.max_flow);
        Point P101(I.max_lambda,I.min_mu,I.max_flow);
        Point P010(I.min_lambda,I.max_mu,I.min_flow);
        Point P001(I.min_lambda,I.min_mu,I.max_flow);
        Point P011(I.min_lambda,I.max_mu,I.max_flow);
        
        // Most effective way to construct a cube(?).
        Polyhedron P;
        
        Halfedge_handle h = P.make_tetrahedron( P100,
                P001,
                P000,
                P010 );
        Halfedge_handle g = h->next()->opposite()->next();
        P.split_edge( h->next());
        P.split_edge( g->next());
        P.split_edge( g);
        h->next()->vertex()->point()     = P101;
        g->next()->vertex()->point()     = P011;
        g->opposite()->vertex()->point() = P110;
        Halfedge_handle f = P.split_facet( g->next(),
                g->next()->next()->next());
        Halfedge_handle e = P.split_edge( f);
        e->vertex()->point() = P111;
        P.split_facet( e, f->next()->next());
                
        Nef_polyhedron Cube(P);
        Solution = Cube;
        
        add_vertices();
        verify_all_vertices();     
        add_facets();
        
        // Setup the initial polyhedron evaluated at the boundary.
        const breakpoint p1(I.min_lambda,I.min_mu,0);       
        Evaluate(&p1);
    }

    
    bool point_on_polyhedron(const Point p)
    {
        Vertex_const_iterator vi;
        CGAL_forall_vertices(vi, Solution)
        {
            Point pc = vi->point();
            
            if ( (pc.x() == p.x()) && (pc.y() == p.y()) )
                return true;
            
        }
        return false;
    }
    
       
    // Calculate a solution in point (lambda,mu) 
    // cut the corresponding plane with the polyhedron
    // a*x+b*y + d = flow
    // If we let the flow be z we 
    // get the affine coordinates as
    // a*x + b*y - z + d = 0.
    // -----
    // For increased speed we use double in the max flow computations
    // and the calculate the solutions flow in FT.
    // -----
    void Evaluate(const breakpoint *p)
    {
          // Check if point still lies on the polyhedron
        if (point_on_polyhedron(p->get_point()))
        {
            
            FT lambda_FT = p->x();
            FT mu_FT =  p->y();

            // mexPrintf("Evaluting point  (%f,%f,%f) ", 
            // CGAL::to_double(p->x()),
            // CGAL::to_double(p->y()),
            // CGAL::to_double(p->z()));
            captype flow_rounded;

            RT a_RT(0), b_RT(0);
            FT d_FT(0), flow_FT(0);

            vector<bool> local_solution;
            
            // Solves using floats
           captype lambda = CGAL::to_double(lambda_FT);
           captype mu = CGAL::to_double(mu_FT);

           
           // Create the graph
           auto_ptr<graphtype> graph;
           graph.reset( new graphtype(n, pairwise->numel(), error_function) );
           graph->add_node(n);
                
           // Adding Unary and pairwise cost
           for (int i = 0; i < n; i++)
           {
               graph->add_tweights(i, 
                          unary->get_value(i) 
                       +  slope1->get_value(i)*lambda
                       + slope2->get_value(i)*mu, 0);
           }
                        
           for (int i = 0; i < pairwise->M; i++)
           {
               int start = pairwise->get_value(i, 0);
               int end = pairwise->get_value(i, 1);
               captype cost = pairwise->get_value(i, 2);
               captype rev_cost = pairwise->get_value(i,3);

               graph->add_edge(start, end, cost, rev_cost);
           }
                                        
           // Solve the max-flow instance.
           flow_rounded =  graph->maxflow(false);

           for (int i = 0; i < n; ++i)
           {
               if (graph->what_segment(i))
               {
                   a_RT += slope1->get_value(i);
                   b_RT += slope2->get_value(i);
               }

               local_solution.push_back(graph->what_segment(i));
           }          
  
           // Calculate exact flow
           for (int i = 0; i < n; i++)
           {
               if (local_solution.at(i))
               {
                flow_FT +=  unary->get_value(i);
               }
           }
           
           flow_FT +=  FT(a_RT)*lambda_FT;
           flow_FT +=  FT(b_RT)*mu_FT;

           for (int i = 0; i < pairwise->M; i++)
           {
               int start = pairwise->get_value(i, 0);
               int end = pairwise->get_value(i, 1);


               if ( (!local_solution.at(start)) && (local_solution.at(end)) )
               {
                   flow_FT += pairwise->get_value(i, 2);
               }
               
               if ( (local_solution.at(start)) && (!local_solution.at(end)) )
               {
                   flow_FT += pairwise->get_value(i, 3);
               }
           }

            // Debug:
            // mexPrintf("Flow % f , Numerical Diff: %f \n", flow_rounded, flow_rounded-CGAL::to_double(flow_FT));
            int nof = Solution.number_of_vertices();
 
            Point new_point(lambda_FT,mu_FT, flow_FT);
            
            if (inside_polytope(new_point))
            {
                Direction orthogonal(-a_RT,-b_RT,1);
                Plane cutting_plane(new_point,orthogonal);
                
                startTime();
                Solution =  Solution.intersection(cutting_plane, Nef_polyhedron::CLOSED_HALFSPACE);               

                add_vertices();
                add_facets(local_solution);
            } 

        } // end to see if points need to be evaluated
        
        verify_point(p);
    }
    
    // debug function   
    void VERBOSE_polyhedron()
    {
        mexPrintf("The following points make up the boundary of the polyhedron \n");
        Vertex_const_iterator vi;
        CGAL_forall_vertices(vi, Solution)
        {
            Point p = vi->point();
            mexPrintf("(x,y,z) = (%f,%f,%f) \n", CGAL::to_double(p.x()),
                    CGAL::to_double(p.y()),
                    CGAL::to_double(p.z()) );
        }
    }
        
    void add_vertices()
    {
        Vertex_const_iterator vi;
        CGAL_forall_vertices(vi, Solution)
        {
            breakpoints.insert(vi->point());
        }
         
        // //Debug
        // mexPrintf("Listing the points in the que \n");
        
        // std:set<breakpoint>::iterator it;
        
        // for (it = breakpoints.begin();
        //      it != breakpoints.end();
        //      it++)
        // {
        //     if (!it->is_verified())
        //     {
            
        //     mexPrintf("Point = (%f,%f,%f) \n", CGAL::to_double(it->x()), 
        //                                        CGAL::to_double(it->y()), 
        //                                        CGAL::to_double(it->z()));
        //     }
        // }
    }
    
    // Default solution only zeros
    void add_facets()
    {
        vector<bool> new_solution;
        for (int i = 0; i < n; ++i)
        {
            new_solution.push_back(0);
        }
        
        add_facets(new_solution);
    }
    
    void add_facets(vector<bool> new_solution)
    {
        Halffacet_const_iterator hi;
        set<solution_facet>::iterator it;
        
        int nof = solution_facets.size();
                
        CGAL_forall_facets(hi, Solution)
        {                
           solution_facet test(hi->plane()); 
           it = solution_facets.find(test);

           if (it == solution_facets.end())
           {
               solution_facet full_solution(hi->plane(), new_solution);
               solution_facets.insert(full_solution);       
           }
        }
        
        if (nof == solution_facets.size())
            mexWarnMsgTxt("No solution added \n");
                   
        return;
    }
    
    // Only used for debugging
    void echo_unverfied()
    {
        std::set<breakpoint>::iterator it;
        int num_unverified = 0;

        for (it = breakpoints.begin();
             it != breakpoints.end();
             it++)
        {
            if (!it->is_verified())
            {
                mexPrintf("DEBUG: Unverifed point (%f,%f,%f) \n", 
                            CGAL::to_double(it->x()),
                            CGAL::to_double(it->y()),
                            CGAL::to_double(it->z()));
            }
        }
        
        return;        
    }
    
    
    // Add all vertices in the polyhedron and verify them.
    void verify_all_vertices()
    {
        std::set<breakpoint>::iterator it;
        int num_unverified = 0;

        for (it = breakpoints.begin();
             it != breakpoints.end();
             it++)
        {
            it->verify();
        }

        return; 
    }

    int number_of_unverified()
    {
        std::set<breakpoint>::iterator it;
        int num_unverified = 0;

        for (it = breakpoints.begin();
             it != breakpoints.end();
             it++)
        {
            if (!it->is_verified())
            {
                ++num_unverified;
            }
        }
        
        return num_unverified;
    }

    // debug function
    void echo_facets()
    {
        Halffacet_const_iterator hci;
        
        mexPrintf("All plane equations \n");
        
        CGAL_forall_facets(hci,Solution)
        {
            mexPrintf("Plane (%f,%f,%f,%f) \n", 
                        CGAL::to_double(hci->plane().a()),
                        CGAL::to_double(hci->plane().b()),
                        CGAL::to_double(hci->plane().c()),
                        CGAL::to_double(hci->plane().d()));
        }
                
    }

    bool any_unverified()
    {
        std::set<breakpoint>::iterator it;
        for (it = breakpoints.begin();
             it != breakpoints.end();
             it++)
        {
            if (!it->is_verified())
                return true;
        }  
       
       return false;
    }

    // s a pointer to one unverified point
    const breakpoint* get_unverified()
    {
        std::set<breakpoint>::iterator it;
        for (it = breakpoints.begin();
             it != breakpoints.end();
             it++)
        {
            if (!it->is_verified())
                return &(*it);
        }  

        mexErrMsgTxt("This line should not be reached \n");
    }

    // No built in functionality for this?
    int number_of_vertices_on_facet(const Halffacet * f)
    {
        int num = 0;

        
        for(Halffacet_cycle_const_iterator fc = f->facet_cycles_begin();
            fc != f->facet_cycles_end();
            ++fc)
        {
            if ( fc.is_shalfedge() )
            {
            Nef_polyhedron::SHalfedge_const_handle h = fc;
            Nef_polyhedron::SHalfedge_around_facet_const_circulator hc(h), he(hc);

            // all vertex coordinates in facet cycle
            CGAL_For_all(hc,he)
              {
                 ++num;
              }
            }
        }

        return num;
    }

    mxArray* get_plane(const solution_facet* sf)
    {
        // Preallocate a MATLAB matrix
        mxArray* plane_matrix = mxCreateDoubleMatrix(1, 4, mxREAL);
        double *v_ptr = mxGetPr(plane_matrix);
        
        v_ptr[0] = CGAL::to_double(sf->a());
        v_ptr[1] = CGAL::to_double(sf->b());
        v_ptr[2] = CGAL::to_double(sf->c());
        v_ptr[3] = CGAL::to_double(sf->d());
       
        return plane_matrix;
    }
    
    // Create with the solution assoicated with halffacet
    mxArray* get_solution(const solution_facet* sf)
    {
        // Preallocate a MATLAB matrix
        mxArray* solution_matrix = mxCreateDoubleMatrix(1, n, mxREAL);
        double *v_ptr = mxGetPr(solution_matrix);
          
        // Write the solution the Matlab format
        int index = 0; 
        
        vector<bool>::const_iterator vi;
       
        for (vi = sf->begin();
             vi != sf->end();
             vi++)
        {
            double val = 0;
            if (*vi) val = 1,
                    
            v_ptr[index]  = val;
            ++index;
        }

       
        
        return solution_matrix;
    }
    
    // Create a matrix consisting of all vertices for a facet
    mxArray* get_vertices(const Halffacet * f)
    {
        // Preallocate a MATLAB matrix
        int number_of_vertices = number_of_vertices_on_facet(f);
        mxArray* vertices_matrix = mxCreateDoubleMatrix(3, number_of_vertices, mxREAL);
        double *v_ptr = mxGetPr(vertices_matrix);

        int index = 0;

        
//         Vertex_const_iterator vc;
//         CGAL_forall_vertices(vc, f)
//         {
//             Point p = vc->point();
//             v_ptr[index]  = CGAL::to_double(p.x());
//             v_ptr[index+1] = CGAL::to_double(p.y());
//             v_ptr[index+2] = CGAL::to_double(p.z());
//             index += 3;
//         }
        
         for(Halffacet_cycle_const_iterator fc = f->facet_cycles_begin();
            fc != f->facet_cycles_end();
            ++fc)
          {

            if ( fc.is_shalfedge() )
                {
                Nef_polyhedron::SHalfedge_const_handle h = fc;
                Nef_polyhedron::SHalfedge_around_facet_const_circulator hc(h), he(hc);

                // all vertex coordinates in facet cycle
                CGAL_For_all(hc,he)
                  {
                    Nef_polyhedron::SVertex_const_handle v = hc->source();
                    const Nef_polyhedron::Point_3& point = v->source()->point();

                    v_ptr[index]  = CGAL::to_double(point.x());
                    v_ptr[index+1] = CGAL::to_double(point.y());
                    v_ptr[index+2] = CGAL::to_double(point.z());
                    index += 3;
                  }

              } else 
              {
                mexPrintf("not halfedge? \n");
              }
            
          }

        return vertices_matrix;
    }

    bool is_intial_box(const Halffacet * f)
    {
        Plane p = f->plane();
        
        double a = CGAL::to_double(p.a());
        double b = CGAL::to_double(p.b());
        double c = CGAL::to_double(p.c());
        double d = CGAL::to_double(p.d());
                         
        double pa = abs(a);
        double pb = abs(b);
        double pc = abs(c);
        double pd = abs(d);
        
        // Initial infinitesimal box
        if (pd ==1)
        {
            if ((pa == 1) && (pb == 0) && (pc == 0)) return true;
            if ((pa == 0) && (pb == 1) && (pc == 0)) return true;
            if ((pa == 0) && (pb == 0) && (pc == 1)) return true;
        }
 
        // Sides
        if (pc == 0)
        {
            if ((pa == 0) && (pb == 1))
            {
                if (d = -1) return true;
                if (d = 1) return true;
                if (d = boundary.min_lambda) return true;
                if (d = boundary.max_lambda) return true;
                if (d = boundary.min_mu) return true;
                if (d = boundary.max_mu) return true;
            }
            if ((pa == 1) && (pb == 0)) {
                if (d = -1) return true;
                if (d = 1) return true;
                if (d = boundary.min_lambda) return true;
                if (d = boundary.max_lambda) return true;
                if (d = boundary.min_mu) return true;
                if (d = boundary.max_mu) return true;
            }
        }
        
        // Bottom
        if ((pa == 0) && (pb == 0) && (pc == 1) && (pd == -min_flow)) return true;
          
        return false;
    }
   
    // Counts the number of solutions
    // by going through all facets in the
    // current polyhedron and match them with
    // any entries in the solution set.
    int number_of_solutions()
    {     
        int num_solutions = 0;
        Halffacet_const_iterator f;
        
        CGAL_forall_facets(f, Solution)
        {
            if (is_intial_box(&(*f))) continue;
            
            set<solution_facet>::iterator it;
            solution_facet test(f->plane());
            
            it = solution_facets.find(test);
            
            if (it != solution_facets.end()) num_solutions++;
        }

        return num_solutions;
    }
    
    
    mxArray* get_complete_solution()
    {              
        if (!Solution.is_simple()) {
            mexWarnMsgTxt("Solution polyhedron is not simple. \n");
        }

        const char * struct_names [] = {(char *)"polygon",(char *)"solution",(char *)"tanget_plane"};
        mxArray * struct_ptr = mxCreateStructMatrix(number_of_solutions(), 1, 3, struct_names);

        // Loop over each edge on each Halffacet
        int s = 0;

        Halffacet_const_iterator f;
        CGAL_forall_facets(f, Solution) 
        {
            if (is_intial_box(&(*f))) continue;

            //Check that the solution is nonempty
            const solution_facet* sfp;
            if ( sfp = find_solution_facet(&(*f)) )
            {    
                mxSetField(struct_ptr, s, "polygon", get_vertices(&(*f)));
                mxSetField(struct_ptr, s, "solution", get_solution(sfp));
                mxSetField(struct_ptr, s, "tanget_plane", get_plane(sfp));               
                s++;    
            }
         }
          
         return struct_ptr;
    }
    
    const solution_facet* find_solution_facet(const Halffacet * f)
    {
        // Find the solution based on plane equation
        set<solution_facet>::iterator it;
        solution_facet test(f->plane()); 
        
        it = solution_facets.find(test);        
        if (it == solution_facets.end())
        {
            
//          mexWarnMsgTxt("Solution for this facet not stored. \n");
            mexPrintf("Error \n");

                        
            mexPrintf("Test plane has equation (a,b,c,d) = (%f,%f,%f,%f) \n", CGAL::to_double(f->plane().a()),
                                                              CGAL::to_double(f->plane().b()),
                                                              CGAL::to_double(f->plane().c()),
                                                              CGAL::to_double(f->plane().d()));
            
            for (it = solution_facets.begin();
                 it != solution_facets.end();
                 ++it)
            {
               mexPrintf("A plane in the list (a,b,c,d) = (%f,%f,%f,%f) \n", CGAL::to_double(it->a()),
                                                              CGAL::to_double(it->b()),
                                                              CGAL::to_double(it->c()),
                                                              CGAL::to_double(it->d()));              
            }
            
            // returning whatever
            return &(*solution_facets.begin());
            
            // return false;
            
            
        } 
        
        return &(*it);
    }
};

void mexFunction(int nlhs,          /* number of expected outputs */
        mxArray        *plhs[],     /* mxArray output pointer array */
        int            nrhs,        /* number of inputs */
        const mxArray  *prhs[]      /* mxArray input pointer array */)
{

    ASSERT(nrhs >= 4);
    ASSERT(nrhs <= 6);
    
    ASSERT(nlhs == 1);

    // Read all input
    int curarg = 0;
    unary = new matrixtype(prhs[curarg++]);
    slope1 = new matrixtype(prhs[curarg++]);
    slope2 = new matrixtype(prhs[curarg++]);
    pairwise = new matrixtype(prhs[curarg++]);
    
    n = unary->numel();
    
    Interval I;
        
    if (nrhs >= 4)
    {
        matrix<double> given_interval(prhs[curarg++]);
        
        I.min_lambda = given_interval(0);
        I.max_lambda = given_interval(1);
        I.min_mu = given_interval(2);
        I.max_mu = given_interval(3);
        
    } else 
    {
        I.min_lambda = -100;
        I.max_lambda = 100;
        I.min_mu = -100;
        I.max_mu = 100;
 
    }
    
    int max_iter;
    
    if (nrhs == 6)
    {
        matrix<double> given_max_iter(prhs[curarg++]);
        max_iter = given_max_iter(0);
    }   else {
        max_iter = 1000;
    }
    
    // Solve
    Polytope P(I);    

    // Iterations
    for (int itr = 1; itr < max_iter; itr++)
    {      
       if (!P.iterate())
           break;

       if (itr == max_iter -1)
           mexWarnMsgTxt("Reached maximum number of iterations \n");
     }
    
    // Output
    plhs[0] = P.get_complete_solution();
}