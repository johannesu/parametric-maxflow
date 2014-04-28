// Solves max-flow in a single point

#define TIMING
#include <cstdlib>
#include <stdexcept>
#include <iostream>
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

typedef double tcaptype; // May be negative
typedef double captype; // Always positive
typedef double flowtype; //

typedef Graph<captype, tcaptype, flowtype> graphtype;
typedef graphtype* graph_ptr;

typedef double lambdatype;

void error_function(const char *msg)
{
  throw runtime_error(msg);
}


void mexFunction(int nlhs,          /* number of expected outputs */
        mxArray        *plhs[],     /* mxArray output pointer array */
        int            nrhs,    /* number of inputs */
        const mxArray  *prhs[]    /* mxArray input pointer array */)
{
  int n;

  ASSERT(nrhs == 5);
  ASSERT(nlhs == 2);

  // Read all input
  int curarg = 0;

  matrix<double> unary(prhs[curarg++]);
  matrix<double> slope1(prhs[curarg++]);
  matrix<double> slope2(prhs[curarg++]);
  matrix<double> pairwise(prhs[curarg++]);
  matrix<double> point(prhs[curarg++]);

  n = unary.numel();

  double lambda = point(0);
  double mu = point(1);

  // Creating output
  matrix<double> Flow(1);
  matrix<double> Solution(n);

  // Solves using floats
  captype flow_floats;

 // Create the graph
 auto_ptr<graphtype> graph;
 graph.reset( new graphtype(n, pairwise.numel(), error_function) );
 graph->add_node(n);

 // Adding unary and pairwise cost
 for (int i = 0; i < n; i++)
 {
     graph->add_tweights(i,
                unary(i)
             + slope1(i)*lambda
             + slope2(i)*mu, 0);
 }

 for (int i = 0; i < pairwise.M; i++)
 {
     int start = pairwise(i, 0);
     int end = pairwise(i, 1);
     captype cost = pairwise(i, 2);
     captype rev_cost = pairwise(i,3);

     graph->add_edge(start, end, cost, rev_cost);
 }

 flow_floats =  graph->maxflow(false);
 for (int i = 0; i < n; ++i)
 {
     double val = 0;
     if (graph->what_segment(i))
         val = 1;

     Solution(i) = val;
 }

  Flow(0) = flow_floats;
  plhs[0] = Solution;
  plhs[1] = Flow;
}