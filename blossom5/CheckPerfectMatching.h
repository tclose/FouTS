#ifndef __blossom5_CheckPerfectMatching_h__
#define __blossom5_CheckPerfectMatching_h__

#include "blossom5/PerfectMatching.h"

// in the functions below, 'edges' is an array of size 2*edge_num (edge e = (edges[2*e],edges[2*e+1])),
//                         'weights' is an array of size edge_num

// checks complementary slackness conditions.
// returns 0 if success.
// returns 1 if complementary slackness conditions are violated (then the amount of violation is printed - could potentially happen for double's)
// returns 2 if the blossom tree structure is incorrect (or inconsistent with primal solution)
int CheckPerfectMatchingOptimality(int node_num, int edge_num, int* edges, int* weights,
                                   PerfectMatching* pm, PerfectMatching::REAL threshold =
                                           (PerfectMatching::REAL) (1e-10));

double ComputePerfectMatchingCost(int node_num, int edge_num, int* edges, int* weights,
                                  PerfectMatching* pm);

#endif
