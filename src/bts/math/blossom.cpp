/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Mar 20, 2011.

 This file is part of MRtrix.

 MRtrix is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 MRtrix is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

 */

#include <map>
#include <set>

#include "bts/common.h"
/*
include "bts/math/blossom.h"

include "blossom5/PerfectMatching.h"
include "blossom5/CheckPerfectMatching.h"

const double WEIGHTS_RANGE = 1e6;

namespace BTS {
    
    namespace Math {
        
        namespace Blossom {
            
            std::vector<Edge> match_pairs(const std::vector<Edge>& edges, bool check_match,
                                          bool fractional_jumpstart, int dual_greedy_update_option,
                                          double dual_LP_threshold, bool update_duals_before,
                                          bool update_duals_after, double single_tree_threshold,
                                          bool verbose) {
                
                if (!edges.size())
                    throw Exception(
                            "Edges vector passed to 'Math::Blossom::match_pairs' is empty.");
                
                struct PerfectMatching::Options options;
                
                //---------------------------------------------------------------------------------//
                // Convert weights from double to int with the smallest amount of information lost //
                //---------------------------------------------------------------------------------//
                
                std::set<size_t> nodes;
                
                for (size_t edge_i = 0; edge_i < edges.size(); ++edge_i) {
                    nodes.insert(edges[edge_i].first);
                    nodes.insert(edges[edge_i].second);
                }
                
                if (nodes.size() & 1)
                    throw Exception(
                            "Odd number of nodes (" + str(edges.size())
                            + "), run BTS::Math::Blossom::if_odd_remove_worst() first to remove worst match (in a greedy sense).");
                
                size_t mapped_node = 0;
                std::map<size_t, size_t> mapped_nodes;
                
                for (std::set<size_t>::iterator node_it = nodes.begin(); node_it != nodes.end();
                        ++node_it) {
                    mapped_nodes[*node_it] = mapped_node;
                    ++mapped_node;
                }
                
                int* edges_array = new int[edges.size() * 2];
                
                for (size_t edge_i = 0; edge_i < edges.size(); ++edge_i) {
                    edges_array[edge_i * 2] = mapped_nodes[edges[edge_i].first];
                    edges_array[edge_i * 2 + 1] = mapped_nodes[edges[edge_i].second];
                }
                
                double max_weight = std::max_element(edges.begin(), edges.end())->weight;
                double min_weight = std::min_element(edges.begin(), edges.end())->weight;
                
                double range = max_weight - min_weight;
                
                int* weights_array = new int[edges.size()];
                
                for (size_t edge_i = 0; edge_i < edges.size(); edge_i++) {
                    weights_array[edge_i] = (edges[edge_i].weight - min_weight) * WEIGHTS_RANGE
                            / range;
                }
                
                std::ofstream f("/home/tclose/data/edge_graph.txt");
                f << nodes.size() << " " << edges.size() << std::endl;
                for (size_t edge_i = 0; edge_i < edges.size(); ++edge_i)
                    f << edges_array[edge_i * 2] << " " << edges_array[edge_i * 2 + 1] << " "
                      << weights_array[edge_i] << std::endl;
                
                PerfectMatching *pm = new PerfectMatching(nodes.size(), edges.size());
                for (size_t edge_i = 0; edge_i < edges.size(); ++edge_i)
                    pm->AddEdge(edges_array[edge_i * 2], edges_array[edge_i * 2 + 1],
                            weights_array[edge_i]);
                
                pm->options = options;
                pm->Solve();
                
                if (check_match) {
                    
                    int error = CheckPerfectMatchingOptimality(nodes.size(), edges.size(),
                            edges_array, weights_array, pm);
                    
                    if (error)
                        throw Exception(
                                "Blossom pair matching failed (error code " + str(error) + ").");
                    
                }
                
                std::vector<Edge> output_pairs;
                
                for (size_t edge_i = 0; edge_i < edges.size(); ++edge_i)
                    if (pm->GetSolution(edge_i))
                        output_pairs.push_back(edges[edge_i]);
                
                delete pm;
                delete weights_array;
                delete edges_array;
                
                return output_pairs;
                
            }
            
            //If there are an odd number of nodes then remove the node with the maximum minimum weight.
            std::vector<Edge> make_even(const std::vector<Edge>& edges) {
                
                std::map<size_t, double> node_min_weights;
                
                for (size_t edge_i = 0; edge_i < edges.size(); ++edge_i) {
                    
                    double weight = edges[edge_i].weight;
                    
                    for (size_t node_i = 0; node_i < 2; ++node_i) {
                        
                        size_t node = edges[edge_i].node(node_i);
                        
                        if (!node_min_weights.count(node))
                            node_min_weights[node] = weight;
                        else {
                            if (weight < node_min_weights[node])
                                node_min_weights[node] = weight;
                        }
                        
                    }
                    
                }
                
                size_t num_nodes = node_min_weights.size();
                
                std::vector<Edge> even_edges(edges);
                
                //If odd number of nodes, remove the 'worst' one.
                if (num_nodes & 1) {
                    
                    bool even_found = false;
                    
                    std::vector<std::pair<size_t, double> > sorted_min_weights;
                    
                    for (std::map<size_t, double>::iterator min_weight_it =
                            node_min_weights.begin(); min_weight_it != node_min_weights.end();
                            ++min_weight_it)
                        sorted_min_weights.push_back(
                                std::pair<size_t, double>(min_weight_it->first,
                                        min_weight_it->second));
                    
                    std::sort(sorted_min_weights.begin(), sorted_min_weights.end(),
                            CompareSecond<size_t, double>());
                    
                    for (std::vector<std::pair<size_t, double> >::reverse_iterator mw_it =
                            sorted_min_weights.rbegin(); mw_it != sorted_min_weights.rend();
                            mw_it++) {
                        
                        size_t node = mw_it->first;
                        
                        for (std::vector<Edge>::iterator edge_it = even_edges.begin();
                                edge_it != even_edges.end();)
                            if ((edge_it->first == node) || (edge_it->second == node))
                                edge_it = even_edges.erase(edge_it);
                            else
                                ++edge_it;
                        
                        if (calc_num_nodes(even_edges) & 1)
                            even_edges = edges;
                        else {
                            even_found = true;
                            break;
                        }
                        
                    }
                    
                    if (!even_found)
                        throw Exception(
                                "An even number of nodes could not be found using the 'make_even' algorithm. Suggest manual removal.");
                    
                }
                
                return even_edges;
                
            }
            
            size_t calc_num_nodes(const std::vector<Edge>& edges) {
                
                std::set<size_t> nodes;
                
                for (std::vector<Edge>::const_iterator edge_it = edges.begin();
                        edge_it != edges.end(); ++edge_it) {
                    nodes.insert(edge_it->first);
                    nodes.insert(edge_it->second);
                }
                
                return nodes.size();
                
            }
        
        }
    
    }

}
*/

