/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Mar 16, 2011.

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

#ifndef __bts_math_blossom_h__
#define __bts_math_blossom_h__

#include <vector>
#include "bts/common.h"

#include "bts/math/blossom/edge.h"

namespace FTS {
    
    namespace Math {
        
        namespace Blossom {
            
            std::vector<Edge> match_pairs(const std::vector<Edge>& edges, bool check_match = false,
                                          bool fractional_jumpstart = true,
                                          int dual_greedy_update_option = 0,
                                          double dual_LP_threshold = 0.0, bool update_duals_before =
                                                  false,
                                          bool update_duals_after = false,
                                          double single_tree_threshold = 1.0, bool verbose = true);
            
            std::vector<Edge> make_even(const std::vector<Edge>& edges);
            
            size_t calc_num_nodes(const std::vector<Edge>& edges);
        
        }
    }

}

#endif /* __bts_math_blossom_h__ */
