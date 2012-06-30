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

#ifndef __bts_math_blossom_edge_h__
#define __bts_math_blossom_edge_h__

namespace BTS {

	namespace Math {

	  namespace Blossom {

	    class Edge {

	      public:

	        Edge() {}

	        Edge(size_t first, size_t second, double weight)
	          : first(first), second(second), weight(weight) {}

	      public:

	        size_t first;
	        size_t second;
	        double weight;

                bool  has_node(size_t index) const
                  { return (first == index) || (second == index); }

	        size_t& node(size_t index) { assert(index < 2); return !index ? first : second; }
                const size_t& node(size_t index) const { assert(index < 2); return !index ? first : second; }

                bool operator< (const Edge& e) const { return weight < e.weight; }

	    };


	    inline std::ostream& operator<<(std::ostream& stream, const Edge& e) {

	      stream << "{" << e.first << "," << e.second << "} = " << e.weight;

	      return stream;

	    }
	  }

	}

}

#endif /*  */
