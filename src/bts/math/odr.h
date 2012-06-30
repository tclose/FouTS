/*
    Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

    Written by Thomas G Close, 5/05/09.

    This file is part of Bayesian Tractlet Sampling (BTS).

    BTS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BTS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BTS.  If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef __bts_math_odr_h__
#define __bts_math_odr_h__

#include "math/matrix.h"
#include "math/vector.h"

#include "bts/math/svd.h"
#include "bts/coord.h"


#include "bts/triple.h"

//#define _DEBUG_VERBOSE_

namespace BTS {
  namespace Math {

    class Line {
    
      public:
        Coord centre;
        Coord direction;
    
        friend std::ostream& operator<< (std::ostream& stream, const Line line)
        {
          std::cout << "centre = " << line.centre << std::endl;
          std::cout << "direction = " << line.direction << std::endl;
          return (stream);
        }
    
    };

    // Find the centroid of a set of points.
    inline Coord             centroid(const std::vector<Coord>& points) {
    
      Coord sum(0.0,0.0,0.0);
    
      for (std::vector<Coord>::const_iterator it = points.begin(); it != points.end(); it++)
        sum += *it;
    
      sum /= points.size();
    
      return sum;
    
    }
    
    
    // Find the orthogonal-distance regressor line through a set of points.
    inline Line            odr(const std::vector<Coord>& points) {
    
      Line line;
      
      line.centre = centroid(points);
    
      MR::Math::Matrix<double> matrix(points.size(), 3);
      
      for (size_t point_i = 0; point_i < points.size(); point_i++) {
    
        matrix(point_i,X) = points[point_i][X] - line.centre[X];
        matrix(point_i,Y) = points[point_i][Y] - line.centre[Y];
        matrix(point_i,Z) = points[point_i][Z] - line.centre[Z];                
      
      }
    
    
      USV usv = svd(matrix);
//      
//      std::cout << "V = [";      
//      for (size_t row_i = 0; row_i < usv.v.rows(); row_i++) {
//        for (size_t col_i = 0; col_i < usv.v.columns(); col_i++) 
//         std::cout << usv.v(row_i,col_i) << " "; 
//       std::cout << "; ";         
//      }
//      std::cout << "]" << std::endl;
      
      line.direction.set(usv.v(X,0), usv.v(Y,0), usv.v(Z,0));
//      line.direction.set(usv.v(max_sigma_i,X), usv.v(max_sigma_i,Y), usv.v(max_sigma_i,Z));      
      return line;
    
    }



  }
}

#endif
