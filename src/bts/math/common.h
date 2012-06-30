/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 21/07/2010.

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

#ifndef __bts_math_common_h__
#define __bts_math_common_h__

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include "math/vector.h"
#include "math/matrix.h"

namespace BTS {

	namespace Math {


	  inline int        sign(double x)
      { return (x > 0) - (x < 0); }


	  template <typename T> MR::Math::Matrix<T>             outer(const MR::Math::Vector<T>& v1, const MR::Math::Vector<T>& v2) {

	    MR::Math::Matrix<T> answer(v1.size(), v2.size());

	    for (size_t row_i = 0; row_i < v1.size(); row_i++)
	      for (size_t col_i = 0; col_i < v2.size(); col_i++)
	        answer(row_i,col_i) = v1[row_i] * v2[col_i];

	    return answer;

	  }

	  inline MR::Math::Matrix<double>                       operator* (const MR::Math::Matrix<double>& matrix, double scalar) {

	    MR::Math::Matrix<double> answer (matrix);

	    return answer *= scalar;

	  }

	  inline MR::Math::Matrix<double>                       operator* (double scalar, const MR::Math::Matrix<double>& matrix) {

      return matrix * scalar;

    }

	  MR::Math::Matrix<double>                              random_rotation(gsl_rng* rand_gen, size_t ndims);


    std::vector<bool>                                     binary_string(size_t size, size_t i);



    //! Returns a single-line string representation of a MRtrix vector in a format readable by the MATLAB software package.
    std::string                                           matlab_str(const MR::Math::Vector<double>& v);


    //! Returns a single-line string representation of a MRtrix matrix in a format readable by the MATLAB software package.
    std::string                                           matlab_str(const MR::Math::Matrix<double>& m);


    //! returns the trace of the matrix.
    double                                                trace(const MR::Math::Matrix<double>& m);

	}

}

#endif /* __bts_math_math_h__ */
