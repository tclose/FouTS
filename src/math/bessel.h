/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 27/06/08.

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



/*

 Implementation based on the GSL (http://www.gnu.org/software/gsl/)

*/

#ifndef __math_bessel_h__
#define __math_bessel_h__

#include <limits>

#include "math/chebyshev.h"

namespace MR
{
  namespace Math
  {
    namespace Bessel
    {

      extern const double coef_aI0[];
      extern const double coef_bI0[];
      extern const double coef_cI0[];

      extern const double coef_aI1[];
      extern const double coef_bI1[];
      extern const double coef_cI1[];


      //* Compute the scaled regular modified cylindrical Bessel function of zeroth order exp(-|x|) I_0(x). */
      /** Implementation based on the GSL (http://www.gnu.org/software/gsl/) */
      template <typename T> inline T I0_scaled (const T x)
      {
        assert (x >= 0.0);
        if (x*x < 4.0*std::numeric_limits<T>::epsilon()) return (1.0-x);
        if (x <= 3.0) return (exp (-x) * (2.75 + Chebyshev::eval (coef_aI0, 11, -1.0, 1.0, x*x/4.5-1.0)));
        if (x <= 8.0) return ( (0.375 + Chebyshev::eval (coef_bI0, (sizeof (T) >4?20:13), -1.0, 1.0, (48.0/x-11.0) /5.0)) /sqrt (x));
        return ( (0.375 + Chebyshev::eval (coef_cI0, (sizeof (T) >4?21:11), -1.0, 1.0, 16.0/x-1.0)) /sqrt (x));
      }

      //* Compute the scaled regular modified cylindrical Bessel function of first order exp(-|x|) I_1(x). */
      /** Implementation based on the GSL (http://www.gnu.org/software/gsl/) */
      template <typename T> inline T I1_scaled (const T x)
      {
        assert (x >= 0.0);
        if (x == 0.0) return (0.0);
        if (x*x < 8.0*std::numeric_limits<T>::epsilon()) return (0.5*x);
        if (x <= 3.0) return (x * exp (-x) * (0.875 + Chebyshev::eval (coef_aI1, 10, -1.0, 1.0, x*x/4.5-1.0)));
        if (x <= 8.0) {
          T b = (0.375 + Chebyshev::eval (coef_bI1, (sizeof (T) >4?20:11), -1.0, 1.0, (48.0/x-11.0) /5.0)) / sqrt (x);
          return (x > 0.0 ? b : -b);
        }
        T b = (0.375 + Chebyshev::eval (coef_cI1, (sizeof (T) >4?21:9), -1.0, 1.0, 16.0/x-1.0)) / sqrt (x);
        return (x > 0.0 ? b : -b);
      }

    }
  }
}

#endif

