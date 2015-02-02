/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 22/05/09.

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

#ifndef __math_eigen_h__
#define __math_eigen_h__

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include "math/matrix.h"
#include "math/permutation.h"

namespace MR
{
  namespace Math
  {
    namespace Eigen
    {

      /** @addtogroup linalg
        @{ */

      /** @defgroup eigen Eigensystems
        @{ */

      //! Eigenvalue decomposition for real symmetric matrices
      /** \note the diagonal and lower part of \a A will be destroyed during the decomposition */
      template <typename T> class Symm
      {
        public:
          Symm (size_t n) : work (gsl_eigen_symm_alloc (n)) { }
          ~Symm () {
            gsl_eigen_symm_free (work);
          }

          Vector<T>& operator () (Vector<T>& eval, Matrix<T>& A) {
            assert (A.rows() == A.columns());
            assert (A.rows() == eval.size());
            assert (A.rows() == work->size);
            int status = gsl_eigen_symm (A.gsl(), eval.gsl(), work);
            if (status) throw Exception (std::string ("eigenvalue decomposition failed: ") + gsl_strerror (status));
            return (eval);
          }
        protected:
          gsl_eigen_symm_workspace* work;
      };



      //! Eigenvalue and eigenvector decomposition for real symmetric matrices
      /** \note the diagonal and lower part of \a A will be destroyed during the decomposition */
      template <typename T> class SymmV
      {
        public:
          SymmV (size_t n) : work (gsl_eigen_symmv_alloc (n)), N (n) { }
          SymmV (const SymmV& SV) : work (gsl_eigen_symmv_alloc (SV.size())), N (SV.size()) { }
          ~SymmV () {
            gsl_eigen_symmv_free (work);
          }

          Vector<T>& operator () (Vector<T>& eval, Matrix<T>& A, Matrix<T>& evec) {
            assert (A.rows() == A.columns());
            assert (A.rows() == evec.rows());
            assert (evec.rows() == evec.columns());
            assert (A.rows() == eval.size());
            assert (A.rows() == work->size);
            int status = gsl_eigen_symmv (A.gsl(), eval.gsl(), evec.gsl(), work);
            if (status) throw Exception (std::string ("eigenvalue decomposition failed: ") + gsl_strerror (status));
            return (eval);
          }

          size_t size () const {
            return N;
          }

        protected:
          gsl_eigen_symmv_workspace* work;
          const size_t N;
      };


      //! Eigenvalue sorting
      inline Vector<double>& sort (Vector<double>& eval)
      {
        gsl_sort_vector (eval.gsl());
        return (eval);
      }
      //! Eigenvalue sorting
      template <typename T> inline Vector<T>& sort (Vector<T>& eval, Matrix<T>& evec)
      {
        gsl_eigen_symmv_sort (eval.gsl(), evec.gsl(), GSL_EIGEN_SORT_VAL_ASC);
        return (eval);
      }

      /** @} */
      /** @} */
    }
  }
}

#endif







