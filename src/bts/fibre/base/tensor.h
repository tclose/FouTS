/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 31/08/2010.

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

#ifndef __bts_fibre_base_tensor_h__
#define __bts_fibre_base_tensor_h__

#include "math/matrix.h"

namespace BTS {

	namespace Fibre {

	  namespace Base {

	    template <typename T> class TensorWriter;

      template <typename T> class Tensor : public MR::Math::Matrix<double> {

        //Public nested classes and typedefs
        public:

          typedef T RowType;
          typedef TensorWriter<Tensor> Writer;

        //Protected member variables
        protected:

          T tmpl; // Template object on which the tensor is based.

        //Public member functions
        public:

          Tensor() {}

          Tensor(const T& fibre_template)
            : MR::Math::Matrix<double>(fibre_template.vsize(), fibre_template.vsize()), tmpl (fibre_template) { tmpl.zero(); }


          Tensor (const T& fibre_template, const MR::Math::Matrix<double>::View& view)
            : MR::Math::Matrix<double>(view), tmpl(fibre_template) { tmpl.zero(); }


          Tensor(const Tensor& t)
            : MR::Math::Matrix<double>(t), tmpl(t.tmpl) {}


          Tensor& operator=(const Tensor& t)
            { MR::Math::Matrix<double>::operator=(t); tmpl = t.tmpl; return *this; }

          Tensor& operator=(const MR::Math::Matrix<double>& m)
            { assert(m.rows() == vsize()); assert(m.columns() == vsize()); MR::Math::Matrix<double>::operator=(m); return *this; }

          Tensor& operator=(double c)
            { MR::Math::Matrix<double>::operator=(c); return *this; }


          Tensor& reset(const T& fibre_template) {
            tmpl = fibre_template;
            tmpl.zero();
            MR::Math::Matrix<double>::resize(tmpl.vsize(), tmpl.vsize());
            zero();
            return *this;
          }


          Tensor& zero()
            { MR::Math::Matrix<double>::operator=(0.0); return *this; }


          Tensor& invalidate()
            { MR::Math::Matrix<double>::operator=(NAN); return *this; }


          T                              get_template() const
            { return tmpl; }


          double                         vsize() const
            { return tmpl.vsize(); }


          Tensor&                        operator+=  (double c)
            { MR::Math::Matrix<double>::operator+=(c); return *this; }


          Tensor&                        operator-=  (double c)
            { MR::Math::Matrix<double>::operator-=(c); return *this; }


          Tensor&                        operator*=  (double c)
            { MR::Math::Matrix<double>::operator*=(c); return *this; }


          Tensor&                        operator/=  (double c)
            { MR::Math::Matrix<double>::operator/=(c); return *this; }


          Tensor                         operator+ (double c) const
            { Tensor answer (*this); answer += c; return answer;}


          Tensor                         operator- (double c) const
            { Tensor answer (*this); answer -= c; return answer;}


          Tensor                         operator* (double c) const
            { Tensor answer (*this); answer *= c; return answer;}


          Tensor                         operator/ (double c) const
            { Tensor answer (*this); answer /= c; return answer;}


          Tensor&                        operator+=  (const Tensor& t)
            { MR::Math::Matrix<double>::operator+=(t); return *this; }


          Tensor&                        operator-=  (const Tensor& t)
            { MR::Math::Matrix<double>::operator+=(t); return *this; }


          Tensor                         operator+  (const Tensor& t) const
            { Tensor answer (*this); answer += t; return answer;}


          Tensor                         operator-  (const Tensor& t) const
            { Tensor answer (*this); answer -= t; return answer;}


          using MR::Math::Matrix<double>::operator+=;


      };


      template <typename T> Tensor<T>                         operator+ (double c, const Tensor<T>& t)
        { return t + c; }

      template <typename T> Tensor<T>                         operator* (double c, const Tensor<T>& t)
        { return t * c; }

	  }

	}

}


#endif /* __bts_fibre_base_tensor_h__ */
