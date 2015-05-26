/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 08/08/2010.

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

//#ifndef __bts_fibre_tractlet_set_tensor_h__
//#define __bts_fibre_tractlet_set_tensor_h__
//
//#include "bts/fibre/set/tensor.h"
//#include "bts/fibre/tractlet/set.h"
//#include "bts/fibre/tensor/writer2.h"
//
//namespace FTS {
//
//	namespace Fibre {
//
//#ifndef STRUCTURED_SET_TENSOR
//    class Tractlet::Set::Tensor : public MR::Math::Matrix<double> {
//
//          public:
//
//            typedef Fibre::Tensor::Writer2<Tractlet::Set> Writer;
//
//          public:
//
//            Tensor ()  {}
//
//            Tensor(const Set& default_value)
//              : MR::Math::Matrix<double>(default_value.vsize(), default_value.vsize()) { zero(); }
//
//            Tensor(const Tensor& t)
//              : MR::Math::Matrix<double>(t) {}
//
//            ~Tensor() {}
//
//
//            Tensor&                     operator= (const Tensor& t)
//              { this->MR::Math::Matrix<double>::operator=(t); return *this; }
//
//
//            size_t                        matrix_length() const
//              { return this->rows(); }
//
//
//            const MR::Math::Matrix<double>&    matrix() const
//              { return *this; }
//
//            Tensor&                     zero()
//              { MR::Math::Matrix<double>::operator=(0.0); return *this; }
//
//            Tensor&                     invalidate()
//              { MR::Math::Matrix<double>::operator=(NAN); return *this; }
//
//            void                        from_matrix(const MR::Math::Matrix<double>& matrix)
//              { MR::Math::Matrix<double>::operator=(matrix); }
//
//            Tensor&                     operator+=  (const MR::Math::Matrix<double>& tensor)
//              { MR::Math::Matrix<double>::operator+=(tensor); return *this; }
//
//            Tensor&                     operator-=  (const MR::Math::Matrix<double>& tensor)
//              { MR::Math::Matrix<double>::operator-=(tensor); return *this; }
//
//            Tensor&                     operator*= (const MR::Math::Matrix<double>& tensor)
//              { MR::Math::Matrix<double>::operator*=(tensor); return *this; }
//
//            Tensor&                     operator*= (double c)
//              { MR::Math::Matrix<double>::operator*=(c); return *this; }
//
//            Tensor&                     operator/= (const MR::Math::Matrix<double>& tensor)
//              { MR::Math::Matrix<double>::operator/=(tensor); return *this; }
//
//            Tensor&                     operator/= (double c)
//              { MR::Math::Matrix<double>::operator/=(c); return *this; }
//
//            Tensor                      operator+  (const MR::Math::Matrix<double>& t) const
//              { Tensor answer (*this); answer += t; return answer; }
//
//            Tensor                      operator-  (const MR::Math::Matrix<double>& t) const
//              { Tensor answer (*this); answer -= t; return answer; }
//
//            Tensor                      operator* (const MR::Math::Matrix<double>& t) const
//              { Tensor answer (*this); answer *= t; return answer; }
//
//            Tensor                      operator* (double c) const
//              { Tensor answer (*this); answer *= c; return answer; }
//
//            Tensor                      operator/ (const MR::Math::Matrix<double>& t) const
//              { Tensor answer (*this); answer /= t; return answer; }
//
//            Tensor                      operator/ (double c) const
//              { Tensor answer (*this); answer /= c; return answer; }
//
//        };
//
////        std::ostream&      operator<<(std::ostream& stream, const Strand::Set::Tensor& tensor) {
////
////         stream << "base_intensity: " << tensor.base_intensity << std::endl;
////
////         for (size_t elem_i = 0; elem_i < tensor.size(); elem_i++)
////           stream << tensor[elem_i] << std::endl;
////
////         return stream;
////
////        }
//
//#else
//	  class Tractlet::Set::Tensor : public Tractlet::Set::Tensor_tpl<Tractlet::Set> {
//
//      public:
//
//      typedef Fibre::Tensor::Writer<Tractlet::Set> Writer;
//
//      public:
//
//      Tensor ()  {}
//
//      Tensor(const Set& default_value)
//       : Tensor_tpl<Set>(default_value) {}
//
//      Tensor(const Tensor_tpl<Set>& t)
//       : Tensor_tpl<Set>(t) {}
//
//      ~Tensor() {}
//
//
//      Tensor&                     operator= (const Tensor_tpl<Set>& t)
//       { Tensor_tpl<Set>::operator=(t); return *this; }
//
//
//      size_t                        matrix_length() const
//       { return Flattener< Set >::length(*this); }
//
//
//      MR::Math::Matrix<double>    matrix() const
//       { return Flattener< Set >::flatten(*this); }
//
//
//      Tensor                      from_matrix(const MR::Math::Matrix<double>& matrix)
//       { return Flattener< Set >::unflatten(matrix, *this); }
//
//      Tensor&                     operator+=  (const Tensor& tensor)
//       { Tensor_tpl<Set>::operator+=(tensor); return *this; }
//
//      Tensor&                     operator-=  (const Tensor& tensor)
//       { Tensor_tpl<Set>::operator-=(tensor); return *this; }
//
//      Tensor&                     operator*= (const Tensor& tensor)
//       { Tensor_tpl<Set>::operator*=(tensor); return *this; }
//
//      Tensor&                     operator*= (double c)
//       { Tensor_tpl<Set>::operator*=(c); return *this; }
//
//      Tensor&                     operator/= (const Tensor& tensor)
//       { Tensor_tpl<Set>::operator/=(tensor); return *this; }
//
//      Tensor&                     operator/= (double c)
//       { Tensor_tpl<Set>::operator/=(c); return *this; }
//
//      Tensor                      operator+  (const Tensor& t) const
//       { Tensor answer (*this); answer += t; return answer; }
//
//      Tensor                      operator-  (const Tensor& t) const
//       { Tensor answer (*this); answer -= t; return answer; }
//
//      Tensor                      operator* (const Tensor& t) const
//       { Tensor answer (*this); answer *= t; return answer; }
//
//      Tensor                      operator* (double c) const
//       { Tensor answer (*this); answer *= c; return answer; }
//
//      Tensor                      operator/ (const Tensor& t) const
//       { Tensor answer (*this); answer /= t; return answer; }
//
//      Tensor                      operator/ (double c) const
//       { Tensor answer (*this); answer /= c; return answer; }
//
//	  };
//
//    std::ostream&      operator<<(std::ostream& stream, const Tractlet::Set::Tensor& tensor) {
//
//     stream << "base_intensity: " << tensor.base_intensity << std::endl;
//
//     for (size_t elem_i = 0; elem_i < tensor.size(); elem_i++)
//       stream << tensor[elem_i] << std::endl;
//
//     return stream;
//
//    }
//
//#endif
//
//
//	}
//
//}
//
//#endif /* __bts_fibre_set_tensor_h__ */
