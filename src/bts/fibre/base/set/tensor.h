/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 11/08/2010.

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

#ifndef __bts_fibre_set_tensor_h__
#define __bts_fibre_set_tensor_h__

#include "bts/fibre/base/set.h"

//
//namespace BTS {
//
//	namespace Fibre {
//
//	  namespace Base {
//
//      template <typename T> template <typename U> class Set<T>::Tensor_tpl : public std::map<std::string,std::string>::Intrinsic_tpl<U> {
//
//        public:
//
//          typedef typename U::template Tensor_tpl<U> Element;
//          typedef Set<T> ColumnClass;
//          typedef typename std::vector< typename U::template Tensor_tpl<U> >::iterator iterator;
//          typedef typename std::vector< typename U::template Tensor_tpl<U> >::const_iterator const_iterator;
//
//        public:
//
//          typename T::Set base_intensity;
//
//        //Protected member variables;
//        protected:
//
//          std::vector<Element> elems;
//
//        //Public member functions
//        public:
//
//          Tensor_tpl<U> ()  { init(); }
//
//          Tensor_tpl<U>(const Set& default_value);
//
//          Tensor_tpl<U>(const Tensor_tpl<U>& t)
//            : base_intensity(t.base_intensity), elems(t.elems) { init(); }
//
//          ~Tensor_tpl<U>() {}
//
//          Tensor_tpl<U>&                              operator= (const Tensor_tpl<U>& t)
//            { base_intensity = t.base_intensity; elems = t.elems; return *this; }
//
//          Element&                                    operator[] (int idx)
//            { return (elems[idx]); }
//
//          const Element&                              operator[] (int idx) const
//            { return (elems[idx]); }
//
//          iterator                                    begin()       { return elems.begin(); }
//          const_iterator                              begin() const { return elems.begin(); }
//
//          iterator                                    end()         { return elems.end(); }
//          const_iterator                              end() const   { return elems.end(); }
//
//
//          Tensor_tpl<U>&                              operator+=  (const Tensor_tpl<U>& tensor)
//            { base_intensity += tensor.base_intensity; for (size_t elem_i = 0; elem_i < size(); elem_i++) operator[](elem_i) += tensor[elem_i]; return *this; }
//
//          Tensor_tpl<U>&                              operator-=  (const Tensor_tpl<U>& tensor)
//              { base_intensity -= tensor.base_intensity; for (size_t elem_i = 0; elem_i < size(); elem_i++) operator[](elem_i) -= tensor[elem_i]; return *this; }
//
//          Tensor_tpl<U>&                              operator*= (const Tensor_tpl<U>& tensor)
//              { base_intensity *= tensor.base_intensity; for (size_t elem_i = 0; elem_i < size(); elem_i++) operator[](elem_i) *= tensor[elem_i]; return *this; }
//
//          Tensor_tpl<U>&                              operator*= (double c)
//              { base_intensity *= c; for (size_t elem_i = 0; elem_i < size(); elem_i++) operator[](elem_i) *= c; return *this; }
//
//          Tensor_tpl<U>&                              operator/= (const Tensor_tpl<U>& tensor)
//              { base_intensity /= tensor.base_intensity; for (size_t elem_i = 0; elem_i < size(); elem_i++) operator[](elem_i) /= tensor[elem_i]; return *this; }
//
//          Tensor_tpl<U>&                              operator/= (double c)
//            { base_intensity /= c; for (size_t elem_i = 0; elem_i < size(); elem_i++) operator[](elem_i) /= c; return *this; }
//
//          Tensor_tpl<U>                               operator+  (const Tensor_tpl<U>& t) const
//            { Tensor_tpl<U> answer (*this); answer += t; return answer; }
//
//          Tensor_tpl<U>                               operator-  (const Tensor_tpl<U>& t) const
//            { Tensor_tpl<U> answer (*this); answer -= t; return answer; }
//
//          Tensor_tpl<U>                               operator* (const Tensor_tpl<U>& t) const
//            { Tensor_tpl<U> answer (*this); answer *= t; return answer; }
//
//          Tensor_tpl<U>                               operator* (double c) const
//            { Tensor_tpl<U> answer (*this); answer *= c; return answer; }
//
//          Tensor_tpl<U>                               operator/ (const Tensor_tpl<U>& t) const
//            { Tensor_tpl<U> answer (*this); answer /= t; return answer; }
//
//          Tensor_tpl<U>                               operator/ (double c) const
//            { Tensor_tpl<U> answer (*this); answer /= c; return answer; }
//
//          Tensor_tpl<U>&                              zero()
//            { base_intensity.zero(); for (iterator it = begin(); it != end(); it++) it->zero(); return *this; }
//
//          Tensor_tpl<U>&                              invalidate()
//            { base_intensity.invalidate(); for (iterator it = begin(); it != end(); it++) it->invalidate(); return *this; }
//
//          size_t                                        size() const
//            { return elems.size(); }
//
//          void                                        push_back(const Element& elem)
//            { elems.push_back(elem); }
//
//          void                                        clear()
//            { base_intensity.clear(); elems.clear(); }
//
//
//          Set<T>                                      get_diagonal() const;
//
//
//          void                                        set_diagonal(const Set<T>& diagonal);
//
//
//        protected:
//
//          void init() {
//
//            if (Fibre::Object::include_intrinsic_intensities)
//              this->intrinsic_props["base_intensity"] = &base_intensity;
//            else
//              base_intensity.set(1.0);
//
//          }
//
//
//      };
//
//	  }
//
//	}
//
//}

#endif /* __bts_fibre_set_tensor_h__ */
