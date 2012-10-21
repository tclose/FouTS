/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 07/08/2010.

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

#ifndef __bts_fibre_strand_tensor_h__
#define __bts_fibre_strand_tensor_h__


#include "bts/fibre/strand.h"
#include "bts/fibre/strand/tensor.h"
#include "bts/fibre/base/tensor.h"
#include "bts/fibre/strand/basic_section/tensor.h"



namespace BTS {

  namespace Fibre {

    class Strand::Tensor : public Base::Tensor<Strand> {


      public:

//        typedef Base<Strand>::Tensor::Writer Writer;

      //Public member functions
      public:

        Tensor() {}

        Tensor(const Strand& strand_template)
          : Base::Tensor<Strand>(strand_template) {}


        Tensor (const Strand& strand_template, const MR::Math::Matrix<double>::View& view)
          : Base::Tensor<Strand>(strand_template, view) {}


        Tensor(const Tensor& t)
          : Base::Tensor<Strand>(t) {}


        Tensor& operator=(const Tensor& t)
          { Base::Tensor<Strand>::operator=(t); return *this; }

        Tensor& operator=(const MR::Math::Matrix<double>& m)
          { Base::Tensor<Strand>::operator=(m); return *this; }

        Tensor& operator=(double c)
          { Base::Tensor<Strand>::operator=(c); return *this; }

        size_t        degree() const
          { return tmpl.degree(); }

        size_t        vsize() const
          { return tmpl.vsize(); }

        Strand        row(size_t idx)
          { return Strand(degree(), MR::Math::Matrix<double>::row(idx), tmpl.props, 0); }

        const Strand  row(size_t idx) const
          { return Strand(degree(), MR::Math::Matrix<double>::row(idx), tmpl.props, 0); }

        Strand        column(size_t idx)
          { return Strand(degree(), MR::Math::Matrix<double>::column(idx), tmpl.props, 0); }

        const Strand  column(size_t idx) const
          { return Strand(degree(), MR::Math::Matrix<double>::column(idx), tmpl.props, 0); }

        bool          has_var_acs() const
        { return tmpl.has_var_acs(); }

        Strand        acs()
          { assert(has_var_acs()); return Strand(degree(), MR::Math::Matrix<double>::row(vsize() + tmpl.prop_index(Strand::ALPHA_PROP)), tmpl.props, 0); }


        double&       operator()(size_t row_degree_i, size_t row_dim_i, size_t col_degree_i, size_t col_dim_i)
         { return MR::Math::Matrix<double>::operator()(row_degree_i * 3 + row_dim_i, col_degree_i * 3 + col_dim_i); }


        double        operator()(size_t row_degree_i, size_t row_dim_i, size_t col_degree_i, size_t col_dim_i) const
         { return MR::Math::Matrix<double>::operator()(row_degree_i * 3 + row_dim_i, col_degree_i * 3 + col_dim_i); }


        Strand        operator[](size_t idx)
          { return Strand(degree(), MR::Math::Matrix<double>::row(idx), tmpl.props, 0); }

        const Strand  operator[](size_t idx) const
          { return Strand(degree(), MR::Math::Matrix<double>::row(idx), tmpl.props, 0); }

        void          add_section_hessian(const Strand& strand, const Section& section, const BasicSection& gradient, const BasicSection::Tensor& hessian);


    };

  }

}

#endif /* __bts_fibre_strand_tensor_h__ */
