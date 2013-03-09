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

#ifndef __bts_fibre_tractlet_tensor_h__
#define __bts_fibre_tractlet_tensor_h__

#include "bts/fibre/tractlet.h"
#include "bts/fibre/base/tensor.h"
#include "bts/fibre/strand/basic_section/tensor.h"

namespace BTS {
    
    namespace Fibre {
        
        class Tractlet::Tensor: public Base::Tensor<Tractlet> {
                
            public:
                
//        typedef Base::TensorWriter<Tractlet::Tensor> Writer;
                
//Public member functions
            public:
                
                Tensor() {
                }
                
                Tensor(const Tractlet& tractlet_template)
                        : Base::Tensor<Tractlet>(tractlet_template) {
                }
                
                Tensor(const Tractlet& tractlet_template,
                       const MR::Math::Matrix<double>::View& view)
                        : Base::Tensor<Tractlet>(tractlet_template, view) {
                }
                
                Tensor(const Tensor& t)
                        : Base::Tensor<Tractlet>(t) {
                }
                
                Tensor& operator=(const Tensor& t) {
                    Base::Tensor<Tractlet>::operator=(t);
                    return *this;
                }
                
                Tensor& operator=(const MR::Math::Matrix<double>& m) {
                    Base::Tensor<Tractlet>::operator=(m);
                    return *this;
                }
                
                Tensor& operator=(double c) {
                    Base::Tensor<Tractlet>::operator=(c);
                    return *this;
                }
                
                size_t degree() const {
                    return tmpl.degree();
                }
                
                size_t vsize() const {
                    return tmpl.vsize();
                }
                
                Tractlet row(size_t idx) {
                    return Tractlet(degree(), MR::Math::Matrix<double>::row(idx), tmpl.props, 0);
                }
                
                const Tractlet row(size_t idx) const {
                    return Tractlet(degree(), MR::Math::Matrix<double>::row(idx), tmpl.props, 0);
                }
                
                Tractlet row(size_t ax_i, size_t degree_i, size_t dim_i) {
                    return Tractlet(degree(),
                            MR::Math::Matrix<double>::row((ax_i * degree() + degree_i) * 3 + dim_i),
                            tmpl.props, 0);
                }
                
                const Tractlet row(size_t ax_i, size_t degree_i, size_t dim_i) const {
                    return Tractlet(degree(),
                            MR::Math::Matrix<double>::row((ax_i * degree() + degree_i) * 3 + dim_i),
                            tmpl.props, 0);
                }
                
                Tractlet column(size_t idx) {
                    return Tractlet(degree(), MR::Math::Matrix<double>::column(idx), tmpl.props, 0);
                }
                
                const Tractlet column(size_t idx) const {
                    return Tractlet(degree(), MR::Math::Matrix<double>::column(idx), tmpl.props, 0);
                }
                
                bool has_var_acs() const {
                    return tmpl.has_var_acs();
                }
                
                Tractlet acs() {
                    assert(has_var_acs());
                    return Tractlet(degree(),
                            MR::Math::Matrix<double>::row(
                                    vsize() + tmpl.prop_index(Object::ALPHA_PROP)), tmpl.props, 0);
                }
                
                double& operator()(size_t row_ax_i, size_t row_degree_i, size_t row_dim_i,
                                   size_t col_ax_i, size_t col_degree_i, size_t col_dim_i) {
                    return MR::Math::Matrix<double>::operator()(
                            row_ax_i * tmpl.degree() * 3 + row_degree_i * 3 + row_dim_i,
                            col_ax_i * tmpl.degree() * 3 + col_degree_i * 3 + col_dim_i);
                }
                
                double operator()(size_t row_ax_i, size_t row_degree_i, size_t row_dim_i,
                                  size_t col_ax_i, size_t col_degree_i, size_t col_dim_i) const {
                    return MR::Math::Matrix<double>::operator()(
                            row_ax_i * tmpl.degree() * 3 + row_degree_i * 3 + row_dim_i,
                            col_ax_i * tmpl.degree() * 3 + col_degree_i * 3 + col_dim_i);
                }
                
                Tractlet operator[](size_t idx) {
                    return Tractlet(degree(), MR::Math::Matrix<double>::row(idx), tmpl.props, 0);
                }
                
                const Tractlet operator[](size_t idx) const {
                    return Tractlet(degree(), MR::Math::Matrix<double>::row(idx), tmpl.props, 0);
                }
                
                void add_section_hessian(const Tractlet& tractlet, const Tractlet::Section& section,
                                         const Strand::BasicSection& gradient,
                                         const Strand::BasicSection::Tensor& hessian);
                
        };
    
    }

}

#include "bts/fibre/base/tensor_writer.h"

#endif /* __bts_fibre_tractlet_tensor_h__ */
