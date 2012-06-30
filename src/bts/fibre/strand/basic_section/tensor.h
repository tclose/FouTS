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
//
#ifndef __bts_fibre_strand_section_tensor_h__
#define __bts_fibre_strand_section_tensor_h__

#include "bts/fibre/strand/basic_section.h"

#include "bts/fibre/base/tensor.h"
#include "bts/coord/tensor.h"

namespace BTS {

  namespace Fibre {

    class Strand::BasicSection::Tensor : public Base::Tensor<Strand::BasicSection> {


      public:

        Tensor()
          : Base::Tensor<Strand::BasicSection>(Strand::BasicSection()) {}

        Tensor(const MR::Math::Matrix<double>::View& view)
          : Base::Tensor<Strand::BasicSection>(Strand::BasicSection(), view) {  assert(view.rows() == 7); assert(view.columns() == 7); }

        Strand::BasicSection 	      position(size_t dim_i)
          { return Strand::BasicSection(row(dim_i), tmpl.props); }

        const Strand::BasicSection 	position(size_t dim_i) const
          { return Strand::BasicSection(row(dim_i), tmpl.props); }

        Strand::BasicSection 	      tangent(size_t dim_i)
          { return Strand::BasicSection(row(dim_i + 3), tmpl.props); }

        const Strand::BasicSection 	tangent(size_t dim_i) const
          { return Strand::BasicSection(row(dim_i + 3), tmpl.props); }

        Strand::BasicSection 	      intensity()
          { return Strand::BasicSection(row(6), tmpl.props); }

        const Strand::BasicSection 	intensity() const
          { return Strand::BasicSection(row(6), tmpl.props); }

        Coord::Tensor 	            position()
          { return Coord::Tensor(sub((size_t)0,(size_t)2,(size_t)0,(size_t)2)); }

        const Coord::Tensor 	      position() const
          { return Coord::Tensor(sub((size_t)0,(size_t)2,(size_t)0,(size_t)2)); }

        Coord::Tensor 	            tangent()
          { return Coord::Tensor(sub((size_t)3,(size_t)5,(size_t)3,(size_t)5)); }

        const Coord::Tensor tangent() const
          { return Coord::Tensor(sub((size_t)3,(size_t)5,(size_t)3,(size_t)5)); }


        Coord::Tensor 	            position_tangent()
          { return Coord::Tensor(sub((size_t)0,(size_t)2,(size_t)3,(size_t)5)); }

        const Coord::Tensor 	      position_tangent() const
          { return Coord::Tensor(sub((size_t)0,(size_t)2,(size_t)3,(size_t)5)); }


        Coord::Tensor 	            tangent_position()
          { return Coord::Tensor(sub((size_t)3,(size_t)5,(size_t)0,(size_t)2)); }

        const Coord::Tensor 	      tangent_position() const
          { return Coord::Tensor(sub((size_t)3,(size_t)5,(size_t)0,(size_t)2)); }

        void                        unnormalise_hessian(const Triple<double>& vox_lengths);
    };

  }

}

#endif /* __bts_fibre_strand_section_tensor_h__ */
