/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Aug 16, 2010.

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

#include "bts/fibre/strand/basic_section/tensor.h"

namespace BTS {

  namespace Fibre {

    void    Strand::BasicSection::Tensor::unnormalise_hessian(const Triple<double>& vox_lengths) {

      Triple<double> vox_lengths2 (MR::Math::pow2(vox_lengths[X]), MR::Math::pow2(vox_lengths[Y]),MR::Math::pow2(vox_lengths[Z]));

      position(X).position() /= vox_lengths2;
      position(Y).position() /= vox_lengths2;
      position(Z).position() /= vox_lengths2;

      position(X).tangent() /= vox_lengths[X];
      position(Y).tangent() /= vox_lengths[Y];
      position(Z).tangent() /= vox_lengths[Z];

      tangent(X).position() /= vox_lengths;
      tangent(Y).position() /= vox_lengths;
      tangent(Z).position() /= vox_lengths;

    }

  }

}
