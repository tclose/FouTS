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
//
//#ifndef __bts_fibre_strand_section_tensor_cpp_h__
//#define __bts_fibre_strand_section_tensor_cpp_h__
//
//namespace FTS {
//
//	namespace Fibre {
//
//
//    template <typename T> std::ostream& operator<< (std::ostream& stream , const Strand::BasicSection::Tensor_tpl<T>& seg) {
//
//      stream << "Intensity: " << seg.intensity << std::endl;
//      stream << "Position: " << seg.position << std::endl;
//      stream << "Tangent: " << seg.tangent << std::endl;
//
//      return stream;
//
//    }
//
//
//
//    inline Strand::BasicSection::Tensor::Tensor (const Strand::Section& section) {
//
//      position[X] = section;
//      position[Y] = section;
//      position[Z] = section;
//
//      tangent[X] = section;
//      tangent[Y] = section;
//      tangent[Z] = section;
//
//      intensity = section;
//
//      init();
//
//    }
//
//
//
//    inline Coord::Tensor    Strand::BasicSection::Tensor::position_hessian() const {
//
//      Coord::Tensor pos_hess;
//
//      pos_hess[X] = position[X].position;
//      pos_hess[Y] = position[Y].position;
//      pos_hess[Z] = position[Z].position;
//
//      return pos_hess;
//
//    }
//
//    inline Coord::Tensor    Strand::BasicSection::Tensor::tangent_hessian() const {
//
//      Coord::Tensor tan_hess;
//
//      tan_hess[X] = tangent[X].tangent;
//      tan_hess[Y] = tangent[Y].tangent;
//      tan_hess[Z] = tangent[Z].tangent;
//
//      return tan_hess;
//
//    }
//
//    inline Coord::Tensor    Strand::BasicSection::Tensor::cross_hessian() const {
//
//      Coord::Tensor cross_hess;
//
//      cross_hess[X] = position[X].tangent;
//      cross_hess[Y] = position[Y].tangent;
//      cross_hess[Z] = position[Z].tangent;
//
//      return cross_hess;
//
//    }
//
//
//    inline void    Strand::BasicSection::Tensor::unnormalise_hessian(const Triple<double>& vox_lengths) {
//
//      Triple<double> vox_lengths2 (MR::Math::pow2(vox_lengths[X]), MR::Math::pow2(vox_lengths[Y]),MR::Math::pow2(vox_lengths[Z]));
//
//      position[X].position /= vox_lengths2;
//      position[Y].position /= vox_lengths2;
//      position[Z].position /= vox_lengths2;
//
//      position[X].tangent /= vox_lengths[X];
//      position[Y].tangent /= vox_lengths[Y];
//      position[Z].tangent /= vox_lengths[Z];
//
//      tangent[X].position /= vox_lengths;
//      tangent[Y].position /= vox_lengths;
//      tangent[Z].position /= vox_lengths;
//
//    }
//
//	}
//
//}
//
//#endif /* __bts_fibre_strand_section_tensor_cpp_h__ */
