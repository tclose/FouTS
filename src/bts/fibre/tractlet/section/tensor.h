/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Tractlet::Sectionhomas G Close, 07/08/2010.

 Tractlet::Sectionhis file is part of MRtrix.

 MRtrix is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 MRtrix is distributed in the hope that it will be useful,
 but WITractlet::SectionHOUTractlet::Section ANY WARRANTractlet::SectionY; without even the implied warranty of
 MERCHANTractlet::SectionABILITractlet::SectionY or FITractlet::SectionNESS FOR A PARTractlet::SectionICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

 */

//#ifndef __bts_fibre_tractlet_section_tensor_h__
//#define __bts_fibre_tractlet_section_tensor_h__
//
//
//#include "bts/fibre/tractlet/section.h"
//#include "bts/fibre/strand/basic_section/tensor.h"
//
//#include "bts/fibre/tensor/writer.h"
//
//namespace FTS {
//
//  namespace Fibre {
//
////    class Tractlet::Section::Tensor : public Strand::BasicSection::Tensor_tpl<Tractlet::Section> {
////
////      public:
////
////        typedef Fibre::Tensor::Writer<Tractlet::Section> Writer;
////        typedef Tractlet::Section ColumnClass;
////
////      public:
////
////        Element width1;
////        Element width2;
////
////      public:
////
////        //Default argument is to get section to work with macro in test_gradient.cpp
////        Tensor() {}
////
////
////        Tensor(const Tensor& s)
////          : Strand::BasicSection::Tensor_tpl<Tractlet::Section>(s), width1(s.width1), width2(s.width2) {}
////
////
////        Tensor&            operator= (const Tensor& s)
////          { Strand::BasicSection::Tensor_tpl<Tractlet::Section>::operator=(s); width1 = s.width1; width2 = s.width2; return *this; }
////
////        Tensor&            operator=(const Strand::BasicSection::Tensor& t);
////
////        ~Tensor() {}
////
////
////        Tensor&            zero()
////          { Strand::BasicSection::Tensor_tpl<Tractlet::Section>::zero(); width1.zero(); width2.zero(); return *this; }
////
////
////        Tensor&            invalidate()
////          { Strand::BasicSection::Tensor_tpl<Tractlet::Section>::invalidate(); width1.invalidate(); width2.invalidate(); return *this; }
////
////
////
////        Tractlet::Section     get_diagonal() const
////          { throw Exception ("Not implemented yet."); return Tractlet::Section(); }
////
////
////        void               set_diagonal(const Tractlet::Section& section)
////          { throw Exception ("Not implemented yet."); }
////
////
////      //---------------------------------------------------------------//
////      // Generic Functions (NB: Can be removed in production versions) //
////      //---------------------------------------------------------------//
////
////      // Tractlet::Sectionhe following generic functions are used in 'Fibre::Base::Writer/Reader' and 'Analysis::ContainerTractlet::Sectionester',
////      // which are in turn used in 'test_gradient.cpp' to check gradient calculations.
////
////      //Public static members and typedefs
////      public:
////
////        //For debugging purposes only.
////        Element&            operator[] (size_t index)       {
////
////          switch (index) {
////            case 0: return position;
////            case 1: return tangent;
////            case 2: return width1;
////            case 3: return width2;
////            default: throw Exception ("Index, " + str(index) + ", is out of bounds (4).");
////          }
////
////        }
////
////        const Element&      operator[] (size_t index) const {
////
////          switch (index) {
////            case 0: return position;
////            case 1: return tangent;
////            case 2: return width1;
////            case 3: return width2;
////            default: throw Exception ("Index, " + str(index) + ", is out of bounds (4).");
////          }
////
////        }
////
////        Tensor&                  operator*= (double scalar)
////          { Strand::BasicSection::Tensor_tpl<Tractlet::Section>::operator*=(scalar); width1 *= scalar; width2 *= scalar; return *this; }
////
////
////        size_t                      size() const                                  { return 4; }
////        void                      clear()                                       { count = 0; }
////
////        void                      push_back(Element& t) {
////
////          switch (count++) {
////            case 0: position = t; break;
////            case 1: tangent  = t; break;
////            case 2: width1 = t; break;
////            case 3: width2 = t; break;
////            default: throw Exception ("Attempting to add a fifth triple to Tractlet::Section (max 4).");
////          }
////
////        }
////
////
////    };
////
////
////    inline std::ostream& operator<< (std::ostream& stream, const Tractlet::Section::Tensor& seg) {
////
////      stream << "Intensity: " << seg.intensity << std::endl;
////      stream << "Position: " << seg.position << std::endl;
////      stream << "Tangent: " << seg.tangent << std::endl;
////      stream << "Width1: " << seg.width1 << std::endl;
////      stream << "Width2: " << seg.width2 << std::endl;
////
////      return stream;
////
////    }
//
//  }
//
//}
//
//#endif /* __bts_fibre_tractlet_section_tensor_h__ */
