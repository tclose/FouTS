/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

 This file is part of Fourier Tract Sampling (FouTS).

 FouTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FouTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_fibre_strand_basicsection_h__
#define __bts_fibre_strand_basicsection_h__

#include <assert.h>
#include <vector>
#include <string>

#include "bts/triple.h"

#include "bts/fibre/base/object.h"
#include "bts/fibre/strand.h"

namespace FTS {
    
    namespace Fibre {
        
        class Strand::BasicSection: public Base::Object {
                
            public:
                
                class Tensor;
                friend class Tensor;

                typedef Coord Element;
                typedef Fibre::Base::Reader<BasicSection> Reader;
                typedef Fibre::Base::Writer<BasicSection> Writer;

            protected:
                
                const static char* INTENSITY_PROP;
                const static std::vector<const char*> PROPS;

            public:
                const static char* PROPS_LIST[];

                const static std::string FILE_EXTENSION;
                const static Coord FILE_SEPARATOR;

                //Public member functions.
            public:
                
                BasicSection()
                        : Base::Object((size_t) 2, (size_t) 6 + PROPS.size(), PROPS) {
                    Base::Object::set(NAN);
                }
                
                BasicSection(const MR::Math::Vector<double>::View& view,
                             std::vector<const char*>* props)
                        : Base::Object((size_t) 2, view, props) {
                    assert(view.size()==7);
                }
                
                BasicSection(const Strand& strand,
                             const MR::Math::Vector<double>::View& position_coeffs,
                             const MR::Math::Vector<double>::View& tangent_coeffs,
                             double length_fraction)
                        : Base::Object((size_t) 2, (size_t) 3, PROPS)

                {
                    Base::Object::set(NAN);
                    set(strand, position_coeffs, tangent_coeffs, length_fraction);
                }
                
                // Used in derived classes
            protected:
                
                BasicSection(size_t size, const std::vector<const char*>& props)
                        : Base::Object(size, size * 3 + props.size(), props) {
                    set(NAN);
                }
                
                BasicSection(size_t size, const MR::Math::Vector<double>::View& view,
                             std::vector<const char*>* props)
                        : Base::Object(size, view, props) {
                }
                
            public:
                
                BasicSection& set(const Strand& strand,
                                  const MR::Math::Vector<double>::View& position_coeffs,
                                  const MR::Math::Vector<double>::View& tangent_coeffs,
                                  double length_fraction)

                                  {
                    
                    this->intensity() = strand.acs();
                    position() = strand.left_product(position_coeffs);
                    tangent() = strand.left_product(tangent_coeffs) * length_fraction;
                    return *this;
                }
                
                using Base::Object::set;

                virtual ~BasicSection() {
                }
                
                Coord position() {
                    return Coord(MR::Math::Vector<double>::sub(0, 3));
                }
                
                Coord position() const {
                    return Coord(MR::Math::Vector<double>::sub(0, 3));
                }
                
                Coord tangent() {
                    return Coord(MR::Math::Vector<double>::sub(3, 6));
                }
                
                Coord tangent() const {
                    return Coord(MR::Math::Vector<double>::sub(3, 6));
                }
                
                double& position(size_t idx) {
                    return MR::Math::Vector<double>::operator[](idx);
                }
                
                const double& position(size_t idx) const {
                    return MR::Math::Vector<double>::operator[](idx);
                }
                
                double& tangent(size_t idx) {
                    return MR::Math::Vector<double>::operator[](idx + 3);
                }
                
                const double& tangent(size_t idx) const {
                    return MR::Math::Vector<double>::operator[](idx + 3);
                }
                
                double& intensity() {
                    return prop(INTENSITY_PROP);
                }
                
                const double& intensity() const {
                    return prop(INTENSITY_PROP);
                }
                
                BASE_GENERAL_FUNCTIONS(BasicSection)
                ;

                BASE_MULT_DIVIDE_FUNCTIONS(BasicSection)
                ;

                BASE_ADD_SUBTRACT_FUNCTIONS(BasicSection)
                ;

                double length() const {
                    return tangent().norm();
                }
                
                void normalize(const Triple<double>& vox_lengths, const Triple<double>& offsets);

                void unnormalize_gradient(const Triple<double>& vox_lengths);

                MR::Math::Vector<double> vector() const;

                void from_vector(const MR::Math::Vector<double>& vector);

                friend std::ostream& operator<<(std::ostream& stream,
                                                const Strand::BasicSection& seg);
                
                // End generic functions.
        };
        
        inline Strand::BasicSection operator*(double scalar, const Strand::BasicSection& section) {
            return section * scalar;
        }
        
        std::ostream& operator<<(std::ostream& stream, const Strand::BasicSection& seg);
        
        inline void Strand::BasicSection::normalize(const Triple<double>& vox_lengths,
                                                    const Triple<double>& offsets) {
            
            position() -= offsets;
            
            position() /= vox_lengths;
            
        }
        
        inline void Strand::BasicSection::unnormalize_gradient(const Triple<double>& vox_lengths) {
            
            position() /= vox_lengths;
            
        }
    
    }

}

#endif
