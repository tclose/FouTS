/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

 This file is part of Bayesian Tractlet Sampling (BTS).

 BTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 BTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_fibre_tractlet_section_h__
#define __bts_fibre_tractlet_section_h__

#include "bts/fibre/tractlet.h"
#include "bts/fibre/strand/section.h"

namespace BTS {
    
    namespace Fibre {
        
        class Tractlet::Section: public Strand::Section {
                
            public:
                
                //TODO: Work out how to do this more elegantly
//        class Tensor;
                
                typedef Fibre::Base::Reader<Section> Reader;
                typedef Fibre::Base::Writer<Section> Writer;

            public:
                
                const static std::string FILE_EXTENSION;
                const static Coord FILE_SEPARATOR;

            public:
                
                //Reference data. Used in working backwards to calculate gradient.
                double ax1_fraction;    //The fraction of the section along axis1 (from -1 to 1)
                double ax2_fraction;    //The fraction of the section along axis2 (from -1 to 1)
                double width_fraction;

                //TODO: Remove this reference.
                const Tractlet* parent;

            public:
                
                Section(const MR::Math::Vector<double>::View& view, std::vector<const char*>* props)
                        : Strand::Section((size_t) 4, view, props) {
                }
                
                Section(size_t init_precalc_size = 0)
                        : Strand::Section((size_t) 4, PROPS, init_precalc_size), ax1_fraction(NAN), ax2_fraction(
                                  NAN), width_fraction(NAN) {
                }
                
                Section(const Tractlet& tractlet,
                        const MR::Math::Vector<double>::View& position_coeffs,
                        const MR::Math::Vector<double>::View& tangent_coeffs, double ax1_fraction,
                        double ax2_fraction, double length_fraction, double width_fraction,
                        double intensity_scale, size_t init_precalc_size = 0)
                        
                        : Strand::Section(init_precalc_size)

                {
                    
                    set(tractlet, position_coeffs, tangent_coeffs, ax1_fraction, ax2_fraction,
                            length_fraction, width_fraction, intensity_scale);
                    
                }
                
                void set(const Tractlet& tractlet,
                         const MR::Math::Vector<double>::View& position_coeffs,
                         const MR::Math::Vector<double>::View& tangent_coeffs, double ax1_fraction,
                         double ax2_fraction, double length_fraction, double width_fraction,
                         double intensity_scale) {
                    
                    Strand::Section::set(tractlet.acs(), position_coeffs, tangent_coeffs,
                            length_fraction);
                    this->intensity() *= intensity_scale;
                    
                    this->ax1_fraction = ax1_fraction;
                    this->ax2_fraction = ax2_fraction;
                    this->width_fraction = width_fraction;
                    
                    width1() = tractlet[1].left_product(position_coeffs);
                    width2() = tractlet[2].left_product(position_coeffs);
                    
                    position() = tractlet[0].left_product(position_coeffs) + ax1_fraction * width1()
                                 + ax2_fraction * width2();
                    
                    width1() *= width_fraction;
                    width2() *= width_fraction;
                    
                    tangent() = (tractlet[0].left_product(tangent_coeffs)
                            + ax1_fraction * tractlet[1].left_product(tangent_coeffs)
                            + ax2_fraction * tractlet[2].left_product(tangent_coeffs))
                                * length_fraction;
                    
                    this->parent = &tractlet;
                    
                }
                
                Section(const Section& s)
                        : Strand::Section(s), ax1_fraction(s.ax1_fraction), ax2_fraction(
                                  s.ax2_fraction), parent(s.parent) {
                }
                
                Section& operator=(const Section& s) {
                    this->Strand::Section::operator=(s);
                    ax1_fraction = s.ax1_fraction;
                    ax2_fraction = s.ax2_fraction;
                    parent = s.parent;
                    return *this;
                }
                
                ~Section() {
                }
                
                Coord width1() {
                    return Coord(MR::Math::Vector<double>::sub(6, 9));
                }
                
                Coord width1() const {
                    return Coord(MR::Math::Vector<double>::sub(6, 9));
                }
                
                Coord width2() {
                    return Coord(MR::Math::Vector<double>::sub(9, 12));
                }
                
                Coord width2() const {
                    return Coord(MR::Math::Vector<double>::sub(9, 12));
                }
                
                BASE_GENERAL_FUNCTIONS(Section)
                ;

                BASE_MULT_DIVIDE_FUNCTIONS(Section)
                ;

                BASE_ADD_SUBTRACT_FUNCTIONS(Section)
                ;
                
        };
        
        inline Tractlet::Section operator*(double scalar, const Tractlet::Section& section) {
            return section * scalar;
        }
        
        std::ostream& operator<<(std::ostream& stream,
                                 const BTS::Fibre::Tractlet::Section& section);
    
    }
}

#endif
