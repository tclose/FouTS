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



#ifndef __bts_fibre_strand_section_h__
#define __bts_fibre_strand_section_h__

#include "bts/fibre/strand.h"

#include <assert.h>
#include <vector>
#include <string>

#include "bts/fibre/base/writer.h"
#include "bts/fibre/base/reader.h"

//#include "bts/image/observed/buffer.h"



#include "bts/fibre/strand/basic_section.h"
#include "bts/fibre/strand/basic_section/tensor.h"
#include "bts/fibre/base/tensor.h"
#include "bts/coord/tensor.h"


namespace BTS {

  namespace Fibre {

    class Strand::Section : public Strand::BasicSection {

      public:

        typedef Coord Element;

        typedef Fibre::Base::Reader<Section> Reader;
        typedef Fibre::Base::Writer<Section> Writer;

     //Public static members and typedefs
     public:

       const static std::string     FILE_EXTENSION;
       const static Coord           FILE_SEPARATOR;


      // Public member variables.
      public:

        //References used when working backwards to create the gradient.
        MR::Math::Vector<double> position_coeffs;
        MR::Math::Vector<double> tangent_coeffs;
        double length_fraction;

        std::vector<double> precalc_weightings;
        std::vector<Coord> precalc_weight_gradients;
        std::vector<Coord::Tensor> precalc_weight_hessians;

        double precalc_interpolation;
        BasicSection precalc_interp_gradient;
        BasicSection::Tensor precalc_interp_hessian;

        const Strand*  parent;

      public:

        Section(size_t init_precalc_size = 0)
         : length_fraction(NAN),
           precalc_weightings(init_precalc_size),
           precalc_weight_gradients(init_precalc_size),
           precalc_weight_hessians(init_precalc_size) {}


        //TODO: Look into this, this is not the most efficient way to do this but will have to do until image generation
        // is reworked.
        Section(const Strand& strand,
                 const MR::Math::Vector<double>::View& position_coeffs,
                 const MR::Math::Vector<double>::View& tangent_coeffs,
                 double length_fraction,
                 size_t init_precalc_size = 0)

           : BasicSection(strand,
                          position_coeffs,
                          tangent_coeffs,
                          length_fraction),
             precalc_weightings(init_precalc_size),
             precalc_weight_gradients(init_precalc_size),
             precalc_weight_hessians(init_precalc_size)

         { set(strand,position_coeffs,tangent_coeffs,length_fraction); }


        void set(const Strand& strand,
               const MR::Math::Vector<double>::View& position_coeffs,
               const MR::Math::Vector<double>::View& tangent_coeffs,
               double length_fraction)

        {

          BasicSection::set(strand, position_coeffs, tangent_coeffs, length_fraction);

          this->position_coeffs = position_coeffs;
          this->tangent_coeffs = tangent_coeffs;
          this->length_fraction = length_fraction;

          this->parent = &strand;
        }


        Section(const Section& s)
         : BasicSection(s),
           position_coeffs(s.position_coeffs),
           tangent_coeffs(s.tangent_coeffs),
           length_fraction(s.length_fraction),
           precalc_weightings(s.precalc_weightings),
           precalc_weight_gradients(s.precalc_weight_gradients),
           precalc_weight_hessians(s.precalc_weight_hessians),
           precalc_interpolation(s.precalc_interpolation),
           precalc_interp_gradient(s.precalc_interp_gradient),
           precalc_interp_hessian(s.precalc_interp_hessian),
           parent(s.parent) {}


      // Used in derived classes
      protected:

        Section(size_t size, const std::vector<const char*>& props, size_t init_precalc_size)
          : BasicSection(size, props),
            length_fraction(NAN),
             precalc_weightings(init_precalc_size),
             precalc_weight_gradients(init_precalc_size),
             precalc_weight_hessians(init_precalc_size) {}

        Section(size_t size, const MR::Math::Vector<double>::View& view, std::vector<const char*>* props)
          : BasicSection(size, view, props) {}

      public:

        Section&                  operator= (const Section& s) {
          BasicSection::operator=(s);
          length_fraction = s.length_fraction;
          position_coeffs = s.position_coeffs;
          tangent_coeffs = s.tangent_coeffs;
          precalc_weightings = s.precalc_weightings;
          precalc_weight_gradients = s.precalc_weight_gradients;
          precalc_weight_hessians = s.precalc_weight_hessians;
          precalc_interpolation = s.precalc_interpolation;
          precalc_interp_gradient = s.precalc_interp_gradient;
          precalc_interp_hessian = s.precalc_interp_hessian;
          parent = s.parent;
          return *this;
        }


        virtual ~Section() {}

      protected:

         Section(double intensity,
                 const MR::Math::Vector<double>::View& position_coeffs,
                 const MR::Math::Vector<double>::View& tangent_coeffs,
                 double length_fraction) {

           set(intensity, position_coeffs, tangent_coeffs, length_fraction);
         }


         void set(  double intensity,
                   const MR::Math::Vector<double>::View& position_coeffs,
                   const MR::Math::Vector<double>::View& tangent_coeffs,
                   double length_fraction) {

           this->intensity() = intensity;
           this->position_coeffs = position_coeffs;
           this->tangent_coeffs = tangent_coeffs;
           this->length_fraction = length_fraction;
           this->parent = 0;
         }


      public:

         BASE_GENERAL_FUNCTIONS(Section);

         BASE_MULT_DIVIDE_FUNCTIONS(Section);

         BASE_ADD_SUBTRACT_FUNCTIONS(Section);

      //---------------------------------------------------------------//
      // Generic Functions (NB: Can be removed in production versions) //
      //---------------------------------------------------------------//

      // The following generic functions are used in 'Fibre::Base::Writer/Reader' and 'Analysis::GradientTester',
      // which are in turn used in 'test_gradient.cpp' to check gradient calculations. Can be removed in production version.


//        typedef Fibre::Base::Reader<Strand::Section> Reader;
//        typedef Fibre::Base::Writer<Strand::Section> Writer;


    };


  }
  
}






#endif
