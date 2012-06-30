/*
    Copyright 2012 Brain Research Institute, Melbourne, Australia

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


#ifndef __bts_coord_h__
#define __bts_coord_h__

namespace BTS {

  class Coord;


}

#include "math/vector.h"
#include "bts/triple.h"


namespace BTS {

  class Coord : public MR::Math::Vector<double> {

    public:

      class Tensor;
      typedef Utilities::Reader<Coord> Reader;
      typedef Utilities::Writer<Coord> Writer;

    public:

      const static Coord                  Invalid;
      const static Coord                  Zeros;
      const static Coord                  Halves;
      const static Coord                  Ones;

      static const std::string   FILE_EXTENSION;

    public:

      /*! Generates a basic fibre coordinate derived from a base fibre object of length three and row size 1 */
      Coord()
        : MR::Math::Vector<double>(3) { set(NAN); }

      Coord(double x, double y, double z)
        : MR::Math::Vector<double>(3) { set(x,y,z); }

      Coord(const Triple<double>& t)
        : MR::Math::Vector<double>(3) { set(t[0],t[1],t[2]); }

      Coord(const MR::Math::Vector<double>& v)
        : MR::Math::Vector<double>(v) { assert(v.size() == 3); }

      Coord(const MR::Math::Vector<double>::View& view)
        : MR::Math::Vector<double>(view) { assert(view.size() == 3); }

      bool          operator! () const
        { return (isnan(operator[](X)) || isnan(operator[](Y)) || isnan(operator[](Z))); }

      bool          valid () const
        { return (!(isnan(operator[](X)) || isnan(operator[](Y)) || isnan(operator[](Z)))); }


      Coord&        set (double x_value, double y_value, double z_value)
        { operator[](X) = x_value; operator[](Y) = y_value; operator[](Z) = z_value; return *this; }

      Coord&        set (double value)
        { set(value,value,value); return *this; }

      Coord&        zero ();

      Coord&        invalidate ();

      Coord         abs() const
        { return Coord(fabs(operator[](X)), fabs(operator[](Y)), fabs(operator[](Z)));}

      Coord         invert() const
        { return Coord (1/operator[](X), 1/operator[](Y), 1/operator[](Z)); }

      void          copy_to(double* array)                { array[X] = operator[](X); array[Y] = operator[](Y); array[Z] = operator[](Z); }

      bool          operator== (const Coord& A) const;
      bool          operator!= (const Coord& A) const;

      Coord         operator- () const
        { return (Coord (-operator[](X), -operator[](Y), -operator[](Z))); }

      Coord&        operator= (const Coord& A)
        { operator[](X) = A[X]; operator[](Y) = A[Y]; operator[](Z) = A[Z]; return (*this);  }


      //Scalar addition, subtraction, multiplication and division.
      Coord&        operator+= (double inc)
        { operator[](X) += inc; operator[](Y) += inc; operator[](Z) += inc; return *this; }

      Coord&        operator-= (double inc)
          { operator[](X) -= inc; operator[](Y) -= inc; operator[](Z) -= inc; return *this; }

      Coord&        operator*= (double M)
          { operator[](X)*=M; operator[](Y)*=M; operator[](Z)*=M; return *this; }

      Coord&        operator/= (double M)
          { operator[](X)/=M; operator[](Y)/=M; operator[](Z)/=M; return *this; }


      Coord         operator+ (double inc) const
      { return (Coord (operator[](X)+inc, operator[](Y)+inc, operator[](Z)+inc)); }

      Coord         operator- (double inc) const
      { return (Coord (operator[](X)-inc, operator[](Y)-inc, operator[](Z)-inc)); }

      Coord         operator* (double M) const
          { return (Coord (operator[](X)*M, operator[](Y)*M, operator[](Z)*M)); }

      Coord         operator/ (double M) const
          { return (Coord (operator[](X)/M, operator[](Y)/M, operator[](Z)/M)); }


      //Pairwise addition, subtraction, multiplication and division.
      Coord         operator+ (const Coord& c) const
        { return Coord (operator[](X)+c[X], operator[](Y)+c[Y], operator[](Z)+c[Z]); }

      Coord         operator- (const Coord& c) const
        { return Coord (operator[](X)-c[X], operator[](Y)-c[Y], operator[](Z)-c[Z]); }

      Coord         operator* (const Coord& c) const
        { return Coord (operator[](X)*c[X], operator[](Y)*c[Y], operator[](Z)*c[Z]); }

      Coord         operator/ (const Coord& c) const
        { return Coord (operator[](X)/c[X], operator[](Y)/c[Y], operator[](Z)/c[Z]); }


      Coord&        operator+= (const Coord& c)
        { operator[](X) += c[X]; operator[](Y) += c[Y]; operator[](Z) += c[Z]; return *this; }

      Coord&        operator-= (const Coord& c)
        { operator[](X) -= c[X]; operator[](Y) -= c[Y]; operator[](Z) -= c[Z]; return *this; }

      Coord&        operator*= (const Coord& c)
        { operator[](X) *= c[X]; operator[](Y) *= c[Y]; operator[](Z) *= c[Z]; return *this; }

      Coord&        operator/= (const Coord& c)
        { operator[](X) /= c[X]; operator[](Y) /= c[Y]; operator[](Z) /= c[Z]; return *this; }


      bool              positive() const
        { return (operator[](X) >  0.0) && (operator[](Y) >  0.0) && (operator[](Z) >  0.0); }

      bool              non_negative() const
        { return (operator[](X) >= 0.0) && (operator[](Y) >= 0.0) && (operator[](Z) >= 0.0); }

      bool              non_zero() const
        { return operator[](X) || operator[](Y) || operator[](Z); }


      template <typename U> bool  upper_bounded(U bound) const
        { return (operator[](X) <= (double)bound) && (operator[](Y) <= (double)bound) && (operator[](Z) <= (double)bound); }

      template <typename U> bool  lower_bounded(U bound) const
        { return (operator[](X) >= (double)bound) && (operator[](Y) >= (double)bound) && (operator[](Z) >= (double)bound); }

      template <typename U> bool  upper_bounded(const Triple<U>& bound) const
        { return (operator[](X) <= (double)bound[X]) && (operator[](Y) <= (double)bound[Y]) && (operator[](Z) <= (double)bound[Z]); }

      template <typename U> bool  lower_bounded(const Triple<U>& bound) const
        { return (operator[](X) >= (double)bound[X]) && (operator[](Y) >= (double)bound[Y]) && (operator[](Z) >= (double)bound[Z]); }

      double                 dot (const Coord& A) const
        { return (operator[](X)*A[X] + operator[](Y)*A[Y] + operator[](Z)*A[Z]); }

      double                 dot (const Coord& t, Coord& gradient) const;

      double                 angle(const Coord& t) const;

      Tensor                 outer(const Coord& t) const;

      double                 norm2 () const
        { return (operator[](X)*operator[](X) + operator[](Y)*operator[](Y) + operator[](Z)*operator[](Z)); }

      double                 norm () const
        { return (sqrt (norm2())); }

      double                 norm2 (Coord& gradient) const;

      double                 norm (Coord& gradient) const;

      double                 norm2 (Coord& gradient, Tensor& hessian) const;

      double                 norm (Coord& gradient, Tensor& hessian) const;

      const Coord&           normalise ()
        { double nrm = norm(); operator[](X) /= nrm; operator[](Y) /= nrm; operator[](Z) /= nrm; return (*this);  }

      Coord                  cross (const Coord& A) const
        { return (Coord (operator[](Y)*A[Z]-operator[](Z)*A[Y], operator[](Z)*A[X]-operator[](X)*A[Z], operator[](X)*A[Y]-operator[](Y)*A[X]));  }

      Coord                  min_axis() const
        { Coord axis(0.0,0.0,0.0); if (fabs(operator[](X)) < fabs(operator[](Y))) { if (fabs(operator[](X)) < fabs(operator[](Z))) axis[X] = 1.0; else axis[Z] = 1.0; } else if (fabs(operator[](Y)) < fabs(operator[](Z))) axis[Y] = 1.0; else axis[Z] = 1.0; return axis; }

      Coord         	       sign() const
        { Coord s; for (size_t i = 0; i < 3; ++i) s[i] = operator[](i) > 0 ? 1 : (operator[](i) < 0 ? -1 : 0); return s; }

      double                 vsize() const
        { return 3; }

      void                   clear()
      { // do nothing, overrides MR::Math::Vector<double> clear to make sure it doesn't resize the vector.
        // TODO: Reimplement Coord class so that it doesn't need to use MR::Math::Vector<double> and remove redundencies like this one.
      }


  };

  std::ostream&  operator<< (std::ostream& stream , const Coord& c);


  inline Coord  operator* (double M, Coord c) { return Coord(M * c[X], M * c[Y], M * c[Z]); }
  inline Coord  operator+ (double inc, Coord c) { return Coord(inc + c[X], inc + c[Y], inc + c[Z]); }

}

#endif /* __fibre_coord_h__ */
