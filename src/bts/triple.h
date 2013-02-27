/*
    Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

    Written by J-Donald Tournier, 27/06/08.

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


#include <map>

#include "bts/common.h"
#include "bts/coord.h"

#ifndef __bts_triple_h__
#define __bts_triple_h__

#include <gsl/gsl_randist.h>

#include "mrtrix.h"
#include "math/math.h"
#include "math/vector.h"
#include "bts/utilities/reader.h"
#include "bts/utilities/writer.h"


namespace BTS {


  template <typename T> class Triple {

    public:

      typedef Utilities::Reader<Triple<T> > Reader;
      typedef Utilities::Writer<Triple<T> > Writer;

      template <typename U> class Tensor_tpl; //Dummy class used to get template functions in Hessian Tester to work;
      typedef T Element;

      const static Triple<T>                  Invalid;
      const static Triple<T>                  Zeros;
      const static Triple<T>                  Halves;
      const static Triple<T>                  Ones;

      const static std::string                FILE_EXTENSION;


    protected:
      T p[3];

    public:

      Triple<T> ()
        { invalidate(); }

      Triple<T> (T value)
        { set(value); }

      Triple<T> (T x_value, T y_value, T z_value)
        { set (x_value, y_value, z_value); }

      Triple<T> (const Coord& coord);

      bool              operator! () const                    { return (isnan(p[X]) || isnan(p[Y]) || isnan(p[Z])); }
      bool              valid () const                        { return (!(isnan(p[X]) || isnan(p[Y]) || isnan(p[Z]))); }

      T*                get ()                                { return (p);  }
      const T*          get () const                          { return (p); }
      void              copy_to(double* array)                { array[X] = (double)p[X]; array[Y] = (double)p[Y]; array[Z] = (double)p[Z]; }
      void              copy_to(size_t* array)                { array[X] = (size_t)p[X]; array[Y] = (size_t)p[Y]; array[Z] = (size_t)p[Z]; }

      Triple<T>&        set (T x_value, T y_value, T z_value) { p[X] = x_value; p[Y] = y_value; p[Z] = z_value; return *this; }
      Triple<T>&        set (const T point[3])                { set(point[X], point[Y], point[Z]); return *this; }
      Triple<T>&        set (T value)                         { set(value,value,value); return *this; }

      Triple<T>&        zero ();
      Triple<T>&        invalidate ();

      Triple<T>         abs() const                           { return Triple<T>(fabs(p[X]), fabs(p[Y]), fabs(p[Z]));}

      Triple<T>         invert() const                        { return Triple<T> (1/p[X], 1/p[Y], 1/p[Z]); }

      T&                operator[] (int idx)                  { return (p[idx]);  }
      const T&          operator[] (int idx) const            { return (p[idx]); }

      bool              operator== (const Triple<T>& A) const;
      bool              operator!= (const Triple<T>& A) const;

      Triple<T>         operator- () const                    { return (Triple<T> (-p[X], -p[Y], -p[Z])); }

      Triple<T>&        operator= (const Triple<T>& A)        { p[X] = A[X]; p[Y] = A[Y]; p[Z] = A[Z]; return (*this);  }


      MR::Math::Vector<double> vector() const
        { MR::Math::Vector<double> vec(3); vec[X] = p[X]; vec[Y] = p[Y]; vec[Z] = p[Z]; return vec; }


      void              from_vector(const MR::Math::Vector<double>& vec)
        { assert(vec.size() == 3); p[X] = vec[X]; p[Y] = vec[Y]; p[Z] = vec[Z]; }

      //Scalar addition, subtraction, multiplication and division.
      Triple<T>&        operator+= (T inc)                    { p[X] += inc; p[Y] += inc; p[Z] += inc; return *this; }
      Triple<T>&        operator-= (T inc)                    { p[X] -= inc; p[Y] -= inc; p[Z] -= inc; return *this; }
      Triple<T>&        operator*= (T M)                      { p[X]*=M; p[Y]*=M; p[Z]*=M; return *this; }
      Triple<T>&        operator/= (T M)                      { p[X]/=M; p[Y]/=M; p[Z]/=M; return *this; }


      Triple<T>         operator+ (T M) const                 { return (Triple<T> (p[X]+M, p[Y]+M, p[Z]+M)); }
      Triple<T>         operator- (T M) const                 { return (Triple<T> (p[X]-M, p[Y]-M, p[Z]-M)); }
      Triple<T>         operator* (T M) const                 { return (Triple<T> (p[X]*M, p[Y]*M, p[Z]*M)); }
      Triple<T>         operator/ (T M) const                 { return (Triple<T> (p[X]/M, p[Y]/M, p[Z]/M)); }

      //Pairwise addition, subtraction, multiplication and division.
      Triple<T>         operator+ (const Triple<T>& A) const  { return Triple<T> (p[X]+A[X], p[Y]+A[Y], p[Z]+A[Z]); }
      Triple<T>         operator- (const Triple<T>& A) const  { return Triple<T> (p[X]-A[X], p[Y]-A[Y], p[Z]-A[Z]); }
      Triple<T>         operator* (const Triple<T>& A) const  { return Triple<T> (p[X]*A[X], p[Y]*A[Y], p[Z]*A[Z]); }
      Triple<T>         operator/ (const Triple<T>& A) const  { return Triple<T> (p[X]/A[X], p[Y]/A[Y], p[Z]/A[Z]); }

      Triple<T>&        operator+= (const Triple<T>& A)       { p[X] += A[X]; p[Y] += A[Y]; p[Z] += A[Z]; return *this; }
      Triple<T>&        operator-= (const Triple<T>& A)       { p[X] -= A[X]; p[Y] -= A[Y]; p[Z] -= A[Z]; return *this; }
      Triple<T>&        operator*= (const Triple<T>& A)       { p[X] *= A[X]; p[Y] *= A[Y]; p[Z] *= A[Z]; return *this; }
      Triple<T>&        operator/= (const Triple<T>& A)       { p[X] /= A[X]; p[Y] /= A[Y]; p[Z] /= A[Z]; return *this; }


//      //Scalar addition, subtraction, multiplication and division by doubles.
//      template <typename U> Triple<T>&        operator+= (U inc);
//      template <typename U> Triple<T>&        operator-= (U inc);
//      template <typename U> Triple<T>&        operator*= (U M);
//      template <typename U> Triple<T>&        operator/= (U M);
//
//
//      template <typename U> Triple<T>         operator+ (U M) const;
//      template <typename U> Triple<T>         operator- (U M) const;
//      template <typename U> Triple<T>         operator* (U M) const;
//      template <typename U> Triple<T>         operator/ (U M) const;


      bool              positive() const                      { return (p[X] >  0.0) && (p[Y] >  0.0) && (p[Z] >  0.0); }
      bool              non_negative() const                  { return (p[X] >= 0.0) && (p[Y] >= 0.0) && (p[Z] >= 0.0); }
      bool              non_zero() const                      { return p[X] || p[Y] || p[Z]; }

      template <typename U> bool  upper_bounded(U bound) const          { return (p[X] <= (T)bound) && (p[Y] <= (T)bound) && (p[Z] <= (T)bound); }
      template <typename U> bool  lower_bounded(U bound) const          { return (p[X] >= (T)bound) && (p[Y] >= (T)bound) && (p[Z] >= (T)bound); }

      template <typename U> bool  upper_bounded(const Triple<U>& bound) const          { return (p[X] <= (T)bound[X]) && (p[Y] <= (T)bound[Y]) && (p[Z] <= (T)bound[Z]); }
      template <typename U> bool  lower_bounded(const Triple<U>& bound) const          { return (p[X] >= (T)bound[X]) && (p[Y] >= (T)bound[Y]) && (p[Z] >= (T)bound[Z]); }

      T                 dot (const Triple<T>& A) const        { return (p[X]*A[X] + p[Y]*A[Y] + p[Z]*A[Z]); }
      T                 dot (const Triple<T>& t, Triple<T>& gradient) const;

      T                 angle(const Triple<T>& t) const;

      T                 norm2 () const                        { return (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]); }
      T                 norm () const                         { return (sqrt (norm2())); }


      const Triple<T>&  normalise ()                          { MR::Math::normalise (p); return (*this);  }


      Triple<T>         cross (const Triple<T>& A) const      { return (Triple<T> (p[Y]*A[Z]-p[Z]*A[Y], p[Z]*A[X]-p[X]*A[Z], p[X]*A[Y]-p[Y]*A[X]));  }

      Triple<T>         min_axis() const                      { Triple<T> axis(0.0,0.0,0.0); if (fabs(p[X]) < fabs(p[Y])) { if (fabs(p[X]) < fabs(p[Z])) axis[X] = 1.0; else axis[Z] = 1.0; } else if (fabs(p[Y]) < fabs(p[Z])) axis[Y] = 1.0; else axis[Z] = 1.0; return axis; }


      Triple<T>         sign() const                          { Triple<T> s; for (size_t i = 0; i < 3; ++i) s[i] = p[i] > 0 ? 1 : (p[i] < 0 ? -1 : 0); return s; }


      //These functions and static members are only required for Fibre::Base::Reader/Writer< Triple<double> > functions.
      void              clear()                               {}

      size_t              size() const                        { return 3; }

      void              resize(size_t dummy)                  { assert(dummy == 3); }

    protected:

      explicit Triple<T> (const T point[3])                   { set (point);  }

  };

#include "bts/coord.h"
#include "bts/common.h"


  template <> inline Triple<int>&                       Triple<int>::zero() {
    memset (p, 0, sizeof(p)); return *this;
  }

  template <> inline Triple<double>&                    Triple<double>::zero() {
    memset (p, 0, sizeof(p)); return *this;
  }

  template <> inline Triple<float>&                     Triple<float>::zero() {
    memset (p, 0, sizeof(p)); return *this;
  }

  template <> inline Triple<size_t>&                    Triple<size_t>::zero() {
    memset (p, 0, sizeof(p)); return *this;
  }


  template <> inline Triple<int>&                       Triple<int>::invalidate() {
    p[X] = p[Y] = p[Z] = std::numeric_limits<int>::min(); return *this;
  }


  template <> inline Triple<size_t>&                    Triple<size_t>::invalidate() {
    p[X] = p[Y] = p[Z] = std::numeric_limits<size_t>::max(); return *this;
  }


  template <> inline Triple<double>&                    Triple<double>::invalidate() {
    p[X] = p[Y] = p[Z] = NAN; return *this;
  }

  template <> inline Triple<float>&                     Triple<float>::invalidate() {
    p[X] = p[Y] = p[Z] = NAN; return *this;
  }

  template<typename T> std::ostream&  operator<< (std::ostream& stream , const Triple<T>& pt)
  {
    stream << "[ " << pt[X] << ", " << pt[Y] << ", " << pt[Z] << " ]";
    return (stream);
  }

  template <typename T> Triple<T>&                      Triple<T>::zero() {

    p[X].zero();
    p[Y].zero();
    p[Z].zero();

    return *this;
  }

  template <typename T> Triple<T>&                      Triple<T>::invalidate() {

    p[X].invalidate();
    p[Y].invalidate();
    p[Z].invalidate();

    return *this;
  }

  template <typename T> bool              Triple<T>::operator== (const Triple<T>& A) const
      { return (p[X] == A[X]) && (p[Y] == A[Y]) && (p[Z] == A[Z]); };

  template <typename T> bool              Triple<T>::operator!= (const Triple<T>& A) const
      { return (p[X] != A[X]) || (p[Y] != A[Y]) || (p[Z] != A[Z]); };


  template<> inline bool                        Triple<double>::operator== (const Triple<double>& A) const
      { return (memcmp (p, A.p, sizeof(p)) == 0); };

  template<> inline bool                        Triple<double>::operator!= (const Triple<double>& A) const
      { return (memcmp (p, A.p, sizeof(p))); };

  template<> inline bool                        Triple<int>::operator== (const Triple<int>& A) const
      { return (memcmp (p, A.p, sizeof(p)) == 0); };

  template<> inline bool                        Triple<int>::operator!= (const Triple<int>& A) const
      { return (memcmp (p, A.p, sizeof(p))); };

  template<> inline bool                        Triple<size_t>::operator== (const Triple<size_t>& A) const
      { return (memcmp (p, A.p, sizeof(p)) == 0); };

  template<> inline bool                        Triple<size_t>::operator!= (const Triple<size_t>& A) const
      { return (memcmp (p, A.p, sizeof(p))); };

  template<> inline bool                        Triple<float>::operator== (const Triple<float>& A) const
      { return (memcmp (p, A.p, sizeof(p)) == 0); };

  template<> inline bool                        Triple<float>::operator!= (const Triple<float>& A) const
      { return (memcmp (p, A.p, sizeof(p))); };



  template <typename T> const std::string              Triple<T>::FILE_EXTENSION = "trp";

  template<typename T> inline Triple<T>  operator* (T M, const Triple<T>& P)
      { return Triple<T> (P[X]*M, P[Y]*(T)M, P[Z]*(T)M); }

  template<typename T> inline Triple< Triple<T> >  operator* (T M, const Triple< Triple<T> >& t)
      { return Triple< Triple<T> > (t[X] * M, t[Y] * M, t[Z] * M); }


  inline Triple<double>   operator+ (const Triple<int>& P, double M)   { return Triple<double> (((double)P[X])+M, ((double)P[Y])+M, ((double)P[Z])+M); }
  inline Triple<double>   operator+ (const Triple<size_t>& P, double M)  { return Triple<double> (((double)P[X])+M, ((double)P[Y])+M, ((double)P[Z])+M); }
  inline Triple<int>      operator+ (const Triple<size_t>& P, int M)    { return Triple<int> (((int)P[X])+M, ((int)P[Y])+M, ((int)P[Z])+M); }
  inline Triple<int>      operator+ (int M, const Triple<size_t>& P)    { return Triple<int> (((int)P[X])+M, ((int)P[Y])+M, ((int)P[Z])+M); }


  inline Triple<double>   operator+ (double M, const Triple<int>& P)   { return Triple<double> (((double)P[X])+M, ((double)P[Y])+M, ((double)P[Z])+M); }
  inline Triple<double>   operator+ (double M, const Triple<size_t>& P)  { return Triple<double> (((double)P[X])+M, ((double)P[Y])+M, ((double)P[Z])+M); }

  inline Triple<double>   operator+ (int M, const Triple<double>& P)   { return Triple<double> (P[X]+(double)M, P[Y]+(double)M, P[Z]+(double)M); }

  inline Triple<double>   operator+ (const Triple<double>& f, const Triple<int>& i)    { return Triple<double>(f[X] + (double)i[X],  f[Y] + (double)i[Y],   f[Z] + (double)i[Z] ); }
  inline Triple<double>   operator+ (const Triple<double>& f, const Triple<size_t>& ui)  { return Triple<double>(f[X] + (double)ui[X], f[Y] + (double)ui[Y],  f[Z] + (double)ui[Z] ); }

  inline Triple<double>   operator+ (const Triple<int>& i, const Triple<double>& f)    { return Triple<double>(f[X] + (double)i[X],  f[Y] + (double)i[Y],   f[Z] + (double)i[Z] ); }
  inline Triple<double>   operator+ (const Triple<size_t>& ui, const Triple<double>& f)  { return Triple<double>(f[X] + (double)ui[X], f[Y] + (double)ui[Y],  f[Z] + (double)ui[Z] ); }

  inline Triple<double>   operator- (const Triple<int>& P, double M)   { return Triple<double> (((double)P[X])-M, ((double)P[Y])-M, ((double)P[Z])-M); }
  inline Triple<double>   operator- (const Triple<size_t>& P, double M)  { return Triple<double> (((double)P[X])-M, ((double)P[Y])-M, ((double)P[Z])-M); }


  inline Triple<double>   operator- (const Triple<double>& f, const Triple<int>& i)    { return Triple<double>(f[X] - (double)i[X],  f[Y] - (double)i[Y],   f[Z] - (double)i[Z] ); }
  inline Triple<double>   operator- (const Triple<double>& f, const Triple<size_t>& ui)  { return Triple<double>(f[X] - (double)ui[X], f[Y] - (double)ui[Y],  f[Z] - (double)ui[Z] ); }

  inline Triple<double>   operator- (const Triple<int>& i, const Triple<double>& f)    { return Triple<double>(f[X] - (double)i[X],  f[Y] - (double)i[Y],   f[Z] - (double)i[Z] ); }
  inline Triple<double>   operator- (const Triple<size_t>& ui, const Triple<double>& f)  { return Triple<double>(f[X] - (double)ui[X], f[Y] - (double)ui[Y],  f[Z] - (double)ui[Z] ); }


  inline Triple<double>   operator* (const Triple<int>& P, double M)   { return Triple<double> (((double)P[X])*M, ((double)P[Y])*M, ((double)P[Z])*M); }
  inline Triple<double>   operator* (const Triple<size_t>& P, double M)  { return Triple<double> (((double)P[X])*M, ((double)P[Y])*M, ((double)P[Z])*M); }

  inline Triple<double>   operator* (double M, const Triple<int>& P)   { return Triple<double> (((double)P[X])*M, ((double)P[Y])*M, ((double)P[Z])*M); }
  inline Triple<double>   operator* (double M, const Triple<size_t>& P)  { return Triple<double> (((double)P[X])*M, ((double)P[Y])*M, ((double)P[Z])*M); }

  inline Triple<double>   operator* (int M, const Triple<double>& P)   { return Triple<double> (P[X]*(double)M, P[Y]*(double)M, P[Z]*(double)M); }

  inline Triple<double>   operator* (const Triple<double>& f, const Triple<int>& i)    { return Triple<double>(f[X] * (double)i[X],  f[Y] * (double)i[Y],   f[Z] * (double)i[Z] ); }
  inline Triple<double>   operator* (const Triple<double>& f, const Triple<size_t>& ui)  { return Triple<double>(f[X] * (double)ui[X], f[Y] * (double)ui[Y],  f[Z] * (double)ui[Z] ); }

  inline Triple<double>   operator* (const Triple<int>& i, const Triple<double>& f)    { return Triple<double>(f[X] * (double)i[X],  f[Y] * (double)i[Y],   f[Z] * (double)i[Z] ); }
  inline Triple<double>   operator* (const Triple<size_t>& ui, const Triple<double>& f)  { return Triple<double>(f[X] * (double)ui[X], f[Y] * (double)ui[Y],  f[Z] * (double)ui[Z] ); }


  inline Triple<double>   operator/ (const Triple<double>& f, const Triple<int>& i)    { return Triple<double>(f[X] / (double)i[X],  f[Y] / (double)i[Y],   f[Z] / (double)i[Z] ); }
  inline Triple<double>   operator/ (const Triple<double>& f, const Triple<size_t>& ui)  { return Triple<double>(f[X] / (double)ui[X], f[Y] / (double)ui[Y],  f[Z] / (double)ui[Z] ); }

  inline Triple<double>   operator/ (const Triple<int>& i, const Triple<double>& f)    { return Triple<double>(f[X] / (double)i[X],  f[Y] / (double)i[Y],   f[Z] / (double)i[Z] ); }
  inline Triple<double>   operator/ (const Triple<size_t>& ui, const Triple<double>& f)  { return Triple<double>(f[X] / (double)ui[X], f[Y] / (double)ui[Y],  f[Z] / (double)ui[Z] ); }

//  template <typename T> Triple<T> operator*(double M, const Triple<T>& t)              { return t * M; }


  //Scalar addition, subtraction, multiplication and division by doubles.

//
//  template <> template <> inline Triple<int>&         Triple<int>::operator+=<double> (double inc)     { throw Exception("Triple<int> is forbidden from scalar addition with  double to avoid accidental use."); }
//  template <> template <> inline Triple<int>&         Triple<int>::operator-=<double> (double inc)     { throw Exception("Triple<int> is forbidden from scalar subtraction with double to avoid accidental use."); }
//  template <> template <> inline Triple<int>&         Triple<int>::operator*=<double> (double M)       { throw Exception("Triple<int> is forbidden from scalar multiplication with double to avoid accidental use."); }
//  template <> template <> inline Triple<int>&         Triple<int>::operator/=<double> (double M)       { throw Exception("Triple<int> is forbidden from scalar division with scalar double to avoid accidental use."); }
//
//  template <> template <> inline Triple<int>          Triple<int>::operator+<double> (double M) const  { throw Exception("Triple<int> is forbidden from scalar addition with double to avoid accidental use."); }
//  template <> template <> inline Triple<int>          Triple<int>::operator-<double> (double M) const  { throw Exception("Triple<int> is forbidden from scalar subtraction with double to avoid accidental use."); }
//  template <> template <> inline Triple<int>          Triple<int>::operator*<double> (double M) const  { throw Exception("Triple<int> is forbidden from scalar multiplication with double to avoid accidental use."); }
//  template <> template <> inline Triple<int>          Triple<int>::operator/<double> (double M) const  { throw Exception("Triple<int> is forbidden from scalar division with scalar double to avoid accidental use."); }
//
//  template <> template <> inline Triple<size_t>&        Triple<size_t>::operator+=<double> (double inc)     { throw Exception("Triple<size_t> is forbidden from scalar addition with  double to avoid accidental use."); }
//  template <> template <> inline Triple<size_t>&        Triple<size_t>::operator-=<double> (double inc)     { throw Exception("Triple<size_t> is forbidden from scalar subtraction with double to avoid accidental use."); }
//  template <> template <> inline Triple<size_t>&        Triple<size_t>::operator*=<double> (double M)       { throw Exception("Triple<size_t> is forbidden from scalar multiplication with double to avoid accidental use."); }
//  template <> template <> inline Triple<size_t>&        Triple<size_t>::operator/=<double> (double M)       { throw Exception("Triple<size_t> is forbidden from scalar division with scalar double to avoid accidental use."); }
//
//  template <> template <> inline Triple<size_t>         Triple<size_t>::operator+<double> (double M) const  { throw Exception("Triple<size_t> is forbidden from scalar addition with double to avoid accidental use."); }
//  template <> template <> inline Triple<size_t>         Triple<size_t>::operator-<double> (double M) const  { throw Exception("Triple<size_t> is forbidden from scalar subtraction with double to avoid accidental use."); }
//  template <> template <> inline Triple<size_t>         Triple<size_t>::operator*<double> (double M) const  { throw Exception("Triple<size_t> is forbidden from scalar multiplication with double to avoid accidental use."); }
//  template <> template <> inline Triple<size_t>         Triple<size_t>::operator/<double> (double M) const  { throw Exception("Triple<size_t> is forbidden from scalar division with scalar double to avoid accidental use."); }
//
//
//  template <typename T> template <typename U> Triple<T>&        Triple<T>::operator+= (U inc)     { p[X] += inc; p[Y] += inc; p[Z] += inc; return *this; }
//  template <typename T> template <typename U> Triple<T>&        Triple<T>::operator-= (U inc)     { p[X] -= inc; p[Y] -= inc; p[Z] -= inc; return *this; }
//  template <typename T> template <typename U> Triple<T>&        Triple<T>::operator*= (U M)       { p[X]*=M; p[Y]*=M; p[Z]*=M; return *this; }
//  template <typename T> template <typename U> Triple<T>&        Triple<T>::operator/= (U M)       { p[X]/=M; p[Y]/=M; p[Z]/=M; return *this; }
//
//  template <typename T> template <typename U> Triple<T>         Triple<T>::operator+ (U M) const  { return (Triple<T> (p[X]+M, p[Y]+M, p[Z]+M)); }
//  template <typename T> template <typename U> Triple<T>         Triple<T>::operator- (U M) const  { return (Triple<T> (p[X]-M, p[Y]-M, p[Z]-M)); }
//  template <typename T> template <typename U> Triple<T>         Triple<T>::operator* (U M) const  { return (Triple<T> (p[X]*M, p[Y]*M, p[Z]*M)); }
//  template <typename T> template <typename U> Triple<T>         Triple<T>::operator/ (U M) const  { return (Triple<T> (p[X]/M, p[Y]/M, p[Z]/M)); }


  template<typename T> inline T              dist2 (const Triple<T>& a, const Triple<T>& b)   { return ((a-b).norm2()); }
  template<typename T> inline T              dist (const Triple<T>& a, const Triple<T>& b)    { return ((a-b).norm()); }




  template <typename T> T                Triple<T>::dot (const Triple<T>& t, Triple<T>& gradient) const {

#ifndef GRADIENT_NOT_REQUIRED
      gradient = *this;
#endif

      return dot(t);

  }


  template <typename T> T                Triple<T>::angle(const Triple<T>& t) const {

    return MR::Math::acos(dot(t)/(norm() * t.norm()));

  }



  inline Triple<double>                  rand_triple(double scale, gsl_rng *rand_gen) {

    Triple<double> triple(0.0,0.0,0.0);

    for (size_t dim_i = 0; dim_i < 3; dim_i++)
      triple[dim_i] = (double)gsl_ran_gaussian(rand_gen, scale);

    return triple;

  }


  template <typename T> Triple<T>                   parse_triple(const std::string& s) {

    Triple<T> triple;

    size_t start_number = 0;
    size_t end_number;

    size_t dim_i = 0;

    start_number = s.find_first_of("0123456789.-");

    while (start_number != std::string::npos && (dim_i < 3) ) {

      end_number = s.find_first_not_of("0123456789.-", start_number);

      std::string num_string = s.substr(start_number, end_number-start_number);


     if (num_string.find_first_of('-', 1) != std::string::npos  || num_string.find('.', num_string.find('.')+1) != std::string::npos)  // If there is a '-' character anywhere other than the first character, or if there is more than one decimal triple then throw an error.
        throw Exception("Coordinate, " + num_string + " at dim: " + str(dim_i+1) + ", is invalid");//

      triple[dim_i] = to<T>(num_string);

      start_number = s.find_first_of("0123456789.-", end_number);
      dim_i++;

    }

    if (dim_i != 3)
      throw Exception("Incorrect number of values in sequence '" + s +"' for conversion to triple.");

    return triple;

  }

}


#endif
