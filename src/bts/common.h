/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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


#ifndef __bts_common_h__
#define __bts_common_h__


//Defines the parameters
#define COMMON_PARAMETERS \
  Option ("comments", "Description of the fibres/image, and/or the purpose of the performed experiment.")  \
     + Argument ("comments", "").type_text() \

  
//Loads the parameters into variables
#define SET_COMMON_PARAMETERS \
  std::string cmn_comments                = ""; \
\
   Options cmn_opt = get_options("comments"); \
  if (cmn_opt.size()) \
    cmn_comments = cmn_opt[0][0].c_str(); \



//Adds the parameters to the properties to be saved with the data.
#define ADD_COMMON_PROPERTIES(properties) \
  properties["comments"]        = cmn_comments; \


#include <string>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <linux/limits.h>

#include "types.h"
#include "mrtrix.h"


#undef X
#undef Y
#undef Z
#undef DW
#undef AZ
#undef EL
#undef RAD
#undef INV_M_PI
#undef SQRT_3
#undef SQRT_2
#undef M_PI_ON_2


namespace BTS {

  //Standard dimension order
  const size_t X = 0;
  const size_t Y = 1;
  const size_t Z = 2;
  const size_t DW = 3;

  //Standard spherical dimension order
  const size_t EL = 0;
  const size_t AZ = 1;
  const size_t RAD = 2;

  //Common mathematical constants
  const double INV_M_PI            = 0.318309886;
  const double GAUSS_NORM          = 0.398942280;
  const double JUST_UNDER_PI_ON_2  = 1.5706;
  const double SQRT_3              = 1.732050807;
  const double SQRT_2              = 1.414213562;
  const double SQRT_LOG_2          = 0.832554611;
  const double LOG_2               = 0.693147181;
  const double LOG_PI              = 1.144729886;

  const std::string PARAM_DIRECTORY = "/home/tclose/data/";

  const double LARGE_FLOAT  = 1e10;  //A large double still with a bit of room to move.
  const double SMALL_FLOAT  = 1e-10; //A small double still with a bit of room to move.
  const size_t  LARGE_INT     = 1e8;   //A large number still with a bit of room to move.

  typedef MR::Exception Exception;

//  using   std::cout;
//  using   std::endl;

}

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}


#include <iomanip>
#include <vector>
#include "file/path.h"
#include "exception.h"
#include "math/math.h"
#include "math/vector.h"


namespace BTS {

  template<typename T, typename U> class CompareFirst { public: bool operator()(const std::pair<T,U>& p1, const std::pair<T,U>& p2) { return p1.first < p2.first; }};
  template<typename T, typename U> class CompareSecond { public: bool operator()(const std::pair<T,U>& p1, const std::pair<T,U>& p2) { return p1.second < p2.second; }};

  template <typename T> T               min2 (T const x, T const y) throw ()        { return (x < y)? x: y;}


  template <typename T> T               max2 (T const x, T const y) throw ()        { return (x > y)? x: y;}


  template <typename T> void            zero(T& element) {
    element.zero();
  }


  inline void                           zero(float& element) {
    element = 0.0;
  }


  inline void                           zero(double& element) {
    element = 0.0;
  }


  inline void                           zero(int& element) {
    element = 0;
  }


  inline void                           zero(size_t& element) {
    element = 0;
  }


  template <typename T> void            negate(T& value) {
    value.negate();
  }


  inline void                           negate(float& value) {
    value = -value;
  }


  inline void                           negate(double& value) {
    value = -value;
  }


  inline void                           negate(int& value) {
    value = -value;;
  }


  inline void                           negate(size_t& value) {
    value = -value;;
  }


  template <typename T> void            invalidate(T& element) {
    element.invalidate();
  }


  inline void                           invalidate(float& element) {
    element = NAN;
  }


  inline void                           invalidate(double& element) {
    element = NAN;
  }


  inline void                           invalidate(int& element) {
    element = std::numeric_limits<int>::max();
  }


  template <typename T> bool            is_empty(const T& element) {
    return !element.size();
  }


  inline bool                           is_empty(double dble) {
    return isnan(dble);
  }



  void                                  breakpoint();


  double                                rand_double();


  std::vector<double>                   parse_double_vector(const std::string& str);


  std::vector<int>                      parse_int_vector(const std::string& str);


  std::vector<size_t>                   parse_size_t_vector(const std::string& str);


  int                                   num_dec_places(int number);


  int                                   num_precision_places(double number);

  
  template <typename T> std::string str(const T& value)  {
    std::stringstream ss;
    ss << value;
    return ss.str();
  }


  inline std::string                    str (int number, size_t dec_places) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(dec_places + (number < 0)) << number;
    return ss.str();
  }


  inline std::string                    str (double number, size_t dec_places, size_t precision) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << std::setfill('0') << std::setw(dec_places + precision + (precision > 0) + (number < 0)) << number;
    return ss.str();
  }
  

  template <typename T> T               to(const std::string &s) {
    
    //Copied from http://snippets.dzone.com/posts/show/347
  
    T t;
    std::istringstream iss(s);
    iss >> t;
    return t;
  }
  

  inline std::string                    current_datetime(std::string format_string = "%c") {
  
    //Adapted from example at http://www.cplusplus.com/reference/clibrary/ctime/strftime/ .
  
    time_t rawtime;
    struct tm *timeinfo;
    char current_time[200];

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    
    
    strftime (current_time, 200, format_string.c_str(), timeinfo);
    
    return str(current_time);
  
  }

  // Used to sort vectors of const char*
  inline bool                           cstr_sort(const char* a, const char* b) {

    return strcmp(a, b) < 0;

  }


  template <class T>
  struct PrintMap : public std::unary_function<T, void>
  {
      std::ostream& os;
      PrintMap(std::ostream& strm) : os(strm) {}

      void operator()(const T& elem) const
      {
          os << elem.first << ", " << elem.second << std::endl;
      }
  };

  template <typename T> std::vector<T> parse_sequence (const std::string& spec, T last = std::numeric_limits<T>::max()) {
    std::vector<T> V;
    if (!spec.size()) throw Exception ("unsigned integer sequence specifier is empty");
    std::string::size_type start = 0, end;
    T num[3];
    T i = 0;
    try {
     do {
       end = spec.find_first_of (",:", start);
       std::string str (MR::strip (spec.substr (start, end-start)));
       MR::lowercase (str);
       if (str == "end") {
         if (last == std::numeric_limits<T>::max())
           throw Exception ("value of \"end\" is not known in number sequence \"" + spec + "\"");
         num[i] = last;
       }
       else num[i] = to<T> (spec.substr (start, end-start));

       char last_char = end < spec.size() ? spec[end] : '\0';
       if (last_char == ':') { i++; if (i > 2) throw Exception ("invalid number range in number sequence \"" + spec + "\""); }
       else {
         if (i) {
           T inc, last;
           if (i == 2) { inc = num[1]; last = num[2]; }
           else { inc = 1; last = num[1]; }
           if (inc * (last - num[0]) < 0) inc = -inc;
           for (; ( inc > 0 ? num[0] <= last : num[0] >= last ) ; num[0] += inc) V.push_back (num[0]);
         }
         else V.push_back (num[0]);
         i = 0;
       }

       start = end+1;
     } while (end < spec.size());
    }
    catch (Exception& E) {
     throw Exception (E, "can't parse unsigned integer sequence specifier \"" + spec + "\"");
    }

    return (V);

  }


  inline std::string                          cwd() {
     char buff[PATH_MAX];
     if (!getcwd(buff, PATH_MAX))
       throw Exception ("Could not get the cwd directory");
     std::string work_dir(buff);
     return work_dir;
  }

  // A gdb pretty printing hack
  class CoordView {};

}

#endif

