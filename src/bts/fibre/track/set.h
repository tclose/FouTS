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

#ifndef __bts_fibre_track_set_h__
#define __bts_fibre_track_set_h__


#include "bts/common.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track.h"

#include "bts/fibre/base/set.h"
#include "bts/fibre/properties.h"



namespace BTS {

  namespace Fibre {

    class Track::Set : public Base::Set<Track> {

      public:

        typedef Base::SetReader<Track::Set> Reader;
        typedef Base::SetWriter<Track::Set> Writer;

        typedef Base::SetTextReader<Track::Set> TextReader;
        typedef Base::SetTextWriter<Track::Set> TextWriter;

        friend class Object;
        friend class Strand::Set;
        friend class Tractlet::Set;

        friend class Base::SetReader<Track::Set>;
        friend class Base::SetWriter<Track::Set>;

      //Public static constants
      public: 
      
       const static std::string FILE_EXTENSION;
       
       const static std::string ORIGINAL_DEGREE_PROP;       

      protected:

       const static char*          PROPS_LIST[];
    
      //Public methods;
      public:
      

        /*! Initialises a tractlet set with variable tractlet sizes, without base intensity or tractlet intensities
        *
        * @param extended_props The extended properties of the set
        */
        Set (const std::map<std::string,std::string>& extended_props = Properties())
         : Fibre::Base::Set<Track>(std::vector<const char*>(), std::vector<const char*>(), extended_props) {}



        /*! Initialises a track set with variable strand sizes, with properties element propertiess
        *
        * @param props The set properties
        * @param elem_props The element properties
        * @param extended_props The extended properties of the set
        */
        Set (const std::vector<const char*>& props, const std::vector<const char*>& elem_props, const std::map<std::string,std::string>& extended_props = Properties())
         : Fibre::Base::Set<Track>(props, elem_props, extended_props) {}



        /*! Initialises a tractlet set of a preset size with tractlet degree
        *
        * @param size Size of the tractlet set
        * @param num_points The (fixed) number of points per track in the set (can be freed using the 'free_elem_degree()' function)
        * @param extended_props The extended properties of the set
        */
        Set (size_t size, size_t num_points, const std::map<std::string,std::string>& extended_props = Properties())
         : Fibre::Base::Set<Track>(size, num_points, num_points * 3, std::vector<const char*>(), std::vector<const char*>(), extended_props) {}


        /*! Initialises a track set of a preset size with strand degree, with properties and element properties
        *
        * @param size Size of the track set
        * @param props The set properties
        * @param elem_props The element properties
        * @param num_points The (fixed) number of points in the the tracks in the track set (can be freed using the 'free_elem_degree()' function)
        * @param extended_props The extended properties of the set
        */
        Set (size_t size, size_t num_points, const std::vector<const char*>& props, const std::vector<const char*>& elem_props, const std::map<std::string,std::string>& extended_props = Properties())
         : Fibre::Base::Set<Track>(size, num_points, num_points * 3 + elem_props.size(), props, elem_props, extended_props) {}


//
//        Set (const std::map<std::string,std::string>& properties = Properties())
//          : Base::Set<Track>(properties) {}
//
//        Set (size_t num_tcks, const std::map<std::string,std::string>& properties = Properties())
//          : Base::Set<Track>(num_tcks, properties) {}
//
//        Set (size_t num_tcks, size_t degree, const std::map<std::string,std::string>& properties = Properties())
//          : Base::Set<Track>(num_tcks, degree, properties) {}
//
//        Set (size_t num_tcks, const Track& default_value, const std::map<std::string,std::string>& properties = Properties())
//          : Base::Set<Track>(num_tcks, default_value, properties) {}



        Set (const Set& s_set)
          : Base::Set<Track>(s_set) {}

        Set (std::string file_location, size_t num_points = 0, double strands_per_width = -1.0)
          : Fibre::Base::Set<Track>(std::vector<const char*>(), std::vector<const char*>(), std::map<std::string,std::string>())
          { load(file_location, num_points, strands_per_width); }

        Set (std::string file_location, double strands_per_width)
          : Fibre::Base::Set<Track>(std::vector<const char*>(), std::vector<const char*>(), std::map<std::string,std::string>())
          { load(file_location, strands_per_width); }

        Set (const Strand::Set& strands, size_t num_length_sections = 0);


        Set (const Tractlet::Set& t, size_t num_length_sections, size_t num_width_sections);


        Set (const Tractlet::Set& t, size_t num_length_sections, double strands_per_acs);


        BASE_GENERAL_FUNCTIONS(Set);

        BASE_SET_FUNCTIONS(Set, Track);


        Strand::Set                       to_strands(size_t degree) const;

        Tractlet::Set                     to_tractlets(size_t degree) const;

        Tractlet::Set                     to_tractlets(size_t degree, size_t num_tractlets) const;

        Tractlet::Set                     to_tractlets(size_t degree, double width) const;

        Set                               peel(double distance, bool to_cube = false);

        Set                               trim(double sphere_radius);

        void                              remove_short_tracks(size_t min_num_points);

        void                              load(const std::string& path, size_t num_points = 0, double strands_per_width = -1.0);

        void                              load(const std::string& path, double strands_per_width)
          { load(path, 0, strands_per_width); }

        void                              save(const std::string& path, size_t degree = 0, size_t num_tractlets = 0) const;

        
        Set&                              operator+=  (const Set& set)
          { Base::Set<Track>::operator+=(set); return *this; }

        Set&                              operator-=  (const Set& set)
          { Base::Set<Track>::operator-=(set); return *this; }

        Set&                              operator*= (const Set& set)
          { Base::Set<Track>::operator*=(set); return *this; }

        Set&                              operator*= (double c)
          { Base::Set<Track>::operator*=(c); return *this; }


        Set&                              operator/= (const Set& set)
          { Base::Set<Track>::operator/=(set); return *this; }

        Set&                              operator/= (double c)
          { Base::Set<Track>::operator/=(c); return *this; }


        Set                               operator+  (const Set& set) const
          { Set answer (*this); answer += set; return answer; }

        Set                               operator-  (const Set& set) const
          { Set answer (*this); answer -= set; return answer; }

        Set                               operator* (const Set& set) const
          { Set answer (*this); answer *= set; return answer; }
        
        Set                               operator* (double c) const
          { Set answer (*this); answer *= c; return answer; }


        Set                               operator/ (const Set& set) const
          { Set answer (*this); answer /= set; return answer; }

        Set                               operator/ (double c) const
          { Set answer (*this); answer /= c; return answer; }

      private:

        void                              from_strands(const Strand::Set& strands, size_t num_length_sections);

    };



    inline Track::Set                           operator* (double c, Track::Set set)
      { Track::Set answer (set); answer *= c;    return answer; }

  }

}


#endif
