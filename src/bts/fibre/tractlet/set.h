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

#include "bts/fibre/base/set.h"


#ifndef __bts_fibre_tractlet_set_h__
#define __bts_fibre_tractlet_set_h__

#include <vector>

#include "bts/fibre/track.h"
#include "bts/fibre/tractlet.h"
#include "bts/fibre/base/set_reader.h"
#include "bts/fibre/base/set_writer.h"
#include "bts/fibre/properties.h"

#include "bts/common.h"
#include "bts/file.h"


#include "bts/math/munkres.h"

#include "k_means/KMlocal.h"			// k-means algorithms

#define LOOP(op) \
  for (size_t tractlet_i = 0; tractlet_i < size(); ++tractlet_i) { op }

namespace BTS {

  namespace Fibre {

    class Tractlet::Set : public Base::Set<Tractlet> {

      //Public nested classes and typedefs
      public:

        class Walker;
        class Momentum;
        typedef Base::Tensor<Tractlet::Set> Tensor;
        typedef Base::SetReader<Tractlet::Set> Reader;
        typedef Base::SetWriter<Tractlet::Set> Writer;

        typedef Base::SetTextReader<Tractlet::Set> TextReader;
        typedef Base::SetTextWriter<Tractlet::Set> TextWriter;

        friend class Object;
        friend class Base::Tensor<Tractlet::Set>;
        friend class Base::SetReader<Tractlet::Set>;
        friend class Base::SetWriter<Tractlet::Set>;
        friend class Track::Set;
        friend class Strand::Set;

      public:

        std::vector<std::string>&   append_characteristic_keys(std::vector<std::string>& header) const;


      //Public static constants
      public:

        const static std::string FILE_EXTENSION;

        const static char*          LENGTH_EPSILON_PROP;
        const static char*          WIDTH_EPSILON_PROP;
        const static double         LENGTH_EPSILON_DEFAULT;
        const static double         WIDTH_EPSILON_DEFAULT;

      protected:

        const static char*          PROPS_LIST[];
        const static char*          CHARACTERISTIC_PROPS_LIST[];

      //Public methods;
      public:

        /*! Initialises a tractlet set with variable tractlet sizes, without base intensity or tractlet intensities
         *
         * @param extended_props The extended properties of the set
         */
        Set (const std::map<std::string,std::string>& extended_props = Properties())
          : Base::Set<Tractlet>(std::vector<const char*>(), std::vector<const char*>(), extended_props) {}


        /*! Initialises a tractlet set with variable strand sizes, with properties element propertiess
         *
         * @param props The set properties
         * @param elem_props The element properties
         * @param extended_props The extended properties of the set
         */
        Set (const std::vector<const char*>& props, const std::vector<const char*>& elem_props, const std::map<std::string,std::string>& extended_props = Properties())
          : Base::Set<Tractlet>(props, elem_props, extended_props) {}


        /*! Initialises a tractlet set of a preset size with tractlet degree
         *
         * @param size Size of the tractlet set
         * @param degree The (fixed) degree of the the tractlets in the tractlet set (can be freed using the 'free_elem_degree()' function)
         * @param extended_props The extended properties of the set
         */
        Set (size_t size, size_t degree, const std::map<std::string,std::string>& extended_props = Properties())
          : Base::Set<Tractlet>(size, degree, 9 * degree, std::vector<const char*>(), std::vector<const char*>(), extended_props) {}


        /*! Initialises a tractlet set of a preset size with strand degree, with properties and element properties
         *
         * @param size Size of the tractlet set
         * @param props The set properties
         * @param elem_props The element properties
         * @param degree The (fixed) degree of the the tractlets in the tractlet set (can be freed using the 'free_elem_degree()' function)
         * @param extended_props The extended properties of the set
         */
        Set (size_t size, size_t degree, const std::vector<const char*>& props, const std::vector<const char*>& elem_props, const std::map<std::string,std::string>& extended_props = Properties())
          : Base::Set<Tractlet>(size, degree, 9 * degree + elem_props.size(), props, elem_props, extended_props) {}


        Set (std::string file_location, size_t degree = 0, size_t num_tractlets = 0)
          : Base::Set<Tractlet>(std::vector<const char*>(), std::vector<const char*>(), std::map<std::string,std::string>())
          { load(file_location, degree, num_tractlets); }

        Set (const Set& s_set)
          : Base::Set<Tractlet>(s_set) {}

        Set(const Strand::Set& set);

        Set(const Strand::Set& strands, double width);
        
        Set(const Strand::Set& set, size_t num_tractlets, const std::string& cluster_save_location = "");

        Set (const Track::Set&, size_t degree);

        Set (const Track::Set& set, size_t degree, double width);

        Set (const Track::Set& set, size_t degree, size_t num_tractlets, const std::string& cluster_save_location = "");

        BASE_GENERAL_FUNCTIONS(Set);

        BASE_SET_FUNCTIONS(Set);

        size_t                                      max_degree() const;

        size_t                                      degree() const
          { if (var_elem_degrees()) throw Exception ("Tractlet set has variable rows, so its 'degree' is undefined"); return elem_dgree; }

        void                                        set_width_epsilon(double width_epsilon);

        void                                        set_length_epsilon(double length_epsilon);

        void                                        project_to_sphere(double sphere_radius);        

        void                                        load(const std::string& location, size_t degree = 0, size_t num_tractlets = 0);

        void                                        save(const std::string& location, size_t num_length_sections = 0, double strands_per_acs = -1.0) const;

        void                                        save(const std::string& location, double strands_per_acs) const
          { save(location, 0, strands_per_acs); }
        
        //!Converts to a constant number of strands per tractlet
        Strand::Set                                 to_strands(size_t num_length_sections) const;

        //!Converts to a constant number of strands per tractlet cross-sectional area
        Strand::Set                                 to_strands(double strands_per_acs) const;

        Track::Set                                  to_tracks(size_t num_length_sections, size_t num_width_sections) const;
        
        Track::Set                                  to_tracks(size_t num_length_sections, double strands_per_acs) const;

        Tractlet::Set                               peel(double inclusion_radius, size_t num_length_sections, size_t num_width_sections, bool to_cube = false) const;

        Strand::Set                                 backbones() const;

        Tractlet::Set&                              sanitize();

        Set                                         permute(const std::vector<size_t>& indices) const;

        Set                                         smallest_distance_set(const Set& reference)
          { Set smallest; Base::Set<Tractlet>::smallest_distance_set(reference, smallest); return smallest; }

        Set&                                        normalise_densities(double num_points = 100);

        Set&                                        operator+=  (const Set& set)
          { Base::Set<Tractlet>::operator+=(set); return *this; }

        Set&                                        operator-=  (const Set& set)
          { Base::Set<Tractlet>::operator-=(set); return *this; }

        Set&                                        operator*= (const Set& set)
          { Base::Set<Tractlet>::operator*=(set); return *this; }

        Set&                                        operator*= (double c)
          { Base::Set<Tractlet>::operator*=(c); return *this; }

        Set&                                        operator/= (const Set& set)
          { Base::Set<Tractlet>::operator/=(set); return *this; }

        Set&                                        operator/= (double c)
          { Base::Set<Tractlet>::operator/=(c); return *this; }

        Set                                         operator+  (const Set& set) const
          { Set answer (*this); answer += set; return answer; }

        Set                                         operator-  (const Set& set) const
          { Set answer (*this); answer -= set; return answer; }

        Set                                         operator* (const Set& set) const
          { Set answer (*this); answer *= set; return answer; }
        
        Set                                         operator* (double c) const
          { Set answer (*this); answer *= c; return answer; }

        Set                                         operator/ (double c) const
          { Set answer (*this); answer /= c; return answer; }

        Set                                         operator/ (const Set& set) const
          { Set answer (*this); answer /= set; return answer; }

        std::string                                 matlab_str() const;

        void                                        set_characteristics();



      private:

        /*! Initiates the tractlet set from a set of strands by first clustering the strands into clusters and then converting
         * the clusters into tractlets.
         *
         * @param strands The strands to be converted
         */
        void                                        from_strands(const Strand::Set& strands);

        /*! Initiates the tractlet set from a set of strands by first clustering the strands into a number of clusters
         *  specified by the @arg num_clusters parameter, and then converting the clusters into tractlets.
         *
         * @param strands The strands to be converted
         * @param num_clusters The number of new clusters to generate
         */
        void                                        from_strands(const Strand::Set& strands, size_t num_clusters, const std::string& cluster_save_location = "");


        /*! Initiates the tractlet set from a set of strands by converting each strand into a new tractlet with given width
         *
         * @param strands The strands to be converted
         * @param width The width of the new tractlets about the provided strand 'backbones'
         */
        void                                        from_strands(const Strand::Set& strands, double width);

    };


    inline Tractlet::Set                            operator* (double c, Tractlet::Set set)
      { Tractlet::Set answer (set); answer *= c;    return answer; }


    inline Tractlet::Set&                          Tractlet::Set::normalise_densities(double num_points) {
      if (!has_elem_prop(Tractlet::ALPHA_PROP))
        add_elem_prop(Tractlet::ALPHA_PROP);
      LOOP(operator[](tractlet_i).normalise_density(num_points);)
      return *this;
    }


  }

}

#undef LOOP

#endif


