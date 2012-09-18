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


#include "bts/fibre/track.h"

#ifndef __bts_fibre_strand_set_h__
#define __bts_fibre_strand_set_h__


#include "bts/common.h"
#include "bts/file.h"

#include "bts/fibre/strand.h"
#include "bts/fibre/strand/section.h"
#include "bts/fibre/base/set.h"
#include "bts/fibre/tractlet.h"
#include "bts/math/munkres.h"

#include "bts/fibre/base/set_reader.h"
#include "bts/fibre/base/set_writer.h"
#include "bts/fibre/properties.h"

#include "k_means/KMlocal.h"	// k-means algorithms


namespace BTS {

  namespace Fibre {

    class Strand::Set : public Fibre::Base::Set<Strand> {

      //Public Nested classes and typedefs
      public:

        class Walker;
        class Momentum;
        typedef Base::Tensor<Strand::Set> Tensor;
        typedef Base::SetReader<Strand::Set> Reader;
        typedef Base::SetWriter<Strand::Set> Writer;

        typedef Base::SetTextReader<Strand::Set> TextReader;
        typedef Base::SetTextWriter<Strand::Set> TextWriter;

        friend class Object;
        friend class Track::Set;
        friend class Tractlet::Set;
        friend class Base::Tensor<Strand::Set>;
        friend class Base::SetReader<Strand::Set>;
        friend class Base::SetWriter<Strand::Set>;

        typedef Strand Element;


      //Public static members
      public:

        const static std::string FILE_EXTENSION;
        const static std::string ORIGINAL_NUM_POINTS_PROP;
        const static std::string ORIGINAL_NUM_TRACTLETS_PROP;


        // Clustering parameters used in reduce_number.

        const static size_t   NUM_STAGES = 100;          // run for 100 stages
        const static double   MAX_ACCUM_RDL;             // min consec RDL
        const static double   MIN_ACCUM_RDL;             // min accum RDL
        const static int      MAX_RUN_STAGES = 3;        // max run stages
        const static double   INIT_PROB_ACCEPTANCE;      // init. prob. OF ACCEPTance
        const static int      TEMP_RUN_LENGTH = 10;      // temp. RUN LENGTH
        const static double   TEMP_REDUC_FACTOR;         // temp. reduction factor

      protected:

        const static char*          PROPS_LIST[];

  
      //Public static functions
      public:

        //Do nothing, no characteristic properties in Strand::Set class.
        static std::vector<std::string>&  append_characteristic_property_keys(std::vector<std::string>& header) { return header; }

      //Public member functions
      public:

        //TODO: Make the properties of these constructors protected, add all Fibre objects to the 'friends' list of each other, and add constructors with flags to determine whether certain properties are included.

        /*! Initialises a strand set with variable strand sizes, without properties or element properties
         *
         * @param extended_props The extended properties of the set
         */
        Set (const std::map<std::string,std::string>& extended_props = Properties())
          : Fibre::Base::Set<Strand>(std::vector<const char*>(), std::vector<const char*>(), extended_props) {}


        /*! Initialises a strand set with variable strand sizes, with properties element propertiess
         *
         * @param props The set properties
         * @param elem_props The element properties
         * @param extended_props The extended properties of the set
         */
        Set (const std::vector<const char*>& props, const std::vector<const char*>& elem_props, const std::map<std::string,std::string>& extended_props = Properties())
          : Fibre::Base::Set<Strand>(props, elem_props, extended_props) {}


        /*! Initialises a strand set of a preset size with strand degree, without properties or element properties
         *
         * @param size Size of the strand set
         * @param degree The (fixed) degree of the the strands in the strand set (can be freed using the 'free_elem_degree()' function)
         * @param extended_props The extended properties of the set
         */
        Set (size_t size, size_t degree, const std::map<std::string,std::string>& extended_props = Properties())
          : Fibre::Base::Set<Strand>(size, degree, degree * 3, std::vector<const char*>(), std::vector<const char*>(), extended_props) {}


        /*! Initialises a strand set of a preset size with strand degree, with properties and element properties
         *
         * @param size Size of the strand set
         * @param props The set properties
         * @param elem_props The element properties
         * @param degree The (fixed) degree of the the strands in the strand set (can be freed using the 'free_elem_degree()' function)
         * @param extended_props The extended properties of the set
         */
        Set (size_t size, size_t degree, const std::vector<const char*>& props, const std::vector<const char*>& elem_props, const std::map<std::string,std::string>& extended_props = Properties())
          : Fibre::Base::Set<Strand>(size, degree, degree * 3 + elem_props.size(), props, elem_props, extended_props) {}


        Set (std::string file_location, size_t degree = 0, double strands_per_width = -1.0)
          : Fibre::Base::Set<Strand>(std::vector<const char*>(), std::vector<const char*>(), std::map<std::string,std::string>())
          { load(file_location, degree, strands_per_width);}


        Set (std::string file_location, double strands_per_width)
          : Fibre::Base::Set<Strand>(std::vector<const char*>(), std::vector<const char*>(), std::map<std::string,std::string>())
          { load(file_location, strands_per_width);}


        //TODO: Switch this around the other way, so that the to_strands function calls this constructor.
        Set (const Track::Set& tcks, size_t degree = 0);


        Set (const Tractlet::Set& set, size_t num_samples);


        Set (const Tractlet::Set& set, double strands_per_vol);


        Set (const Tractlet& set, size_t num_samples);


        Set (const Tractlet& set, double strands_per_vol);


        //        Set (size_t num_strands, const std::map<std::string,std::string>& properties = Properties())
        //          : Fibre::Base::Set<Strand>(num_strands, properties, 1.0) {}
        //
        //        Set (size_t num_strands, double base_intensity, const std::map<std::string,std::string>& properties = Properties())
        //          : Fibre::Base::Set<Strand>(num_strands, properties, base_intensity) {}
        //
        //        Set (size_t num_strands, size_t degree, double base_intensity = 1.0, const std::map<std::string,std::string>& properties = Properties())
        //          : Fibre::Base::Set<Strand>(num_strands, degree, properties, base_intensity) {}
        //
        //        Set (size_t num_strands, const Strand& default_value, double base_intensity = 1.0, const std::map<std::string,std::string>& properties = Properties())
        //          : Fibre::Base::Set<Strand>(num_strands, default_value, properties, base_intensity) {}


        Set (const Set& s_set)
          : Fibre::Base::Set<Strand>(s_set) {}


        BASE_GENERAL_FUNCTIONS(Set);

        BASE_SET_FUNCTIONS(Set);

        Track::Set                                  to_tracks (size_t num_points) const;

        //!Uses Track::BUNDLE_INDEX property to group strands into bundles and convert the bundles into tractlets.
        Tractlet::Set                               to_tractlets() const;

        //!First clusters strands into bundles then converts the bundles into tractlets.
        //!NB:If num_tractlets argument is less than zero reverts back to to_tractlets() function. \arg num_tractlets number of tractlets to cluster
        Tractlet::Set                               to_tractlets (size_t num_tractlets, const std::string& cluster_save_location = "") const;
        
        //!Creates a new tractlet for each strand \arg width width of the created tractlets.
        Tractlet::Set                               to_tractlets (double width) const;


        Strand::Set&                                cluster(size_t num_clusters);

//        Permutation                                 default_permutation();
        
        Set                                         permute(const std::vector<size_t>& indices) const;

        Set                                         smallest_distance_set(const Set& reference) const
          { Set smallest; this->Fibre::Base::Set<Strand>::smallest_distance_set(reference, smallest); return smallest; }

        Set                                         shortest_span() const;

        Set                                         flip(const std::vector<size_t>& to_flip) const;

        Strand                                      centroid() const;

        Set                                         peel(double distance, size_t num_length_sections, bool to_cube = false) const;

        Set                                         reduce_number(size_t new_num);

        std::vector<size_t>                         bundle(size_t index) const;
        
        Set                                         bundle_backbones() const;
        
        std::vector<size_t>                         bundle_indices() const;
            
        void                                        project_to_sphere(double sphere_radius);

        size_t                                      max_degree() const;

        size_t                                      degree() const
        { if (var_elem_degrees()) throw Exception ("Strand set has variable rows, so its 'degree' is undefined"); return elem_dgree; }


        void                                        load(const std::string& path, size_t degree = 0, double strands_per_width = -1.0);

        void                                        load(const std::string& path, double strands_per_width)
          { load(path, 0, strands_per_width); }

        void                                        save(const std::string& path, size_t num_sections = 0, size_t num_tractlets = 0) const;

        
        Set&                                        operator+=  (const Set& set)                { Fibre::Base::Set<Strand>::operator+=(set); return *this; }
        Set&                                        operator-=  (const Set& set)                { Fibre::Base::Set<Strand>::operator-=(set); return *this; }

        Set&                                        operator*=  (const Set& set)                { Fibre::Base::Set<Strand>::operator*=(set); return *this; }
        Set&                                        operator*= (double c)                       { Fibre::Base::Set<Strand>::operator*=(c); return *this; }

        Set&                                        operator/=  (const Set& set)                { Fibre::Base::Set<Strand>::operator/=(set); return *this; }
        Set&                                        operator/= (double c)                       { Fibre::Base::Set<Strand>::operator/=(c); return *this; }

        Set                                         operator+  (const Set& set) const           { Set answer (*this); answer += set; return answer; }
        Set                                         operator-  (const Set& set) const           { Set answer (*this); answer -= set; return answer; }

        Set                                         operator* (const Set& set) const            { Set answer (*this); answer *= set;  return *this; }
        Set                                         operator* (double c) const                  { Set answer (*this); answer *= c; return answer; }
        
        Set                                         operator/ (const Set& set) const            { Set answer (*this); answer /= set; return *this; }
        Set                                         operator/ (double c) const                  { Set answer (*this); answer /= c; return answer; }
        
        
        inline Set&                                 zero()                                      {  Fibre::Base::Set<Strand>::zero();return *this; }
        
        std::string                                 matlab_str() const;

        void                                        calc_characteristic_properties() {}

      private:

        void from_tractlet(const Tractlet& tractlet, size_t num_width_sections);

    };

    double                                         distance(Strand::Set& strands, Strand::Set& reference, std::vector<size_t>& matched_indices);

    double                                         distance(Strand::Set& strands, Strand::Set& reference, std::vector<size_t>& matched_indices, BTS::Math::Munkres<double>&munkres_algorithm, MR::Math::Matrix<double>& similarity);

    inline Strand::Set                             operator* (double c, Strand::Set set)
      { Strand::Set answer (set); answer *= c;    return answer; }


  }

}


#endif

