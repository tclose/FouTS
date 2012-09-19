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


#include <set>

#include "bts/fibre/strand/set.h"
//#include "bts/fibre/strand/set/permutation.h"

#include "bts/fibre/strand.h"


#include "bts/fibre/track.h"
#include "bts/fibre/track/set.h"

#include "bts/fibre/tractlet.h"
#include "bts/fibre/tractlet/set.h"


#include "bts/math/munkres.h"
#include "bts/math/svd.h"
#include "bts/math/common.h"





#include "bts/fibre/base/reader.cpp.h"
#include "bts/fibre/base/writer.cpp.h"

#include "phantom/shared/strand_collection.h"


namespace BTS {

  namespace Fibre {

    const std::string                           Strand::Set::FILE_EXTENSION = "sst";
    const char*                                 Strand::Set::PROPS_LIST[] = { Base::Set<Strand>::BASE_INTENSITY_PROP, PROPS_LIST_END };


    const double                                Strand::Set::MAX_ACCUM_RDL = 0.1;
    const double                                Strand::Set::MIN_ACCUM_RDL = 0.1;
    const double                                Strand::Set::TEMP_REDUC_FACTOR = 0.95;
    const double                                Strand::Set::INIT_PROB_ACCEPTANCE = 0.5;
    
    const std::string                           Strand::Set::ORIGINAL_NUM_POINTS_PROP = "loaded_num_points";
    const std::string                           Strand::Set::ORIGINAL_NUM_TRACTLETS_PROP = "loaded_num_tractlets";


    size_t                                        MAX_ITERATIONS = 1000;



    Strand::Set::Set (const Track::Set& tcks, size_t degree)
      : Base::Set<Strand>(tcks.size(),
                          degree,
                          3 * degree + select_props<Strand>(*tcks.elem_props).size(),
                          select_props<Set>(*tcks.props),
                          select_props<Strand>(*tcks.elem_props),
                          (tcks.ext_props ? *tcks.ext_props : std::map<std::string,std::string>())) {

      for (size_t tck_i = 0; tck_i < size(); tck_i++) {

        size_t strand_degree = degree;

        if (!degree) {
          if (has_extend_elem_prop(Track::Set::ORIGINAL_DEGREE_PROP))
            strand_degree = get_extend_elem_prop<size_t>(Track::Set::ORIGINAL_DEGREE_PROP, tck_i);
          else
            throw Exception ("Degree must be supplied (and cannot be 0) for conversion from tracks to strands.");
        }

        operator[](tck_i) = Strand(tcks[tck_i], strand_degree); //tcks[tck_i].to_strand(strand_degree);

      }

      for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
        prop(prop_i) = tcks.prop(prop_key(prop_i));

      copy_extend_elem_props(tcks);

    }



    Strand::Set::Set(const Tractlet::Set& set, size_t num_samples)
    : Fibre::Base::Set<Strand>(select_props<Set>(*set.props), select_props<Strand>(*set.elem_props), set.get_extend_props()) {

      add_extend_elem_prop(BUNDLE_INDEX_EPROP, "-1");

      for (size_t tractlet_i = 0; tractlet_i < set.size(); tractlet_i++) {

        Set new_strands (set[tractlet_i], num_samples);

        if (new_strands.size()) {

          size_t old_size = size();

          resize(size() + new_strands.size(), new_strands[0]);

          for (size_t strand_i = 0; strand_i < new_strands.size(); strand_i++) {
            operator[](strand_i + old_size) = new_strands[strand_i];
            set_extend_elem_prop(BUNDLE_INDEX_EPROP, str(tractlet_i), strand_i + old_size);
          }
        }
      }

      for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
        prop(prop_i) = set.prop(prop_key(prop_i));


      if (!set.var_elem_degrees())
        freeze_elem_degree();

    }


    Strand::Set::Set(const Tractlet::Set& set, double strands_per_vol)
      : Fibre::Base::Set<Strand>(select_props<Set>(*set.props), select_props<Strand>(*set.elem_props), set.get_extend_props()) {

      add_extend_elem_prop(BUNDLE_INDEX_EPROP, "-1");

      for (size_t tract_i = 0; tract_i < set.size(); tract_i++) {

        std::string bundle_index = set.has_extend_elem_prop(BUNDLE_INDEX_EPROP) ? set.get_extend_elem_prop(BUNDLE_INDEX_EPROP, tract_i) : str(tract_i);

        std::map<std::string,std::string> prop_row;
        prop_row[BUNDLE_INDEX_EPROP] = bundle_index;

        Set new_strands (set[tract_i], strands_per_vol);

        for (size_t new_strand_i = 0; new_strand_i < new_strands.size(); ++new_strand_i)
          push_back(new_strands[new_strand_i], prop_row);

      }

      for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
        prop(prop_i) = set.prop(prop_key(prop_i));

      if (!set.var_elem_degrees())
        freeze_elem_degree(); //Attempt to freeze element size but ignore if fails.

    }


    Strand::Set::Set(const Tractlet& tractlet, size_t num_width_sections)
      : Fibre::Base::Set<Strand>(std::vector<const char*>(), select_props<Strand>(*tractlet.props), std::map<std::string,std::string>()) {

      from_tractlet(tractlet, num_width_sections);

    }


    Strand::Set::Set(const Tractlet& tractlet, double strands_per_acs)
      : Fibre::Base::Set<Strand>(std::vector<const char*>(), select_props<Strand>(*tractlet.props), std::map<std::string,std::string>()) {

      double num_strands = strands_per_acs * tractlet.acs();

      size_t num_width_sections = (size_t)MR::Math::round(MR::Math::sqrt(num_strands * 4.0 / M_PI));

      from_tractlet(tractlet, num_width_sections);

    }


    void   Strand::Set::from_tractlet (const Tractlet& tractlet, size_t num_width_sections) {

      double width_fraction = 1.0 / (double)num_width_sections;

      size_t count = 0;

      // Loop across cross-section of tractlet. Loop each perturbation axis from a 'width_fraction' from -1 to a width_fraction
      // from +1, in intervals of 2*width_fraction's. It is performed this way because each strand is to represent the
      // tractlet +- a width_fraction radius about it.
      for (double ax1_frac = (-1.0 + width_fraction); ax1_frac < 1.0; ax1_frac += 2.0 * width_fraction) {

        for (double ax2_frac = (-1.0 + width_fraction); ax2_frac < 1.0; ax2_frac += 2.0 * width_fraction) {

          //Include section if coord consisting of the combined fractions along each perturbation axes lies within the
          // unit circle minus a half a width_fraction. This will give the resulting cross-sections an elliptical
          // shape (or circular if axes are equal).
          if (MR::Math::pow2(ax1_frac) + MR::Math::pow2(ax2_frac) <= (1.0 - width_fraction /2.0)) {

            Strand strand (tractlet[0] + tractlet[1] * ax1_frac + tractlet[2] * ax2_frac);
            strand.copy_props(tractlet);

            push_back(strand);

            ++count;
          }

        }

        // Normalise the ACS property of the generated strands by the number generated
        if (has_elem_prop(ACS_SQRT_PROP))
          for (size_t strand_i = 0; strand_i < size(); ++strand_i)
            operator[](strand_i).set_acs(tractlet.acs() / (double)count);

      }

      freeze_elem_degree(); //Freeze element size

    }


    Track::Set                                  Strand::Set::to_tracks(size_t num_points) const {

      return Track::Set(*this, num_points);

    }


    Tractlet::Set                               Strand::Set::to_tractlets () const {

      return Tractlet::Set(*this);

    }



    Tractlet::Set                               Strand::Set::to_tractlets (double width) const {

      return Tractlet::Set(*this, width);

    }


    Tractlet::Set                               Strand::Set::to_tractlets(size_t num_tractlets, const std::string& cluster_save_location) const {

      return Tractlet::Set(*this, num_tractlets,cluster_save_location);

    }





    Strand                                      Strand::Set::centroid() const {

      Fibre::Strand centroid(0);

      for (size_t strand_i = 0; strand_i < size(); ++strand_i)
        centroid += operator[](strand_i);

      centroid /= (double)size();

      return centroid;

    }


    //TODO: fix up distance and flip functions to avoid unecessary copies of data.
    Strand::Set                                 Strand::Set::shortest_span() const {

      Strand::Set strands(*this);

      if (size()) {

        Fibre::Strand centroid;

        centroid = strands[0];

        size_t num_iterations = 0;

        double distance = -1.0;

        std::vector<size_t> flip_pattern(size());

        while (num_iterations < MAX_ITERATIONS) {

          double new_distance = 0;

          for (size_t strand_i = 0; strand_i < size(); ++strand_i)
            new_distance += strands[strand_i].distance(centroid,(bool&)flip_pattern[strand_i]);

          if (new_distance == distance)
            break;
          else
            distance = new_distance;

          centroid.zero();

          for (size_t strand_i = 0; strand_i < size(); ++strand_i)
            if (flip_pattern[strand_i])
              centroid += strands[strand_i].flip();
            else
              centroid += strands[strand_i];

          centroid /= (double)size();

          ++num_iterations;
        }

        if (num_iterations == MAX_ITERATIONS)
          throw Exception ("Could not find shortest span of strand set.");

        strands = strands.flip(flip_pattern);

      }

      return strands;

    }

    Strand::Set                                 Strand::Set::flip(const std::vector<size_t>& to_flip) const {

      if (to_flip.size() != size())
        throw Exception ("Size of 'to_flip' vector (" + str(to_flip.size()) + ") does not match that of strand set (" + str(size()) + ").");

      Strand::Set flipped (*this);

      for (size_t strand_i = 0; strand_i < size(); ++strand_i)
        if (to_flip[strand_i])
          flipped[strand_i] = operator[](strand_i).flip();

      return flipped;

    }




    Strand::Set                                 Strand::Set::peel(double distance, size_t num_length_sections, bool to_cube) const {

      Strand::Set include;

      include.set_extend_prop("peel_distance", str(distance));

      if (to_cube)
        include.set_extend_prop("peel_shape", "cube");
      else
        include.set_extend_prop("peel_shape", "sphere");


      include.add_extend_elem_props(*this);

      for (size_t strand_i = 0; strand_i < size(); ++strand_i) {

        Track tck = operator[](strand_i).to_track(num_length_sections);

        bool in_roi = false;

        for (size_t point_i = 0; point_i < tck.num_points(); point_i++) {

          if (to_cube) {
            if ((abs(tck[point_i][X]) <= distance) && (abs(tck[point_i][Y]) <= distance) && (abs(tck[point_i][Z]) <= distance)) {
              in_roi = true;
              break;
            }
          } else
            if (tck[point_i].norm() <= distance) {
              in_roi = true;
              break;
            }

        }

        if (in_roi)
          include.push_back(operator[](strand_i), get_extend_elem_prop_row(strand_i));

      }

      include.copy_props(*this);
      include.copy_relevant_elem_props(*this);

      return include;
    }


    Strand::Set&                                Strand::Set::cluster(size_t num_clusters) {

      if (num_clusters > size())
        throw Exception ("New number of clusters (" + str(num_clusters) + ") must be less than or equal to original number of strands (" + str(size()) + ")");


      KMterm  term( NUM_STAGES, 0, 0, 0,    // run for 100 stages
                    MAX_ACCUM_RDL,      // min consec RDL
                    MIN_ACCUM_RDL,      // min accum RDL
                    MAX_RUN_STAGES,     // max run stages
                    INIT_PROB_ACCEPTANCE,     // init. prob. of acceptance
                    TEMP_RUN_LENGTH,      // temp. run length
                    TEMP_REDUC_FACTOR);     // temp. reduction factor


      //-------------------------//
      //  k-means cluster strands
      //-------------------------//

      size_t degree = operator[](0).degree();

      int ptDim = degree * 3;   // dimension

      KMdata dataPts(ptDim, size());    // allocate data storage

      for (size_t strand_i = 0; strand_i < size(); strand_i++) {
        for (size_t degree_i = 0; degree_i < operator[](strand_i).degree(); degree_i++) {
          dataPts[strand_i][degree_i * 3 + X] = operator[](strand_i)[degree_i][X];
          dataPts[strand_i][degree_i * 3 + Y] = operator[](strand_i)[degree_i][Y];
          dataPts[strand_i][degree_i * 3 + Z] = operator[](strand_i)[degree_i][Z];
        }
      }



      dataPts.buildKcTree();      // build filtering structure

      KMfilterCenters ctrs(num_clusters, dataPts);   // allocate centers

      KMlocalHybrid kmHybrid(ctrs, term);   // Hybrid heuristic
      ctrs = kmHybrid.execute();

      //--------------------------//
      //  Get cluster information
      //--------------------------//


      KMctrIdxArray closeCtr = new KMctrIdx[dataPts.getNPts()];
      double* sqDist = new double[dataPts.getNPts()];
      ctrs.getAssignments(closeCtr, sqDist);

      if (!has_extend_elem_prop(Track::BUNDLE_INDEX_EPROP))
        add_extend_elem_prop(Track::BUNDLE_INDEX_EPROP, "-1");

      for (size_t strand_i = 0; strand_i < size(); ++strand_i)
        set_extend_elem_prop(Track::BUNDLE_INDEX_EPROP, str(closeCtr[strand_i]), strand_i);

      return *this;

    }





    Strand::Set                                 Strand::Set::reduce_number(size_t new_num_strands) {

      if (new_num_strands == this->size())
        return *this;

      if (new_num_strands < 1)
        throw Exception ("New number of strands must be greater than zero");

      if (new_num_strands > this->size())
        throw Exception ("New number of strands (" + str(new_num_strands) + ") must be less than original number of strands (" + str(this->size()) + ")");


      KMterm	term( NUM_STAGES, 0, 0, 0,		// run for 100 stages
                    MAX_ACCUM_RDL,			// min consec RDL
                    MIN_ACCUM_RDL,			// min accum RDL
                    MAX_RUN_STAGES,			// max run stages
                    INIT_PROB_ACCEPTANCE,			// init. prob. of acceptance
                    TEMP_RUN_LENGTH,			// temp. run length
                    TEMP_REDUC_FACTOR);			// temp. reduction factor


    //-------------------------//
    //  k-means cluster strands
    //-------------------------//

      size_t degree = operator[](0).degree();

      int	ptDim	= degree * 3;		// dimension

      KMdata dataPts(ptDim, size());		// allocate data storage

      for (size_t strand_i = 0; strand_i < size(); strand_i++) {
        for (size_t degree_i = 0; degree_i < operator[](strand_i).degree(); degree_i++) {
          dataPts[strand_i][degree_i * 3 + X] = operator[](strand_i)[degree_i][X];
          dataPts[strand_i][degree_i * 3 + Y] = operator[](strand_i)[degree_i][Y];
          dataPts[strand_i][degree_i * 3 + Z] = operator[](strand_i)[degree_i][Z];
        }
      }



      dataPts.buildKcTree();			// build filtering structure

      KMfilterCenters ctrs(new_num_strands, dataPts);		// allocate centers

      KMlocalHybrid kmHybrid(ctrs, term);		// Hybrid heuristic
      ctrs = kmHybrid.execute();

//    //--------------------------//
//    //  Get cluster information
//    //--------------------------//
//
//
//      KMctrIdxArray closeCtr = new KMctrIdx[dataPts.getNPts()];
//      double* sqDist = new double[dataPts.getNPts()];
//      ctrs.getAssignments(closeCtr, sqDist);


    //-----------------------------------//
    //  Save centre points as new strands
    //-----------------------------------//

      Strand::Set output_strands;

      for (size_t ctr_i = 0; ctr_i < new_num_strands; ctr_i++) {

        Strand strand(degree);

        for (size_t degree_i = 0; degree_i < operator[](0).degree(); degree_i++) {
          for (size_t dim_i = 0; dim_i < 3; dim_i++)
            strand[degree_i][dim_i] = ctrs[ctr_i][degree_i * 3 + dim_i];
        }

        output_strands.push_back(strand);

      }

      output_strands.copy_props(*this);
      output_strands.copy_relevant_elem_props(*this);

      return output_strands;

    }


    std::vector<size_t>                         Strand::Set::bundle(size_t bundle_index) const {

      if (!this->has_extend_elem_prop(Track::BUNDLE_INDEX_EPROP))
        throw Exception ("Bundle index property (" + Track::BUNDLE_INDEX_EPROP + ") not set on strand set.");

      std::vector<size_t> bundle;

      for (size_t strand_i = 0; strand_i < size(); strand_i++) {

        if (get_extend_elem_prop<size_t>(Track::BUNDLE_INDEX_EPROP, strand_i) == bundle_index)
          bundle.push_back(strand_i);

      }

      return bundle;

    }
    
    size_t                                      Strand::Set::max_degree() const {

      size_t max_degree = 0;

      for (size_t strand_i = 0; strand_i < size(); ++strand_i)
        if (operator[](strand_i).degree() > max_degree)
          max_degree = operator[](strand_i).degree();

      return max_degree;

    }


    Strand::Set                                 Strand::Set::bundle_backbones() const {
    

      if (!this->has_extend_elem_prop(Track::BUNDLE_INDEX_EPROP))
        throw Exception ("Bundle index property (" + Track::BUNDLE_INDEX_EPROP + ") not set on strand set.");
     

      Strand::Set backbones;
      backbones.add_extend_elem_prop(Track::BUNDLE_INDEX_EPROP,"-1");

      std::vector<size_t> bundle_indices = this->bundle_indices();

      for (size_t bundle_i = 0; bundle_i < bundle_indices.size(); ++bundle_i) {

        std::vector<size_t> bundle = this->bundle(bundle_indices[bundle_i]);

        Fibre::Strand backbone;

        for (size_t strand_i = 0; strand_i < bundle.size(); ++strand_i)
          backbone += operator[](bundle[strand_i]);

        backbone /= (double)bundle.size();

        std::map<std::string,std::string> properties_row;
        properties_row[Track::BUNDLE_INDEX_EPROP] = str(operator[](bundle_indices[bundle_i]));

        backbones.push_back(backbone, properties_row);


      }

      backbones.copy_props(*this);
      backbones.copy_relevant_elem_props(*this);

      return backbones;
      
    }
    
    
    std::vector<size_t>                         Strand::Set::bundle_indices() const {
    
      if (!this->has_extend_elem_prop(Track::BUNDLE_INDEX_EPROP))
        throw Exception ("Bundle index property (" + Track::BUNDLE_INDEX_EPROP + ") not set on strand set.");
    
      std::vector<size_t> bundle_indices;
      
      for (size_t strand_i = 0; strand_i < size(); ++strand_i)  {
        size_t bundle_i = get_extend_elem_prop<size_t>(Track::BUNDLE_INDEX_EPROP, strand_i);

        if (find(bundle_indices.begin(), bundle_indices.end(), bundle_i) == bundle_indices.end())
          bundle_indices.push_back(bundle_i);
          
      }
          
      return bundle_indices;
    
    }


    void                                        Strand::Set::project_to_sphere(double sphere_radius) {

      for (size_t strand_i = 0; strand_i < size(); ++strand_i)
        operator[](strand_i).project_to_sphere(sphere_radius);

    }



    void                                        Strand::Set::load(const std::string& location, size_t degree, double strands_per_acs) {

      this->clear();

      if (File::has_or_txt_extension<Strand>(location)) {

        Base::Set<Strand>::load(location);

      } else if (File::has_or_txt_extension<Track>(location)) {

        if (degree < 1)
          throw Exception(" '-degree' parameter is required for implicit conversion from track set ('" + location + "').");

        Track::Set tcks (location);
        
        operator=(tcks.to_strands(degree));

        this->add_extend_elem_prop(ORIGINAL_NUM_POINTS_PROP, "0");
        
        for (size_t tck_i = 0; tck_i < tcks.size(); tck_i++)
          this->set_extend_elem_prop(ORIGINAL_NUM_POINTS_PROP, str(tcks[tck_i].num_points()), tck_i);

      } else if (File::has_or_txt_extension<Tractlet>(location)) {

        if (strands_per_acs <= 0.0)
          throw Exception(" '-strands_per_acs' parameter is required for implicit conversion from tractlet set ('" + location + "').");

        Tractlet::Set tractlets (location);

        set_extend_prop(ORIGINAL_NUM_TRACTLETS_PROP, str(tractlets.size()));

        operator=(tractlets.to_strands(strands_per_acs));

      } else
         throw Exception("The extension of file \"" + location + "\" is not a recognised type (\"" + Strand::FILE_EXTENSION + "\" or \"" + Track::FILE_EXTENSION + "\").");

    }


    Strand::Set                                 Strand::Set::permute(const std::vector<size_t>& indices) const {

      if (var_elem_degrees())
        throw Exception ("'permute' is currently not implemented for variable element size sets.");

      Set permuted_set (size(), degree(), *props, *elem_props);

      this->Fibre::Base::Set<Strand>::permute(permuted_set, indices);

      return permuted_set;

    }


    void                                        Strand::Set::save(const std::string& location, size_t num_points, size_t num_tractlets) const {


      if (File::has_or_txt_extension<Strand>(location)) {

        //Remove ORIGINAL_NUM_POINTS_PROP property if it exists as it is no longer relevant.
        if (this->has_extend_elem_prop(ORIGINAL_NUM_POINTS_PROP)) {
          
          Set strands = *this;
          
          strands.remove_extend_elem_prop(ORIGINAL_NUM_POINTS_PROP);
          
          strands.save(location);
        
        } else
        
          Base::Set<Strand>::save(location);
      
      } else if (File::has_or_txt_extension<Track>(location)) {

        Track::Set tcks;
      
        if (num_points != 0) 
        
          tcks = this->to_tracks(num_points);
          
        else if (this->has_extend_elem_prop(ORIGINAL_NUM_POINTS_PROP)) {

          for (size_t strand_i = 0; strand_i < this->size(); strand_i++)  {
            
            size_t loaded_num_points = get_extend_elem_prop<size_t>(ORIGINAL_NUM_POINTS_PROP, strand_i);
            tcks.push_back(this->operator[](strand_i).to_track(loaded_num_points), this->get_extend_elem_prop_row(strand_i));
            
          }
        
        } else
          throw Exception(" \"-num_points\" parameter is required for implicit conversion to track set ('" + location + "').");   

        tcks.remove_extend_elem_prop(ORIGINAL_NUM_POINTS_PROP);

        tcks.save(location);
      
      } else if (File::has_or_txt_extension<Tractlet>(location)) {

        Tractlet::Set tcts;

        std::map<std::string,std::string> properties = get_extend_props();

        std::map<std::string, std::string>::const_iterator num_tract_it = properties.find(ORIGINAL_NUM_TRACTLETS_PROP);

        if (num_tractlets)
          tcts = this->to_tractlets(num_tractlets, File::strip_extension(location) + ".clustered.str");
        else if (num_tract_it != properties.end())
          tcts = this->to_tractlets(to<size_t>(num_tract_it->second));
        else if (has_extend_elem_prop(Track::BUNDLE_INDEX_EPROP))
          tcts = this->to_tractlets();
        else
          throw Exception(" '-num_tractlets' is required for implicit conversion to tractlet set ('" + location + "') if '" + Track::BUNDLE_INDEX_EPROP + "' is not present.");

        tcts.save(location);
      
      } else 
         throw Exception("The extension of file \"" + location + "\" is not a recognised type (\"" + Strand::FILE_EXTENSION + "\" or \"" + Track::FILE_EXTENSION + "\").");      
                
    }  


    std::string                                 Strand::Set::matlab_str() const {

      std::stringstream ss;

      ss << "{ 'b_intens', " << base_intensity() << ", ";

      for (size_t strand_i = 0; strand_i < size(); strand_i++) {

        ss << "{ 'acs', " << operator[](strand_i).acs() << ", [";

        for (size_t degree_i = 0; degree_i < operator[](strand_i).size(); degree_i++) {

          for (size_t dim_i = 0; dim_i < 3; dim_i++)
            ss << operator[](strand_i)[degree_i][dim_i] << " ";

          if (degree_i != 2)
            ss << "; ";

        }

        ss << "] }";

        if (strand_i != (size() -1))
          ss << ", ";

      }

      ss << " }";

      return ss.str();

    }

  }

}
