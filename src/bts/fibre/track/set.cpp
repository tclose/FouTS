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

#include "bts/fibre/track/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/track.h"

#include "bts/fibre/base/reader.cpp.h"
#include "bts/fibre/base/writer.cpp.h"

namespace BTS {
    
    namespace Fibre {
        
        const std::string Track::Set::FILE_EXTENSION = "kst";
        
        const char* Track::Set::PROPS_LIST[] = { Base::Set<Track>::BASE_INTENSITY_PROP,
                                                 PROPS_LIST_END };
        
        const std::string Track::Set::ORIGINAL_DEGREE_PROP = "loaded_degree";
        
        Track::Set::Set(const Strand::Set& strands, size_t num_length_sections)
                : Base::Set<Track>(strands.size(), num_length_sections,
                        num_length_sections * 3 + select_props<Strand>(*strands.elem_props).size(),
                        select_props<Set>(*strands.props),
                        select_props<Strand>(*strands.elem_props), strands.get_extend_props()) {
            
            from_strands(strands, num_length_sections);
            
        }
        
        Track::Set::Set(const Tractlet::Set& tractlets, size_t num_length_sections,
                        size_t num_width_sections)
                : Base::Set<Track>(tractlets.size(), num_length_sections,
                        num_length_sections * 3 + select_props<Strand>(*tractlets.elem_props).size(),
                        select_props<Set>(*tractlets.props),
                        select_props<Strand>(*tractlets.elem_props), tractlets.get_extend_props()) {
            
            Strand::Set strands(tractlets, num_width_sections);
            
            resize(strands.size());
            
            from_strands(strands, num_length_sections);
            
        }
        
        Track::Set::Set(const Tractlet::Set& tractlets, size_t num_length_sections,
                        double strands_per_acs)
                : Base::Set<Track>(tractlets.size(), num_length_sections,
                        num_length_sections * 3 + select_props<Strand>(*tractlets.elem_props).size(),
                        select_props<Set>(*tractlets.props),
                        select_props<Strand>(*tractlets.elem_props), tractlets.get_extend_props()) {
            
            Strand::Set strands(tractlets, strands_per_acs);
            
            resize(strands.size());
            
            from_strands(strands, num_length_sections);
            
        }
        
        void Track::Set::from_strands(const Strand::Set& strands, size_t num_length_sections) {
            
            // Record the original degree of the strands, just in case the set needs to be converted back into strands again
            // afterwards.
            add_extend_elem_prop(Track::Set::ORIGINAL_DEGREE_PROP, "0");
            
            for (size_t strand_i = 0; strand_i < size(); strand_i++) {
                
                operator[](strand_i) = Track(strands[strand_i], num_length_sections);
                set_extend_elem_prop(Track::Set::ORIGINAL_DEGREE_PROP,
                        str(strands[strand_i].degree()), strand_i);
                
            }
            
            for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
                prop(prop_i) = strands.prop(prop_key(prop_i));
            
            copy_extend_elem_props(strands);
            
        }
        
        Strand::Set Track::Set::to_strands(size_t degree) const {
            
            return Strand::Set(*this, degree);
            
        }
        
        Tractlet::Set Track::Set::to_tractlets(size_t degree) const {
            
            return Tractlet::Set(*this, degree);
            
        }
        
        Tractlet::Set Track::Set::to_tractlets(size_t degree, size_t num_tractlets) const {
            
            return Tractlet::Set(*this, degree, num_tractlets);
            
        }
        
        Tractlet::Set Track::Set::to_tractlets(size_t degree, double width) const {
            
            return Tractlet::Set(*this, degree, width);
            
        }
        
        void Track::Set::remove_short_tracks(size_t min_num_points) {
            
            for (size_t tck_i = 0; tck_i < this->size(); tck_i++) {
                
                if (operator[](tck_i).num_points() < min_num_points) {
                    
                    std::cout << std::endl << "Erased tck " << tck_i << " as it only had "
                              << operator[](tck_i).num_points() << " control points." << std::endl;
                    
                    this->erase(tck_i);
                    
                }
                
            }
            
        }
        
        void Track::Set::load(const std::string& location, size_t num_points,
                              double strands_per_acs) {
            
            this->clear();
            
            if (File::has_or_txt_extension<Track>(location))
                
                this->Base::Set<Track>::load(location);
            
            else if (File::has_or_txt_extension<Strand>(location)) {
                
                if (num_points < 1)
                    throw Exception(
                            " '-num_points' parameter is required for implicit conversion from strand set ('" + location
                            + "').");
                
                Strand::Set strands(location);
                
                operator=(strands.to_tracks(num_points));
                
                this->add_extend_elem_prop(ORIGINAL_DEGREE_PROP, "0");
                
                for (size_t strand_i = 0; strand_i < strands.size(); strand_i++)
                    this->set_extend_elem_prop(ORIGINAL_DEGREE_PROP,
                            str(strands[strand_i].degree()), strand_i);
                
            } else if (File::has_or_txt_extension<Tractlet>(location)) {
                
                if (strands_per_acs <= 0.0)
                    throw Exception(
                            " '-strands_per_acs' parameter is required for implicit conversion from tractlet set ('" + location
                            + "').");
                
                Tractlet::Set tractlets(location);
                
                operator=(tractlets.to_tracks(num_points, strands_per_acs));
                
                this->add_extend_elem_prop(ORIGINAL_DEGREE_PROP, "0");
                
                for (size_t tck_i = 0; tck_i < size(); tck_i++) {
                    
                    size_t tract_i = to<size_t>(get_extend_elem_prop(BUNDLE_INDEX_EPROP, tck_i));
                    
                    set_extend_elem_prop(ORIGINAL_DEGREE_PROP, str(tractlets[tract_i].degree()),
                            tck_i);
                }
                
            } else
                throw Exception(
                        "The extension of file \"" + location + "\" is not a recognised type (\""
                        + Track::FILE_EXTENSION + "\" or \"" + Strand::FILE_EXTENSION + "\").");
            
        }
        
        void Track::Set::save(const std::string& location, size_t degree,
                              size_t num_tractlets) const {
            
            if (File::has_or_txt_extension<Track>(location)) {
                
                //Remove ORIGINAL_DEGREE_PROP property if it exists as it is no longer relevant.
                if (this->has_extend_elem_prop(ORIGINAL_DEGREE_PROP)) {
                    
                    Set tcks = *this;
                    
                    tcks.remove_extend_elem_prop(ORIGINAL_DEGREE_PROP);
                    
                    tcks.Base::Set<Track>::save(location);
                    
                } else
                    
                    this->Base::Set<Track>::save(location);
                
            } else if (File::has_or_txt_extension<Strand>(location)) {
                
                Strand::Set strands;
                
                if (degree)
                    strands = this->to_strands(degree);
                
                else if (this->has_extend_elem_prop(ORIGINAL_DEGREE_PROP)) {
                    
                    for (size_t tck_i = 0; tck_i < this->size(); tck_i++) {
                        size_t loaded_degree = to<int>(
                                this->get_extend_elem_prop(ORIGINAL_DEGREE_PROP, tck_i));
                        strands.push_back(this->operator[](tck_i).to_strand(loaded_degree));
                    }
                    
                    strands.copy_extend_elem_props(*this);
                    
                } else
                    throw Exception(
                            " '-degree' is required for implicit conversion to strand set ('" + location
                            + "').");
                
                strands.remove_extend_elem_prop(ORIGINAL_DEGREE_PROP);
                
                strands.save(location);
                
            } else if (File::has_or_txt_extension<Tractlet>(location)) {
                
                Tractlet::Set tractlets;
                
                Strand::Set strands;
                
                if (degree)
                    strands = this->to_strands(degree);
                else if (this->has_extend_elem_prop(ORIGINAL_DEGREE_PROP)) {
                    for (size_t tck_i = 0; tck_i < this->size(); tck_i++) {
                        size_t loaded_degree = to<int>(
                                this->get_extend_elem_prop(ORIGINAL_DEGREE_PROP, tck_i));
                        strands.push_back(this->operator[](tck_i).to_strand(loaded_degree));
                    }
                    strands.copy_extend_elem_props(*this);
                } else
                    throw Exception(
                            " '-degree' is required for implicit conversion to strand set ('" + location
                            + "').");
                
                if (!has_extend_elem_prop(BUNDLE_INDEX_EPROP) && !num_tractlets)
                    throw Exception(
                            " '-num_tractlets' is required for implicit conversion to tractlet set ('" + location
                            + "') if '" + BUNDLE_INDEX_EPROP + "' is not present.");
                
                tractlets = strands.to_tractlets(num_tractlets);
                tractlets.copy_extend_elem_props(*this);
                tractlets.save(location);
                
            } else
                throw Exception(
                        "The extension of file \"" + location + "\" is not a recognised type (\""
                        + Track::FILE_EXTENSION + "\" or \"" + Strand::FILE_EXTENSION + "\").");
            
        }
        
        Track::Set Track::Set::peel(double distance, bool to_cube, bool by_bundles) const {
            
            Track::Set include(get_extend_props());
            
            include.set_extend_prop("peel_distance", str(distance));
            
            if (to_cube)
                include.set_extend_prop("peel_shape", "cube");
            else
                include.set_extend_prop("peel_shape", "sphere");
            
            include.add_extend_elem_props(*this);
            
            std::vector<bool> in_rois(max_bundle_index());
            for (size_t i = 0; i < in_rois.size(); ++i)
                in_rois[i] = false;

            for (size_t tck_i = 0; tck_i < this->size(); tck_i++) {
                
                Track tck = operator[](tck_i);
                size_t in_roi_i;
                if (has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP) && by_bundles)
                    in_roi_i = to<size_t>(get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, tck_i));
                else
                    in_roi_i = tck_i;
                for (size_t point_i = 0; point_i < tck.num_points(); point_i++) {
                    
                    if (to_cube) {
                        if ((abs(tck[point_i][X]) <= distance) && (abs(tck[point_i][Y]) <= distance)
                            && (abs(tck[point_i][Z]) <= distance)) {
                            in_rois[in_roi_i] = true;
                            break;
                        }
                    } else if (tck[point_i].norm() <= distance) {
                        in_rois[in_roi_i] = true;
                        break;
                    }
                    
                }
            }
            
            for (size_t tck_i = 0; tck_i < this->size(); ++tck_i) {
                size_t in_roi_i;
                if (has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP) && by_bundles)
                    in_roi_i = to<size_t>(get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, tck_i));
                else
                    in_roi_i = tck_i;
                if (in_rois[in_roi_i])
                    include.push_back(operator[](tck_i), this->get_extend_elem_prop_row(tck_i));
            }

            include.copy_props(*this);
            
            return include;
        }

        size_t Track::Set::max_bundle_index() const {
            size_t max_bundle_index = 0;
            if (has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP)) {
                for (size_t tck_i = 0; tck_i < size(); ++tck_i) {
                    size_t bundle_index = to<size_t>(
                            get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, tck_i));
                    if (bundle_index > max_bundle_index)
                        max_bundle_index = bundle_index;
                }
            } else
                max_bundle_index = size();
            return max_bundle_index;
        }
    
    }

}
