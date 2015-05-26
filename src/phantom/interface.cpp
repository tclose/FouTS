/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Fourier Tract Sampling (FouTS).

 FouTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FouTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#include <cstdlib>

//#include "bts/cmd.h"

#include "bts/common.h"

#include "phantom/interface.h"

#include "bts/fibre/strand.h"
#include "bts/fibre/track.h"
#include "bts/fibre/track/set.h"

#include "bts/inline_functions.h"

//extern "C" {
#include "phantom/shared/strand_collection.h"
//}

std::vector<FTS::Triple<double> >& generate_pre_points(
        FTS::Fibre::Track::Set& tracks, std::vector<FTS::Triple<double> >& pre_points) {
    
    for (size_t tck_i = 0; tck_i < tracks.size(); ++tck_i) {
        
        FTS::Fibre::Track tck = tracks[tck_i];
        
        if (tck.size() > 1) {
            FTS::Triple<double> pre_point = tck[0] - (tck[1] - tck[0]);
            pre_points.push_back(pre_point);
        } else {
            std::cout << "Warning strand " << tck_i
                      << " has been ommitted as it has less than two control points (" << tck.size()
                      << ")" << std::endl;
        }
        
    }
    
    return pre_points;
    
}

std::vector<FTS::Triple<double> >& generate_post_points(
        FTS::Fibre::Track::Set& tracks, std::vector<FTS::Triple<double> >& post_points) {
    
    for (size_t tck_i = 0; tck_i < tracks.size(); ++tck_i) {
        
        FTS::Fibre::Track tck = tracks[tck_i];
        
        if (tck.size() > 1) {
            FTS::Triple<double> post_point = tck[tck.size() - 1]
                    + (tck[tck.size() - 1] - tck[tck.size() - 2]);
            post_points.push_back(post_point);
        } else {
            std::cout << "Warning strand " << tck_i
                      << " has been ommitted as it has less than two control points (" << tck.size()
                      << ")" << std::endl;
        }
        
    }
    
    return post_points;
}

Strand_collection* convert_mr_to_nfg(Strand_collection* c, FTS::Fibre::Track::Set& tracks,
                                     std::vector<FTS::Triple<double> >& pre_points,
                                     std::vector<FTS::Triple<double> >& post_points) {
    
    if (tracks.size() == 0)
        throw FTS::Exception(
                "Track set supplied to FTSTrix to NFG conversion function does not have any tracks in it.");
    
    if (!tracks.has_extend_elem_prop(FTS::Fibre::Track::RADIUS_PROP))
        tracks.add_extend_elem_prop(FTS::Fibre::Track::RADIUS_PROP, "0.03");
    
    size_t total_num_control_points = 0;
    
    for (size_t tck_i = 0; tck_i < tracks.size(); ++tck_i)
        total_num_control_points += tracks[tck_i].size() - 2;
    
    collection_alloc(c, tracks.size(), total_num_control_points, 0);
    
    size_t control_point_count = 0;
    
    for (size_t strand_i = 0; strand_i < tracks.size(); strand_i++) {
        
        FTS::Fibre::Track tck = tracks[strand_i];
        
        c->num_strand_control_points[strand_i] = tck.size() - 2;    //Start and end points are not included.
                
        c->strand_r[strand_i] = atof(
                tracks.get_extend_elem_prop(FTS::Fibre::Track::RADIUS_PROP, strand_i).c_str());
        
        size_t bundle_i;
        if (tracks.has_extend_elem_prop(FTS::Fibre::Track::BUNDLE_INDEX_EPROP))
            bundle_i =
                    atoi(
                            tracks.get_extend_elem_prop(FTS::Fibre::Track::BUNDLE_INDEX_EPROP,
                                    strand_i).c_str());
        else
            bundle_i = strand_i;
        
        c->bundle_i_of_strand[strand_i] = bundle_i;
        
        pre_points[strand_i].copy_to(&(c->pre_points[strand_i * 3]));
        post_points[strand_i].copy_to(&(c->post_points[strand_i * 3]));
        
        tck[0].copy_to(&(c->start_points[strand_i * 3]));
        
        size_t point_i;
        
        for (point_i = 1; point_i < (tck.size() - 1); point_i++)
            tck[point_i].copy_to(&(c->control_points[(control_point_count + point_i - 1) * 3]));
        
        tck[point_i].copy_to(&(c->end_points[strand_i * 3]));
        
        construct_strand(&(c->strands[strand_i]), strand_i, bundle_i,
                &(c->control_points[control_point_count * 3]), &(c->start_points[strand_i * 3]),
                &(c->end_points[strand_i * 3]), &(c->pre_points[strand_i * 3]),
                &(c->post_points[strand_i * 3]), &(c->segments[control_point_count + 3 * strand_i]),
                c->num_strand_control_points[strand_i], 0.0, c->strand_r[strand_i]);
        
        control_point_count += c->num_strand_control_points[strand_i];
        
    }
    
    c->num_bundles = construct_bundles(c->bundles, c->bundle_i_of_strand, c->num_strands,
            c->strands);
    
    return c;
}

FTS::Fibre::Track::Set& convert_nfg_to_mr(FTS::Fibre::Track::Set& tracks,
                                          std::vector<FTS::Triple<double> >& pre_points,
                                          std::vector<FTS::Triple<double> >& post_points,
                                          Strand_collection* c) {
    
    size_t total_num_control_points = 0;
    
    tracks.clear();
    pre_points.clear();
    post_points.clear();
    
    tracks.add_extend_elem_prop(FTS::Fibre::Track::BUNDLE_INDEX_EPROP, "NaN");
    tracks.add_extend_elem_prop(FTS::Fibre::Track::RADIUS_PROP, "NaN");
    
    for (int strand_i = 0; strand_i < c->num_strands; strand_i++) {
        
        FTS::Fibre::Track tck;
        
        pre_points.push_back(
                FTS::Triple<double>((double) c->pre_points[strand_i * 3 + X],
                        (double) c->pre_points[strand_i * 3 + Y],
                        (double) c->pre_points[strand_i * 3 + Z]));
        post_points.push_back(
                FTS::Triple<double>((double) c->post_points[strand_i * 3 + X],
                        (double) c->post_points[strand_i * 3 + Y],
                        (double) c->post_points[strand_i * 3 + Z]));
        
        tck.push_back(
                FTS::Coord((double) c->start_points[strand_i * 3 + X],
                        (double) c->start_points[strand_i * 3 + Y],
                        (double) c->start_points[strand_i * 3 + Z]));
        
        for (int point_i = 0; point_i < c->num_strand_control_points[strand_i]; point_i++) {
            
            tck.push_back(
                    FTS::Coord((double) c->control_points[total_num_control_points * 3 + X],
                            (double) c->control_points[total_num_control_points * 3 + Y],
                            (double) c->control_points[total_num_control_points * 3 + Z]));
            
            total_num_control_points++;
            
        }
        
        tck.push_back(
                FTS::Coord((double) c->end_points[strand_i * 3 + X],
                        (double) c->end_points[strand_i * 3 + Y],
                        (double) c->end_points[strand_i * 3 + Z]));
        
        std::map<std::string, std::string> properties_row;
        properties_row[FTS::Fibre::Track::BUNDLE_INDEX_EPROP] = FTS::str(
                c->bundle_i_of_strand[strand_i]);
        properties_row[FTS::Fibre::Track::RADIUS_PROP] = FTS::str(c->strand_r[strand_i]);
        
        tracks.push_back(tck, properties_row);
        
    }
    
    return tracks;
    
}

