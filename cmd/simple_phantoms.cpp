/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Written by Thomas G. Close, 04/03/2009.

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

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include "bts/cmd.h"

#include "bts/common.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"

#include "bts/mcmc/proposal/distribution.h"
#include "bts/mcmc/proposal/distribution/gaussian.h"

#include "bts/image/buffer.h"
#include "bts/image/expected/buffer.h"
#include "bts/diffusion/model.h"

#include "bts/inline_functions.h"

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Generates a particular set of basic phantom configurations",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("", "").optional(),
    Argument ("", "").optional(),
    Argument()
};

OPTIONS= {

    Option ("type", "Type of the basic configuration to generate. (0 - crossing, 1 - branching, 2 - kissing).")
    + Argument ("type", "").type_integer (0, 2, 0),

    Option ("degree", "degree of the strands.")
    + Argument ("degree", "").type_integer (1, 2000, 100),

    Option ("curve", "curve of the bends.")
    + Argument ("curve", "").type_float (1e-9, INFINITY, 0.1),

    Option ("num_layers", "number of layers.")
    + Argument ("num_layers", "").type_integer (1, 2000, 100),

    Option ("num_voxels", "number of voxels.")
    + Argument ("num_voxels", "").type_integer (1, 2000, 100),

    Option ("voxel_size", "Size of the voxels.")
    + Argument ("voxel_size", "").type_float (1e-9, INFINITY, 0.1),

    Option ("angle", "Angle of the branch (in degrees)")
    + Argument ("angle", "").type_float (1e-9, INFINITY, 0.1),

    Option ("curv_angle", "Angle of the branch (in degrees)")
    + Argument ("curv_angle", "").type_float (1e-9, INFINITY, 0.1),

    Option()};

EXECUTE {
    
        const size_t CROSSING = 0;
        const size_t BRANCHING = 1;
        const size_t KISSING = 2;
        
        size_t type = CROSSING;
        
        size_t degree = 3;
        double curve = 0.1;
        size_t num_layers = 6;
        double vox_size = 0.15;
        size_t num_voxels = 3;
        double angle = 30;
        double curv_angle = 30;
        
        Options opt = get_options("type");
        if (opt.size())
            type = opt[0][0];
        
        opt = get_options("degree");
        if (opt.size())
            degree = opt[0][0];
        
        opt = get_options("curve");
        if (opt.size())
            curve = opt[0][0];
        
        opt = get_options("num_layers");
        if (opt.size())
            num_layers = opt[0][0];
        
        opt = get_options("num_voxels");
        if (opt.size())
            num_voxels = opt[0][0];
        
        opt = get_options("vox_size");
        if (opt.size())
            vox_size = opt[0][0];
        
        opt = get_options("angle");
        if (opt.size())
            angle = opt[0][0];
        
        opt = get_options("curv_angle");
        if (opt.size())
            curv_angle = opt[0][0];
        
        switch (type) {
            
            case CROSSING: {
                
                double extent = vox_size * (double) num_voxels / 2.0;
                
                double strand_diameter = 2.0 * extent / (double) num_layers;
                
                Fibre::Strand straight(degree, 0.0);
                
                straight[1][1] = M_SQRT2 * extent;
                
                Fibre::Strand curved(degree, 0.0);
                
                curved[1][0] = M_SQRT2 * extent;
                curved[2][1] = curve;
                curved[0][1] = ((double) curve);
                
                Fibre::Strand::Set ismrm;
                
                ismrm.add_extend_elem_prop("bundle_index", "0");
                
                double layer_offset = -extent + strand_diameter / 2.0;
                double col_offset = -extent + strand_diameter / 2.0;
                
                for (size_t layer_i = 0; layer_i < num_layers / 2; layer_i++) {
                    
                    for (size_t col_i = 0; col_i < num_layers; col_i++) {
                        
                        Fibre::Strand straight_offset = straight;
                        Fibre::Strand curved_offset = curved;
                        
                        straight_offset[0][0] += strand_diameter * (double) col_i + col_offset;
                        straight_offset[0][2] += 2.0 * strand_diameter * (double) layer_i
                                + layer_offset;
                        
                        curved_offset[0][1] += strand_diameter * (double) col_i + col_offset;
                        curved_offset[0][2] += 2.0 * strand_diameter * (double) layer_i
                                + layer_offset + strand_diameter;
                        
                        ismrm.push_back(straight_offset);
                        ismrm.set_extend_elem_prop("bundle_index", "0", ismrm.size() - 1);
                        
                        ismrm.push_back(curved_offset);
                        ismrm.set_extend_elem_prop("bundle_index", "1", ismrm.size() - 1);
                        
                    }
                    
                }
                
                ismrm.add_extend_elem_prop("track_radius", str(strand_diameter / 2.0));
                
                ismrm.save(PARAM_DIRECTORY + "ismrm_abstractlet/ismrm-crossing.frr");
                
                break;
                
            }
            case BRANCHING: {
                
                double extent = vox_size * (double) num_voxels / 2.0;
                
                double strand_diameter = 2.0 * extent / (double) num_layers;
                
                Fibre::Strand curve1(degree, 0.0);
                
                double rad_angle = 2.0 * M_PI * (angle / 360.0);
                
                double curve_rad_angle = 2.0 * M_PI * (curv_angle / 360.0);
                
                curve1[0][0] = -curve * cos(curve_rad_angle);
                curve1[0][1] = 2 * curve * sin(curve_rad_angle);
                
                curve1[1][0] = M_SQRT2 * extent * sin(rad_angle) * 2;
                curve1[1][1] = M_SQRT2 * extent * cos(rad_angle) * 2;
                
                curve1[2][0] = curve * cos(curve_rad_angle);
                curve1[2][1] = curve * sin(curve_rad_angle);
                
                Fibre::Strand curve2 = curve1;
                for (size_t degree_i = 0; degree_i < degree; degree_i++) {
                    curve2[degree_i][1] *= -1;
                    
                }
                
                std::cout << curve1 << std::endl;
                std::cout << curve2 << std::endl;
                
                Fibre::Strand::Set ismrm;
                
                ismrm.add_extend_elem_prop("bundle_index", "0");
                
                double layer_offset = -extent + strand_diameter / 2.0;
                double col_offset = -extent + strand_diameter / 2.0;
                
                for (size_t layer_i = 0; layer_i < num_layers / 2; layer_i++) {
                    
                    for (size_t col_i = 0; col_i < num_layers; col_i++) {
                        
                        Fibre::Strand curve2_offset = curve2;
                        Fibre::Strand curve1_offset = curve1;
                        
                        curve1_offset[0][0] += (strand_diameter * (double) col_i + col_offset)
                                * sin(rad_angle);
                        curve1_offset[0][1] += -(strand_diameter * (double) col_i + col_offset)
                                * cos(rad_angle);
                        
                        curve1_offset[0][2] += 2.0 * strand_diameter * (double) layer_i
                                + layer_offset + strand_diameter;
                        
                        curve2_offset[0][0] += (strand_diameter * (double) col_i + col_offset)
                                * sin(rad_angle);
                        curve2_offset[0][1] += (strand_diameter * (double) col_i + col_offset)
                                * cos(rad_angle);
                        
                        curve2_offset[0][2] += 2.0 * strand_diameter * (double) layer_i
                                + layer_offset;
                        
                        ismrm.push_back(curve2_offset);
                        ismrm.set_extend_elem_prop("bundle_index", "0", ismrm.size() - 1);
                        
                        ismrm.push_back(curve1_offset);
                        ismrm.set_extend_elem_prop("bundle_index", "1", ismrm.size() - 1);
                        
                    }
                    
                }
                
                ismrm.add_extend_elem_prop("track_radius", str(strand_diameter / 2.0));
                
                ismrm.save(PARAM_DIRECTORY + "ismrm_abstractlet/ismrm-branching.frr");
                break;
                
            }
            case KISSING: {
                
                break;
            }
                
        }
        
    }
    
