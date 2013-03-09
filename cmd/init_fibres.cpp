/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 18/06/2010.

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

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include "progressbar.h"
#include "bts/cmd.h"

#include "bts/common.h"

#include "bts/image/expected/buffer.h"

#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"

#include "bts/inline_functions.h"

using namespace BTS;

const double BASE_INTENSITY_DEFAULT = 0.0;
const size_t DEGREE_DEFAULT = 3;
const size_t NUM_FIBRES_DEFAULT = 5;
//const float   LENGTH_STDDEV_DEFAULT = 0.02;
//const float   WIDTH_MEAN_DEFAULT= 0.05;
//const float   WIDTH_STDDEV_DEFAULT= 0.005;
const float LENGTH_STDDEV_DEFAULT = 0.3;
const float WIDTH_MEAN_DEFAULT = 0.75;
const float WIDTH_STDDEV_DEFAULT = 0.075;
const float REJECT_RADIUS_DEFAULT = NAN;
const float DEFAULT_REJECT_RADIUS_SCALE = 1.5;
const double EDGE_BUFFER_DEFAULT = 0.075;
const size_t MAX_NUM_ATTEMPTS = 10000;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Initialises Fourier tracts within the image bounds with uniformly distributed orientations, normally distributed "
    "lengths and widths for the initialisation of Bayesian sampling. The default values produce small tracts to allow "
    "the tracts to ""grow"" in the direction of the underlying fibre orientation.",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("output","Output file location"),
    Argument()
};

OPTIONS= {

    Option("acs", "Intensity of the generated fibres.")
    + Argument("acs","").type_float(SMALL_FLOAT,NAN,LARGE_FLOAT),

    Option("base_intensity", "Base intensity of the generated fibres.")
    + Argument("base_intensity","").type_float(SMALL_FLOAT,BASE_INTENSITY_DEFAULT,LARGE_FLOAT),

    Option("degree", "The degree of the generated fibres.")
    + Argument("degree","").type_integer(0,DEGREE_DEFAULT,LARGE_INT),

    Option("num_fibres", "Size of the generated fibre sets.")
    + Argument("num_fibres","").type_integer(0,NUM_FIBRES_DEFAULT,LARGE_INT),

    Option("length_stddev", "Standard deviation of the projection along each (i.e. x,y,z) axes.")
    + Argument("length_stddev","").type_float(SMALL_FLOAT, LENGTH_STDDEV_DEFAULT,LARGE_FLOAT),

    Option("width_mean", "Mean width of the initialised fibres.")
    + Argument("width_mean","").type_float(SMALL_FLOAT, WIDTH_MEAN_DEFAULT,LARGE_FLOAT),

    Option("width_stddev", "Standard deviation of width of the initialised fibres")
    + Argument("width_stddev","").type_float(SMALL_FLOAT, WIDTH_STDDEV_DEFAULT,LARGE_FLOAT),

    Option("reject_radius", "If a centre of a randomly generated fibre lies withinness the radius of a previously "
            "generated strand it is rejected.")
    + Argument("reject_radius","").type_float(SMALL_FLOAT, REJECT_RADIUS_DEFAULT,LARGE_FLOAT),

    Option("edge_buffer", "Restrict the distributed_centres of the generated fibres to be withinness a certain distance "
            "from the edge of the ROI.")
    + Argument("edge_buffer","").type_float(-SMALL_FLOAT, EDGE_BUFFER_DEFAULT,LARGE_FLOAT),

    Option("curve_stddev", "The standard deviation of the curvature parameters.")
    + Argument("curve_stddev", "").type_float(0, 0, LARGE_FLOAT),

    Option("width_epsilon", "The width epsilon of the generated fibres.")
    + Argument("width_epsilon", "").type_float(0.0, Fibre::Tractlet::Set::WIDTH_EPSILON_DEFAULT, LARGE_FLOAT),

    Option("length_epsilon", "The length epsilon of the generated fibres.")
    + Argument("length_epsilon", "").type_float(0.0, Fibre::Tractlet::Set::LENGTH_EPSILON_DEFAULT, LARGE_FLOAT),

    Option ("seed", "Seed for the random generation")
    + Argument ("seed", ""),

    Option ("centre_radius", "Radius from origin that the fibre centres will be initiated on.")
    + Argument("centre_radius","").type_float(SMALL_FLOAT, -1.0,LARGE_FLOAT),

    IMAGE_PARAMETERS,

    Option()};

std::vector<Triple<double> > distributed_centres(const Triple<double>& img_offsets,
                                                 const Triple<double>& roi_extent,
                                                 size_t num_fibres, gsl_rng* rand_gen,
                                                 double reject_radius);

std::vector<Triple<double> > distant_centres(size_t num_fibres, double centre_radius,
                                             gsl_rng* rand_gen);

std::vector<Triple<double> > orientations(size_t num_fibres, gsl_rng* rand_gen,
                                          double length_stddev);

EXECUTE {
    
        std::string output_location = argument[0];
        
        double acs = NAN;
        double base_intensity = BASE_INTENSITY_DEFAULT;
        size_t degree = DEGREE_DEFAULT;
        size_t num_fibres = NUM_FIBRES_DEFAULT;
        double length_stddev = LENGTH_STDDEV_DEFAULT;
        double width_mean = WIDTH_MEAN_DEFAULT;
        double width_stddev = WIDTH_STDDEV_DEFAULT;
        double reject_radius = REJECT_RADIUS_DEFAULT;
        double edge_buffer = EDGE_BUFFER_DEFAULT;
        double centre_radius = -1.0;
        double curve_stddev = 0;
        double width_epsilon = Fibre::Tractlet::Set::WIDTH_EPSILON_DEFAULT;
        double length_epsilon = Fibre::Tractlet::Set::LENGTH_EPSILON_DEFAULT;
        size_t seed = time(NULL);
        
        Options opt;
        
        opt = get_options("acs");
        if (opt.size())
            acs = opt[0][0];
        
        opt = get_options("base_intensity");
        if (opt.size())
            base_intensity = opt[0][0];
        
        opt = get_options("degree");
        if (opt.size())
            degree = opt[0][0];
        
        opt = get_options("num_fibres");
        if (opt.size())
            num_fibres = opt[0][0];
        
        opt = get_options("length_stddev");
        if (opt.size())
            length_stddev = opt[0][0];
        
        opt = get_options("width_mean");
        if (opt.size())
            width_mean = opt[0][0];
        
        opt = get_options("width_stddev");
        if (opt.size())
            width_stddev = opt[0][0];
        
        opt = get_options("reject_radius");
        if (opt.size())
            reject_radius = opt[0][0];
        
        opt = get_options("edge_buffer");
        if (opt.size())
            edge_buffer = opt[0][0];
        
        opt = get_options("centre_radius");
        if (opt.size())
            centre_radius = opt[0][0];
        
        opt = get_options("curve_stddev");
        if (opt.size())
            curve_stddev = opt[0][0];
        
        opt = get_options("width_epsilon");
        if (opt.size())
            width_epsilon = opt[0][0];
        
        opt = get_options("length_epsilon");
        if (opt.size())
            length_epsilon = opt[0][0];
        
        opt = get_options("seed");
        if (opt.size()) {
            std::string seed_string = opt[0][0];
            seed = to<size_t>(seed_string);
        } else
            std::cout << "No random seed provided, using: " << seed << std::endl;
        
        SET_IMAGE_PARAMETERS;
        
        if (!img_offsets.valid())
            img_offsets = Image::Observed::Buffer::default_corner_offset(img_dims, img_vox_lengths);
        
        Triple<double> roi_extent = (img_vox_lengths * img_dims) / 2.0 - edge_buffer;
        
        if (isnan(reject_radius)) {
            double vol = roi_extent[X] * roi_extent[Y] * roi_extent[Z] * 8.0;
            
            reject_radius = MR::Math::pow((3.0 * vol) / (4.0 * M_PI * (double) num_fibres),
                                    1.0 / 3.0)
                            * DEFAULT_REJECT_RADIUS_SCALE;
        }
        
        std::map<std::string, std::string> properties;
        
        ADD_IMAGE_PROPERTIES(properties);
        properties["length_stddev"] = str(length_stddev);
        properties["width_stddev"] = str(width_stddev);
        
        if (centre_radius > 0)
            properties["centre_radius"] = str(centre_radius);
        else
            properties["reject_radius"] = str(reject_radius);
        
        properties["edge_buffer"] = str(edge_buffer);
        properties["seed"] = str(seed);
        
        //Initialise random number generator.
        gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(rand_gen, seed);
        
        std::vector<const char*> props;
        if (base_intensity)
            props.push_back(Fibre::Strand::Set::BASE_INTENSITY_PROP);
        
        std::vector<const char*> elem_props;
        if (acs >= 0.0)
            elem_props.push_back(Fibre::Base::Object::ALPHA_PROP);
        
        if (File::has_or_txt_extension<Fibre::Strand>(output_location)) {
            
            MR::ProgressBar progress_bar("Initialising strands...");
            
            Fibre::Strand::Set strands(num_fibres, degree, props, elem_props);
            strands.zero();
            strands.set_extend_props(properties);
            if (base_intensity)
                strands.set_base_intensity(base_intensity);
            
            std::vector<Triple<double> > centres;
            // Select between initialising the strands on a shell of "centre_radius" or uniformly distributed with an roi
            if (centre_radius > 0)
                centres = distant_centres(num_fibres, centre_radius, rand_gen);
            else
                centres = distributed_centres(img_offsets, roi_extent, num_fibres, rand_gen,
                        reject_radius);
            
            std::vector<Triple<double> > orients = orientations(num_fibres, rand_gen,
                    length_stddev);
            
            for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i) {
                
                strands[strand_i][0] = centres[strand_i];
                strands[strand_i][1] = orients[strand_i];
                
                if (curve_stddev)
                    strands[strand_i][2] = rand_triple(curve_stddev, rand_gen);
                
                if (!isnan(acs))
                    strands[strand_i].set_acs(acs);
                else
                    strands[strand_i].set_acs(1.0);
                
            }
            
            strands.save(output_location);
            
        } else if (File::has_or_txt_extension<Fibre::Tractlet>(output_location)) {
            
            MR::ProgressBar progress_bar("Initialising tractlets...");
            
            Fibre::Tractlet::Set tractlets(num_fibres, degree, props, elem_props);
            tractlets.zero();
            tractlets.set_width_epsilon(width_epsilon);
            tractlets.set_length_epsilon(length_epsilon);
            tractlets.set_extend_props(properties);
            if (base_intensity)
                tractlets.set_base_intensity(base_intensity);
            
            std::vector<Triple<double> > centres;
            if (centre_radius > 0)
                centres = distant_centres(num_fibres, centre_radius, rand_gen);
            else
                centres = distributed_centres(img_offsets, roi_extent, num_fibres, rand_gen,
                        reject_radius);
            
            std::vector<Triple<double> > orients = orientations(num_fibres, rand_gen,
                    length_stddev);
            
            for (size_t tractlet_i = 0; tractlet_i < tractlets.size(); ++tractlet_i) {
                
                tractlets[tractlet_i][0][0] = centres[tractlet_i];
                tractlets[tractlet_i][0][1] = orients[tractlet_i];
                
                if (curve_stddev)
                    tractlets[tractlet_i][0][2] = rand_triple(curve_stddev, rand_gen);
                
                //Get width axes perpendicular to the main orientation and each other.
                Triple<double> ax1 =
                        orients[tractlet_i].cross(orients[tractlet_i].min_axis()).normalise();
                Triple<double> ax2 = orients[tractlet_i].cross(ax1).normalise();
                
                //Scale width axes by width parameter
                tractlets[tractlet_i][1][0] = ax1
                        * (width_mean + gsl_ran_gaussian(rand_gen, width_stddev));
                tractlets[tractlet_i][2][0] = ax2
                        * (width_mean + gsl_ran_gaussian(rand_gen, width_stddev));
                
                if (!isnan(acs))
                    tractlets[tractlet_i].set_acs(acs);
                else
                    tractlets[tractlet_i].normalise_density();
            }
            
            tractlets.save(output_location);
            
        } else
            throw Exception(
                    "Unrecognised file extension '" + File::extension(output_location) + "'.");
        
    }
    
//Randomly generate centre points withinness the region of interest
    std::vector<Triple<double> > distributed_centres(const Triple<double>& img_offsets,
                                                     const Triple<double>& roi_extent,
                                                     size_t num_fibres, gsl_rng* rand_gen,
                                                     double reject_radius) {
        
        std::vector<Triple<double> > points;
        
        for (size_t fibre_i = 0; fibre_i < num_fibres; ++fibre_i) {
            
            Triple<double> centre;
            
            bool accept = false;
            
            for (size_t attempt_i = 0; attempt_i < MAX_NUM_ATTEMPTS; ++attempt_i) {
                
                for (size_t dim_i = 0; dim_i < 3; ++dim_i)
                    centre[dim_i] = gsl_ran_flat(rand_gen, img_offsets[dim_i] - roi_extent[dim_i],
                            img_offsets[dim_i] + roi_extent[dim_i]);
                
                bool reject = false;
                
                for (size_t point_i = 0; point_i < points.size(); ++point_i)
                    if (dist(centre, points[point_i]) < reject_radius) {
                        reject = true;
                        break;
                    }
                
                if (!reject) {
                    accept = true;
                    break;
                }
                
            }
            
            if (accept)
                points.push_back(centre);
            else
                throw Exception(
                        "Could not find centre position for fibre " + str(fibre_i)
                        + " that lies outside the rejection radius of previously generated fibre centres after "
                        + str(MAX_NUM_ATTEMPTS) + ". Suggest reducing the reject_radius ("
                        + str(reject_radius) + ").");
            
        }
        
        return points;
        
    }
    
    std::vector<Triple<double> > distant_centres(size_t num_fibres, double centre_radius,
                                                 gsl_rng* rand_gen) {
        
        std::vector<Triple<double> > points;
        
        for (size_t fibre_i = 0; fibre_i < num_fibres; ++fibre_i) {
            
            Triple<double> point = rand_triple(1.0, rand_gen);
            
            point.normalise();
            
            point *= centre_radius;
            
            points.push_back(point);
            
        }
        
        return points;
        
    }
    
//Randomly generate orientations
    std::vector<Triple<double> > orientations(size_t num_fibres, gsl_rng* rand_gen,
                                              double length_stddev) {
        
        std::vector<Triple<double> > orients;
        
        for (size_t fibre_i = 0; fibre_i < num_fibres; ++fibre_i) {
            
            Triple<double> orient;
            
            for (size_t dim_i = 0; dim_i < 3; ++dim_i)
                orient[dim_i] = gsl_ran_gaussian(rand_gen, length_stddev);
            
            orients.push_back(orient);
            
        }
        
        return orients;
        
    }
