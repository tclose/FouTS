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

#include "bts/cmd.h"
#include "point.h"
#include "progressbar.h"

#include "bts/common.h"
#include "bts/file.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/fibre/base/writer.h"

#include "bts/image/expected/buffer.h"

#include "bts/inline_functions.h"

const static double STEP_SIZE_DEFAULT = 0.05;
const static double PER_ACS_DEFAULT = 10000;
const static size_t DEGREE_DEFAULT = 50;
const static size_t NUM_POINTS_DEFAULT = 100;

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Select a subset of strands",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The strands to be selected from.").type_file (),
    Argument ("output", "The selected strands.").type_file (),
    Argument()
};

OPTIONS= {

    Option ("step_size", "The step size of the generated tracks")
    + Argument ("step_size", "").type_float (0, STEP_SIZE_DEFAULT, LARGE_FLOAT),

    Option ("per_acs","Number of tracks per ACS of tract.")
    + Argument ("per_acs","").type_float (0.0, PER_ACS_DEFAULT, LARGE_FLOAT),

    Option ("degree","The degree of the Fourier tract tracks will be converted"
                     "into before generating the tdi tracks.")
    + Argument ("degree","").type_integer(0, DEGREE_DEFAULT, LARGE_INT),

    Option ("num_points","The number of points to plot before they are resampled.")
    + Argument ("num_points","").type_integer(0, NUM_POINTS_DEFAULT, LARGE_INT),

    Option()};


EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location = argument[1];
        
        size_t degree = DEGREE_DEFAULT;
        double step_size = STEP_SIZE_DEFAULT;
        double per_acs = PER_ACS_DEFAULT;
        size_t num_points = NUM_POINTS_DEFAULT;
        
        Options opt = get_options("degree");
        if (opt.size())
            degree = opt[0][0];

        opt = get_options("step_size");
        if (opt.size())
            step_size = opt[0][0];
        
        opt = get_options("per_acs");
        if (opt.size())
            per_acs = opt[0][0];

        opt = get_options("num_points");
        if (opt.size())
            num_points = opt[0][0];

        if (!File::has_extension<Fibre::Track>(output_location))
            throw Exception("Only supports writing to 'tck' format (output "
                            "provided " + input_location + ").");

        // Get input tractlets
        Fibre::Tractlet::Set input_tracts;
        std::map<std::string, std::string> properties_row;

        if (File::has_extension<Fibre::Tractlet>(input_location)) {

            Fibre::Tractlet::Reader reader(input_location);
            Fibre::Tractlet tract;
            while (reader.next(tract, properties_row))
                input_tracts.push_back(tract);

        } else if (File::has_extension<Fibre::Track>(input_location)) {
            Fibre::Track::Reader reader(input_location);
            Fibre::Track tck;
            while (reader.next(tck, properties_row)) {
                double width = to<double>(properties_row[Fibre::Track::RADIUS_PROP]);
                Fibre::Tractlet tract(tck, degree, width);
                tract.normalise_density();
                input_tracts.push_back(tract);
            }

        } else
            throw Exception("Invalid input extension '" + File::extension(input_location) +
                            "' (only accepts '.tct', or '.tck')");

        Fibre::Properties properties;
        properties["step_size"] = str(step_size);
        properties["per_acs"] = str(per_acs);
        properties["degree"] = str(degree);
        properties["source"] = input_location;

        Fibre::Track::Writer writer(output_location, Fibre::Track::Set(), properties);
        MR::ProgressBar progress_bar("Creating TDI-ready tracks from " +
                                     str(input_tracts.size()) +
                                     " Fourier tracts.", input_tracts.size());

        for (size_t tract_i = 0; tract_i < input_tracts.size(); ++tract_i) {
            const Fibre::Tractlet& tract = input_tracts[tract_i];
            double num_strands = per_acs * tract.acs();
            size_t num_width_sections = (size_t)MR::Math::round(
                                 MR::Math::sqrt(num_strands * 4.0 / M_PI));
            double width_fraction = 1.0 / (double)num_width_sections;
            // Loop across cross-section of tractlet. Loop each perturbation axis from a 'width_fraction' from -1 to a width_fraction
            // from +1, in intervals of 2*width_fraction's. It is performed this way because each strand is to represent the
            // tractlet +- a width_fraction radius about it.
            for (double ax1_frac = (-1.0 + width_fraction); ax1_frac < 1.0; ax1_frac += 2.0 * width_fraction)
                for (double ax2_frac = (-1.0 + width_fraction); ax2_frac < 1.0; ax2_frac += 2.0 * width_fraction)
                    //Include section if coord consisting of the combined fractions along each perturbation axes lies withinness the
                    // unit circle minus a half a width_fraction. This will give the resulting cross-sections an elliptical
                    // shape (or circular if axes are equal).
                    if (MR::Math::pow2(ax1_frac) + MR::Math::pow2(ax2_frac) <= (1.0 - width_fraction / 2.0)) {
                        Fibre::Strand strand(tract.degree());
                        strand.base() = tract[0] + tract[1] * ax1_frac + tract[2] * ax2_frac;
                        Fibre::Track tck(strand, num_points);
                        Fibre::Track resampled_tck = tck.resample(step_size);
                        writer.append(resampled_tck);
                    }
        }
}

