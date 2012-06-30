/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Mar 14, 2011.

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


#include "bts/cmd.h"

#include "bts/common.h"

#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"

#include "bts/inline_functions.h"

using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);


DESCRIPTION = {
  "Map colours to reference set by switching bundle indices to closest match.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("reference","Reference fibres the match fibres will be mapped to.").type_file (),
  Argument ("samples","Fibre samples that stats will collected from.").type_file (),
  Argument ("output","Output filename, if not provided will be same as samples with 'stat' extension.").optional().type_file (),
  Argument()
};


OPTIONS = {

Option() };


EXECUTE {

  std::string reference_location = argument[0];
  std::string samples_location = argument[1];
	std::string output_location = samples_location;

	if (argument.size() == 3)
	  output_location = argument[2].c_str();


	if (File::has_extension<Fibre::Strand>(reference_location)) {

    if (!File::has_extension<Fibre::Strand::Set>(samples_location))
      throw Exception ("Samples location ('" + samples_location + "') does not match reference location ('" + reference_location + "').");

    Fibre::Strand::Set reference(reference_location);

    Fibre::Strand::Set strands;

    Fibre::Strand::Set::Reader reader(samples_location);

    size_t count = 0;

    Fibre::Strand::Set average(reference);

    average.zero();

    MR::ProgressBar progress_bar("Calculating average...", to<size_t>(reader.get_extend_props()["count"]));

    std::ofstream all_out((File::strip_extension(output_location) + ".dist.txt").c_str());
    std::ofstream strand_out((File::strip_extension(output_location) + ".strdist.txt").c_str());
    std::ofstream pos_out((File::strip_extension(output_location) + ".pos.txt").c_str());
    std::ofstream orient_out((File::strip_extension(output_location) + ".orient.txt").c_str());
    std::ofstream len_out((File::strip_extension(output_location) + ".len.txt").c_str());
    std::ofstream curv_out((File::strip_extension(output_location) + ".curv.txt").c_str());


    while (reader.next(strands)) {

      if (strands.size() != reference.size())
        throw Exception ("Size of tractlet set " + str(count) + " (" + str(strands.size()) + ") does not match that of 'reference' set (" + str(reference.size()) + ")");

      strands = strands.smallest_distance_set(reference);

      average += strands;

      all_out << strands.distance(reference) << std::endl;

      for (size_t strand_i = 0; strand_i < strands.size(); ++ strand_i) {

        strand_out << strands[strand_i].distance(reference[strand_i]) << " ";
        pos_out << (strands[strand_i][0] - reference[strand_i][0]).norm() << " ";
        orient_out << strands[strand_i][1].angle(reference[strand_i][1]) << " ";
        len_out << strands[strand_i][1].norm() - reference[strand_i][1].norm() << " ";
        curv_out << (strands[strand_i][2] - reference[strand_i][2]).norm() * SQRT_2 << " ";



      }

      strand_out << std::endl;
      pos_out << std::endl;
      orient_out << std::endl;
      len_out << std::endl;
      curv_out << std::endl;


      ++count;
      ++progress_bar;

    }


    average /= (double)count;

    average.save(File::strip_extension(output_location) + ".avg.str");

    Fibre::Strand::Set variance(reference);

    variance.zero();

    count = 0;

    reader.rewind();

    while (reader.next(strands)) {

      if (strands.size() != reference.size())
        throw Exception ("Size of tractlet set " + str(count) + " (" + str(strands.size()) + ") does not match that of 'reference' set (" + str(reference.size()) + ")");

      Fibre::Strand::Set deviation = strands.smallest_distance_set(reference);

      deviation -= average;

      deviation *= deviation;

      variance += deviation;

      ++count;
      ++progress_bar;

    }

    variance /= (double)count;

    variance.save(File::strip_extension(output_location) + ".var.str");

	} else if (File::has_extension<Fibre::Tractlet>(reference_location)) {

    if (!File::has_extension<Fibre::Tractlet::Set>(samples_location))
      throw Exception ("Samples location ('" + samples_location + "') does not match reference location ('" + reference_location + "').");

    Fibre::Tractlet::Set reference(reference_location);

    Fibre::Tractlet::Set tractlets;

    Fibre::Tractlet::Set::Reader reader(samples_location);

    size_t count = 0;

    Fibre::Tractlet::Set average(reference);

    average.zero();

    MR::ProgressBar progress_bar("Calculating average...", to<size_t>(reader.get_extend_props()["count"]));

    std::ofstream all_out((File::strip_extension(output_location) + ".dist.txt").c_str());
    std::ofstream tractlet_out((File::strip_extension(output_location) + ".tctdist.txt").c_str());
    std::ofstream pos_out((File::strip_extension(output_location) + ".pos.txt").c_str());
    std::ofstream perp_out((File::strip_extension(output_location) + ".perp.txt").c_str());
    std::ofstream orient_out((File::strip_extension(output_location) + ".orient.txt").c_str());
    std::ofstream len_out((File::strip_extension(output_location) + ".len.txt").c_str());
    std::ofstream curv_out((File::strip_extension(output_location) + ".curv.txt").c_str());

    std::ofstream width1_out((File::strip_extension(output_location) + ".width1.txt").c_str());
    std::ofstream width2_out((File::strip_extension(output_location) + ".width2.txt").c_str());

    std::ofstream acs_out((File::strip_extension(output_location) + ".area.txt").c_str());

    std::ofstream rot_out((File::strip_extension(output_location) + ".rot.txt").c_str());

    while (reader.next(tractlets)) {

      if (tractlets.size() != reference.size())
        throw Exception ("Size of tractlet set " + str(count) + " (" + str(tractlets.size()) + ") does not match that of 'reference' set (" + str(reference.size()) + ")");

      tractlets = tractlets.smallest_distance_set(reference);

      average += tractlets;

      all_out << tractlets.distance(reference) << std::endl;

      for (size_t tractlet_i = 0; tractlet_i < tractlets.size(); ++ tractlet_i) {

        tractlet_out << tractlets[tractlet_i].distance(reference[tractlet_i]) << " ";

        pos_out << (tractlets[tractlet_i][0][0] - reference[tractlet_i][0][0]).norm() << " ";
        perp_out << (tractlets[tractlet_i][0][0] - tractlets[tractlet_i][0][0].dot(tractlets[tractlet_i][0][1])/tractlets[tractlet_i][0][1].norm()
                    - reference[tractlet_i][0][0] + reference[tractlet_i][0][0].dot(reference[tractlet_i][0][1])/reference[tractlet_i][0][1].norm()).norm() << " ";
        orient_out << tractlets[tractlet_i][0][1].angle(reference[tractlet_i][0][1]) << " ";
        len_out << tractlets[tractlet_i][0][1].norm() - reference[tractlet_i][0][1].norm() << " ";
        curv_out << (tractlets[tractlet_i][0][2] - reference[tractlet_i][0][2]).norm() * SQRT_2 << " ";

        width1_out << (tractlets[tractlet_i][1][0] - reference[tractlet_i][1][0]).norm() * SQRT_2 << " ";
        width2_out << (tractlets[tractlet_i][2][0] - reference[tractlet_i][2][0]).norm() * SQRT_2 << " ";

        acs_out << tractlets[tractlet_i].acs() - reference[tractlet_i].acs() << " ";

        rot_out << tractlets[tractlet_i].rotation() << " ";

      }

      tractlet_out << std::endl;
      pos_out << std::endl;
      perp_out << std::endl;
      orient_out << std::endl;
      len_out << std::endl;
      curv_out << std::endl;

      width1_out << std::endl;
      width2_out << std::endl;
      acs_out << std::endl;

      rot_out << std::endl;


      ++count;
      ++progress_bar;

    }


    average /= (double)count;

    average.save(File::strip_extension(output_location) + ".avg.tct");

    Fibre::Tractlet::Set variance(reference);

    variance.zero();

    count = 0;

    reader.rewind();

    while (reader.next(tractlets)) {

      if (tractlets.size() != reference.size())
        throw Exception ("Size of tractlet set " + str(count) + " (" + str(tractlets.size()) + ") does not match that of 'reference' set (" + str(reference.size()) + ")");

      Fibre::Tractlet::Set deviation = tractlets.smallest_distance_set(reference);

      deviation -= average;

      deviation *= deviation;

      variance += deviation;

      ++count;
      ++progress_bar;

    }

    variance /= (double)count;

    variance.save(File::strip_extension(output_location) + ".var.tct");

  } else
    throw Exception ("Usupported output file type ('" + output_location + "') for match of type tractlet.");


}
