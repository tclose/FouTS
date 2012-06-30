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
#include "bts/fibre/properties/extended.h"

#include "bts/fibre/base/writer.h"



#include "bts/image/expected/buffer.h"


#include "bts/inline_functions.h"


using namespace BTS;


const double MIN_ACS_DEFAULT = 0.001;

template <typename T> void average_fibres(const std::string& input_location, const std::string& output_location,
                                                                const std::vector<size_t>& include, double min_acs);

template <typename T> double acs(const T& fibre);

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Average a set of fibres sets, inluding the option to remove thin tractlets.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("input", "The fibres to be selected from.").type_file (),
  Argument ("output", "The selected fibres.").type_file (),
  Argument()
};


OPTIONS = {

  Option ("include","The indices of the strands to include")
    + Argument ("include",""),

   Option ("min_acs","Minimum area of tractlet to include")
     + Argument ("include","").type_float (0,LARGE_FLOAT, MIN_ACS_DEFAULT),
   

Option() };


EXECUTE {

  std::string input_location = argument[0];
  std::string output_location = argument[1];

  std::vector<size_t> include;
  double min_acs;

  Options opt = get_options("include");
  if (opt.size())
    include = parse_sequence<size_t> (opt[0][0]);


  if (File::has_extension<Fibre::Strand::Set>(input_location)) {

    if (get_options("min_acs").size())
      throw Exception("Cannot use min_acs option when averaging strand sets (because they have none).");

    min_acs = -1;

  } else if (File::has_extension<Fibre::Tractlet::Set>(input_location)) {

    opt = get_options("min_acs");
    if (opt.size())
      min_acs = opt[0][0];
    else
      min_acs = MIN_ACS_DEFAULT;

  }  else
    throw Exception ("Unsupported extension of file '" + input_location + "'.");



  if (File::has_extension<Fibre::Strand::Set>(input_location)) {

    average_fibres<Fibre::Strand>(input_location, output_location, include, min_acs);

  } else if (File::has_extension<Fibre::Tractlet::Set>(input_location)) {

    average_fibres<Fibre::Tractlet>(input_location, output_location, include, min_acs);

  } else
    throw Exception ("Unsupported extension of file '" + input_location + "'.");


 }




template <typename T> void average_fibres(const std::string& input_location, const std::string& output_location,
                                                                const std::vector<size_t>& include, double min_acs) {


  if (!File::has_extension<T>(output_location))
    throw Exception ("Mismatching extension of output location '" + output_location + "' for input '" + input_location + "'.");


  typename T::Set fibres, avg_fibres;

  typename T::Set::Reader reader(input_location);

  size_t set_count = to<size_t>(reader.get_extend_props()["count"]);

  //Read first tractlets set from input file, save to avg_fibres to set correct size then rewind the reader back to the start.
  if (!reader.next(avg_fibres))
    throw Exception ("No sets loaded from input location '" + input_location + "'.");

  avg_fibres.zero();
  reader.rewind();

  std::vector<double> avg_acss(avg_fibres.size());

  //Variables to loop through tract sets with.
  size_t set_i = 0;
  std::vector<size_t> includes_not_found = include;

  MR::ProgressBar avg_progress_bar("Generating average set...", set_count);

  while (reader.next(fibres)) {

    //Check to see if current set is listed in included list.
    std::vector<size_t>::iterator include_it = find(includes_not_found.begin(), includes_not_found.end(), set_i);

    bool include_set = true;

    if (include.size()) {
      if (include_it == includes_not_found.end())
        include_set = false;
      else
        //Remove matched include from list.
        includes_not_found.erase(include_it);
    }

    if (include_set) {

      if (fibres.size() != avg_fibres.size())
        throw Exception ("Size of set " + str(set_i) + " (" + str(fibres.size()) + ") does not match that of initial set (" + str(avg_fibres.size()) + ").");

      avg_fibres += fibres;

      for (size_t elem_i = 0; elem_i < fibres.size(); ++elem_i)
        avg_acss[elem_i] += acs(fibres[elem_i]);

    }

    ++set_i;
    ++avg_progress_bar;

  }

  if (includes_not_found.size() != 0)
    throw Exception ("Indices, " + str(includes_not_found) + ", were outside the range of the number of loaded sets.");

  size_t num_included_sets;

  if (include.size())
    num_included_sets = include.size();
  else
    num_included_sets = set_i;


  avg_fibres /= (double)num_included_sets;

  typename T::Set output_set(std::vector<const char*>(), reader.elem_prop_keys());

  MR::ProgressBar min_vol_progress_bar ("Excluding fibres that are thinner than average acs...", avg_fibres.size());

  for (size_t elem_i = 0; elem_i < avg_fibres.size(); ++elem_i) {

    avg_acss[elem_i] /= (double)num_included_sets;

    if (avg_acss[elem_i] >= min_acs) {

      output_set.push_back(avg_fibres[elem_i]);

    }


    ++min_vol_progress_bar;
  }


  output_set.save(output_location);

}


template <> double acs<Fibre::Strand>(const Fibre::Strand& fibre) {

  return 0.0;

}


template <> double acs<Fibre::Tractlet>(const Fibre::Tractlet& fibre) {

  return fibre.acs();

}


