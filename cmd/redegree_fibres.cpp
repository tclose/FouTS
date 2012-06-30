/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Jul 26, 2010.

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

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/file.h"

#include "bts/math/svd.h"



#include "bts/inline_functions.h"

using namespace BTS;

const size_t DEGREE_DEFAULT = 3;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);


DESCRIPTION = {
  "redegree_fibres",
  "Changes the degree of the fibres",
  NULL
};

ARGUMENTS = {
  Argument ("input location","").type_file (),
  Argument ("output location","").optional().type_file (),
  Argument()
};


OPTIONS = {

  Option("degree", "The degree of the output fibres..")
   + Argument("degree","").type_integer(0,DEGREE_DEFAULT,LARGE_INT),

Option() };


EXECUTE {


  std::string input_location = argument[0];
  std::string output_location = input_location;

  if (argument.size() == 2)
    output_location = argument[1].c_str();

  size_t degree = DEGREE_DEFAULT;

  Options opt = get_options("degree");
  if (opt.size())
    degree = opt[0][0];

  MR::ProgressBar progress_bar ("Redegreeing fibres to a degree of " + str(degree) + " ...");

  if (File::has_extension<Fibre::Strand>(input_location)) {

    Fibre::Strand::Set strands (input_location);

    for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i)
      strands[strand_i].resize(degree,0.0);

    strands.save(output_location);

  } else if (File::has_extension<Fibre::Tractlet>(input_location)) {

    Fibre::Tractlet::Set tractlets (input_location);

    for (size_t tract_i = 0; tract_i < tractlets.size(); ++tract_i)
      tractlets[tract_i].redegree(degree,0.0);

    tractlets.save(output_location);

  } else
	  throw Exception ("Unrecognised extension '" + File::extension(input_location) + "'.");


}














































