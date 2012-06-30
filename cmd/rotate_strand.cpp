/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Jan 25, 2011.

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

#include "bts/fibre/strand.h"
#include "bts/fibre/strand/set.h"

#include "bts/inline_functions.h"

using namespace BTS;


SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);


const size_t NUM_STEPS_DEFAULT = 100;
const float ANGLE_DEFAULT = M_PI / 2.0;
const std::string AXIS_DEFAULT = "x";

DESCRIPTION = {
  "rotate_strand",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("input_strand","").type_file (),
  Argument ("output strand set","").type_file (),
  Argument()
};


OPTIONS = {

	Option("num_steps","")
	 + Argument("num_steps").type_integer(0,LARGE_INT,NUM_STEPS_DEFAULT),

  Option("angle","")
   + Argument("angle").type_float(0, 2.0 * M_PI, ANGLE_DEFAULT),

  Option("axis","")
   + Argument("axis").type_text(AXIS_DEFAULT.c_str()),

Option() };


EXECUTE {

  std::string input_location = argument[0];
	std::string output_location = argument[1];

	size_t  num_steps = NUM_STEPS_DEFAULT;
	float angle = ANGLE_DEFAULT;
	std::string axis = AXIS_DEFAULT;

  Options opt;

  opt = get_options("num_steps");
  if (opt.size())
		num_steps = opt[0][0];

  opt = get_options("angle");
  if (opt.size())
    angle = opt[0][0];

  opt = get_options("axis");
  if (opt.size())
    axis = opt[0][0].c_str();

  Fibre::Strand::Set strands (input_location);
  Fibre::Strand::Set rotated_strands(strands);

  std::map<std::string,std::string> file_props;
  file_props["method"] = "rotated_strand";

  Fibre::Strand::Set::Writer writer(output_location, strands, file_props);

  MR::Math::Matrix<double> rotation_matrix(3,3);

  MR::ProgressBar progress_bar ("Generating rotated strands...", num_steps);

  for (size_t step_i = 0; step_i < num_steps; ++step_i) {

    double current_angle = ((double)step_i / (double)(num_steps-1)) * angle;
    double cos_angle = MR::Math::cos(current_angle);
    double sin_angle = MR::Math::sin(current_angle);

    rotation_matrix = 0.0;

    if ((axis == "x") || (axis == "X")) {

      rotation_matrix(X,X) = 1.0;
      rotation_matrix(Y,Y) = rotation_matrix(Z,Z) = cos_angle;
      rotation_matrix(Y,Z) = -sin_angle;
      rotation_matrix(Z,Y) = sin_angle;

    } else if ((axis == "y") || (axis == "Y")) {

      rotation_matrix(Y,Y) = 1.0;
      rotation_matrix(X,X) = rotation_matrix(Z,Z) = cos_angle;
      rotation_matrix(X,Z) = sin_angle;
      rotation_matrix(Z,X) = -sin_angle;

    } else if ((axis == "z") || (axis == "Z")) {

      rotation_matrix(Z,Z) = 1.0;
      rotation_matrix(Y,Y) = rotation_matrix(X,X) = cos_angle;
      rotation_matrix(Y,X) = -sin_angle;
      rotation_matrix(X,Y) = sin_angle;

    } else
      throw Exception ("Value for axis (" + axis + ") is not valid (x, y or z)");

    for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i) {

      MR::Math::Vector<double> vector(3), rotated_vector(3);
      Fibre::Strand rotated_strand(3);

      vector[X] = strands[strand_i][1][X];
      vector[Y] = strands[strand_i][1][Y];
      vector[Z] = strands[strand_i][1][Z];

      MR::Math::mult(rotated_vector, rotation_matrix, vector);

      rotated_strands[strand_i][1][X] = rotated_vector[X];
      rotated_strands[strand_i][1][Y] = rotated_vector[Y];
      rotated_strands[strand_i][1][Z] = rotated_vector[Z];


    }


    writer.append(rotated_strands);
    ++progress_bar;

  }



}
