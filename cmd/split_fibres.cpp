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

#include "bts/diffusion/model.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/quartic/buffer.h"

#include "bts/prob/uniform.h"
#include "bts/prob/prior.h"

#include "bts/math/common.h"

#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood/one_sided_gaussian.h"
#include "bts/prob/likelihood/gaussian.h"

#include "bts/image/observed/buffer.h"

#include "bts/image/expected/buffer.h"

#include "bts/math/blossom.h"
#include "bts/math/blossom/edge.h"

#include "bts/inline_functions.h"

using namespace BTS;

const double THRESHOLD_DEFAULT = -0.5;
const size_t NUM_LENGTH_SECTIONS_DEFAULT = 5;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
  "Split Fourier tracts into two to allow them to follow branching structure",
  "",
  NULL
};

ARGUMENTS= {
  Argument ("input", "The strands to be selected from.").type_file (),
  Argument ("output", "The selected strands.").type_file (),
  Argument()
};

OPTIONS= {

  Option ("orientation", "The orientation of the plane to split the Fourier "
      "tract into.",
      "(1,0,0).")
  + Argument ("offsets", "").type_text ("auto"),

  Option()};

EXECUTE {

    std::string obs_image_location = argument[0];
    std::string input_location = argument[1];
    std::string output_location = argument[2];

    Triple<double> orient(1.0, 0.0, 0.0);

    Options opt = get_options("orientation");
    if (opt.size())
      orient = parse_triple<size_t>(std::string(opt[0][0]));

    orient.normalise();

    Fibre::Tractlet::Set tracts(input_location);
    Fibre::Tractlet::Set output_tracts;

    for (size_t tract_i = 0; tract_i < tracts.size(); ++tract_i) {

      Fibre::Tractlet t = tracts[tract_i];
      Fibre::Tractlet t1, t2;

      double disp = MR::Math::abs(orient.dot(t[1][0])) +
                    MR::Math::abs(orient.dot(t[2][0])));

      t[1][0] *= scale;
      t[2][0] *= scale

      t1 = t2 = t;

      t1[0][0] -= disp * orient;
      t2[0][0] -= disp * orient;

      output_tracts.push_back(t1);
      output_tracts.push_back(t2);

    }

    output_tracts.save(output_location);

  }

