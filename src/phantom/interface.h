/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.
 
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

#ifndef __phantom_interface_h__
#define __phantom_interface_h__

#include "bts/fibre/strand/set.h"

//#include "point.h"

//extern "C" {
#include "phantom/shared/strand_collection.h"
//}

std::vector<BTS::Triple<double> >& generate_pre_points(
        BTS::Fibre::Track::Set& strands, std::vector<BTS::Triple<double> >& pre_points);

std::vector<BTS::Triple<double> >& generate_post_points(
        BTS::Fibre::Track::Set& strands, std::vector<BTS::Triple<double> >& post_points);

BTS::Fibre::Track::Set& convert_nfg_to_mr(BTS::Fibre::Track::Set& strands,
                                          std::vector<BTS::Triple<double> >& pre_points,
                                          std::vector<BTS::Triple<double> >& post_points,
                                          Strand_collection* c);

Strand_collection* convert_mr_to_nfg(Strand_collection* c, BTS::Fibre::Track::Set& strands,
                                     std::vector<BTS::Triple<double> >& pre_points,
                                     std::vector<BTS::Triple<double> >& post_points);

#endif
