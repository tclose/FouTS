/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 16/07/2010.

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

#ifndef __bts_mcmc_inline_functions_h__
#define __bts_mcmc_inline_functions_h__


#ifdef __bts_mcmc_proposal_momentum_h__

#include "bts/mcmc/proposal/momentum.cpp.h"

#endif

#ifdef __bts_mcmc_proposal_momentum_weighted_h__

#include "bts/mcmc/proposal/momentum/weighted.cpp.h"

#endif


#ifdef __bts_mcmc_proposal_momentum_weighted_nonseparable_h__

#include "bts/mcmc/proposal/momentum/weighted/non_separable.cpp.h"

#endif

#endif /* __bts_proposal_inline_functions_h__ */
