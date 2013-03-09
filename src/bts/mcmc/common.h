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

#ifndef __bts_mcmc_common_h__
#define __bts_mcmc_common_h__

namespace BTS {
    
    namespace MCMC {
        
        const std::string LOG_PROB_PROP = "log_px";
        const std::string ANNEAL_LOG_PROB_PROP = "anneal_log_px";
        const std::string ACCEPTANCE_RATIO_PROP = "acceptance_ratio";
        const std::string ELAPSED_TIME_PROP = "elapsed_time";
        const std::string H_PROP = "H";
        const std::string PROPOSED_H_PROP = "Proposed H";
        
        typedef std::map<std::string, std::string> SampleProperties;
        typedef std::vector<std::string> SamplePropertyHeader;
    
    }

}

#endif

