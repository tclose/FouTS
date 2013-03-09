/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 24/11/2010.

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

#ifndef __bts_mcmc_naninfexception_h__
#define __bts_mcmc_naninfexception_h__

namespace BTS {
    
    namespace MCMC {
        
        class NanInfException: public MR::Exception {
                
            public:
                
                NanInfException()
                        : MR::Exception("Nan or Inf values found.") {
                }
                
                NanInfException(const NanInfException& e)
                        : MR::Exception(e) {
                }
                
                NanInfException& operator=(const NanInfException& e) {
                    MR::Exception::operator=(e);
                    return *this;
                }
                
                ~NanInfException() {
                }
                
        };
    
    }

}

#endif /* __bts_mcmc_naninfexception_h__ */
