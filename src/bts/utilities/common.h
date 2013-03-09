/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Aug 12, 2010.

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

#ifndef bts_utilities_common_h__
#define bts_utilities_common_h__

#include <map>
#include <string>
#include <vector>

namespace BTS {
    
    namespace Utilities {
        
        typedef std::map<std::string, std::string> Properties;
        
        typedef std::vector<std::string> Header;
    
    }

}

#endif /* bts_utilities_common_h__ */
