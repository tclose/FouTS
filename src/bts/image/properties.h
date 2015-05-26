/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 21/07/2010.

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

#ifndef __bts_image_properties_h__
#define __bts_image_properties_h__

#include <map>

namespace FTS {
    
    namespace Image {
        
        class Properties: public std::map<std::string, std::string> {
                
                //Public static variables, nested classes and typedefs
            public:
                
                //Protected member variables
            protected:
                
                //Public member functions
            public:
                
                //Protected member functions
            protected:
                
        };
    
    }

}

#endif /* __bts_image_properties_h__ */
