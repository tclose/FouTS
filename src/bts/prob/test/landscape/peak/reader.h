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


#ifndef __bts_prob_testpeaks_peak_reader__
#define __bts_prob_testpeaks_peak_reader__
    
namespace BTS {

  namespace Prob {
       
    class Test::Landscape::Peak::Reader {
    
      protected:
        
        std::ifstream fin;
        
      public:
      
        Reader() {}
        
        Reader(const std::string& location)
          { this->open(location); }
        
        ~Reader() { if (fin.is_open()) this->close(); }
        
        void  open(const std::string& location) 
          { fin.open(location.c_str()); }
    
        void  close()
          { fin.close(); }
          
        bool next(Peak& peak) { 
            
          std::string line;  
          bool is_next;
            
          if (std::getline(fin, line)) {  
            
            std::istringstream line_stream(line);

            line_stream >> peak.type >> peak.ndims >> peak.height;

            peak.centre.resize(peak.ndims);
            peak.widths.resize(peak.ndims, peak.ndims);

            for (size_t dim_i = 0; dim_i < peak.ndims; ++dim_i)
              line_stream >> peak.centre[dim_i];

            for (size_t dim_i1 = 0; dim_i1 < peak.ndims; ++dim_i1)
              for (size_t dim_i2 = 0; dim_i2 < peak.ndims; ++dim_i2)
                line_stream >> peak.widths(dim_i1, dim_i2);

//            line_stream >> peak.centre[X] >> peak.centre[Y] >> peak.width[X]  >> peak.width[Y] >> peak.height >> peak.type;
            
            is_next = true;
            
          } else
            is_next = false;
            
                      
          return is_next;
        }
              
    
    };
    
  }
  
}

#endif    
