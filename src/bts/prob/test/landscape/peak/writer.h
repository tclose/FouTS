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


#ifndef __bts_prob_testpeaks_peak_writer__
#define __bts_prob_testpeaks_peak_writer__
    
namespace BTS {

  namespace Prob {
       

    class Test::Landscape::Peak::Writer {
    
        
      protected:
        
        std::ofstream fout;
        
      public:
      
        Writer() {}
        
        Writer(const std::string& location)
          { this->create(location); }
        
        ~Writer() { if (fout.is_open()) this->close(); }
        
        void  create(const std::string& location) 
          { fout.open(location.c_str()); }
    
        void  close()
          { fout.close(); }
          
        void append(const Peak& peak) {

          fout << peak.type << " " << peak.ndims << " " << peak.height;

          for (size_t dim_i = 0; dim_i < peak.ndims; ++dim_i)
            fout << " " << peak.centre[dim_i];

          for (size_t dim_i1 = 0; dim_i1 < peak.ndims; ++dim_i1)
            for (size_t dim_i2 = 0; dim_i2 < peak.ndims; ++dim_i2)
              fout << " " << peak.widths(dim_i1,dim_i2);

          fout << std::endl;

        }
    
    
    };
    
  }
  
}

#endif    
