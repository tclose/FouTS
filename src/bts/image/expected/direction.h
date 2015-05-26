/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Voxelhomas G Close, 5/05/09.

 Voxelhis file is part of FTS.

 FouTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FouTS is distributed in the hope that it will be useful,
 but WIVoxelHOUVoxel ANY WARRANVoxelY; without even the implied warranty of
 MERCHANVoxelABILIVoxelY or FIVoxelNESS FOR A PARVoxelICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_image_expected_direction_h__
#define __bts_image_expected_direction_h__

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            class Direction;
        
        }
    
    }

}

#include "bts/common.h"

#include "bts/image/expected/voxel.h"
#include "bts/diffusion/response.h"

#include "bts/triple.h"
#include "bts/fibre/strand/section.h"

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            class Direction {
                    
                    //Public members
                public:
                    
                    Voxel* voxel;
                    Diffusion::Response* response;

                    //Public methods
                public:
                    
                    Direction(Voxel& voxel, Diffusion::Response& response)
                            : voxel(&voxel), response(&response) {
                    }
                    
                    ~Direction() {
                    }
                    
                    //'signal' functions have explicit Fibre::Strand::Section and Fibre::Tractlet::Section arguments in favour of using
                    // the template directly so that their function signature can be passed as a pointer to GradientTester.
                    double signal(const Fibre::Strand::Section& section) {
                        return signal_tpl(section);
                    }
                    
                    double signal(const Fibre::Tractlet::Section& section) {
                        return signal_tpl(section);
                    }
                    
                    double signal(const Fibre::Strand::Section& section,
                                  Fibre::Strand::BasicSection& gradient) {
                        return signal_tpl(section, gradient);
                    }
                    
                    double signal(const Fibre::Tractlet::Section& section,
                                  Fibre::Strand::BasicSection& gradient) {
                        return signal_tpl(section, gradient);
                    }
                    
                    double signal(const Fibre::Strand::Section& section,
                                  Fibre::Strand::BasicSection& gradient,
                                  Fibre::Strand::BasicSection::Tensor& hessian) {
                        return signal_tpl(section, gradient, hessian);
                    }
                    
                    double signal(const Fibre::Tractlet::Section& section,
                                  Fibre::Strand::BasicSection& gradient,
                                  Fibre::Strand::BasicSection::Tensor& hessian) {
                        return signal_tpl(section, gradient, hessian);
                    }
                    
                    double b_value() const {
                        return response->b_value();
                    }
                    
                    const Diffusion::Encoding& encoding() const {
                        return *response;
                    }
                    
                protected:
                    
                    template<typename T> double signal_tpl(const T& section);

                    template<typename T> double signal_tpl(const T& section,
                                                           Fibre::Strand::BasicSection& gradient);

                    template<typename T> double signal_tpl(const T& section,
                                                           Fibre::Strand::BasicSection& gradient,
                                                           typename T::Tensor& hessian);

                    static void signal_gradient(
                            const Fibre::Strand::Section& section,
                            Fibre::Strand::BasicSection& gradient, double diffusion_weighting,
                            double interpolation_weighting,
                            const Triple<double>& diffusion_gradient,
                            const Fibre::Strand::BasicSection& interpolation_gradient);

                    static void signal_gradient(
                            const Fibre::Tractlet::Section& section,
                            Fibre::Strand::BasicSection& gradient, double diffusion_weighting,
                            double interpolation_weighting,
                            const Triple<double>& diffusion_gradient,
                            const Fibre::Strand::BasicSection& interpolation_gradient);

                    static void signal_hessian(
                            const Fibre::Strand::Section& section,
                            Fibre::Strand::BasicSection::Tensor& hessian,
                            double diffusion_weighting, double interpolation_weighting,
                            const Triple<double>& diffusion_gradient,
                            const Fibre::Strand::BasicSection& interpolation_gradient,
                            const Coord::Tensor& diffusion_hessian,
                            const Fibre::Strand::BasicSection::Tensor& interpolation_hessian);

                    static void signal_hessian(
                            const Fibre::Tractlet::Section& section,
                            Fibre::Strand::BasicSection::Tensor& hessian,
                            double diffusion_weighting, double interpolation_weighting,
                            const Triple<double>& diffusion_gradient,
                            const Fibre::Strand::BasicSection& interpolation_gradient,
                            const Coord::Tensor& diffusion_hessian,
                            const Fibre::Strand::BasicSection::Tensor& interpolation_hessian);

                    friend std::ostream& operator<<(std::ostream& stream,
                                                    const Direction& direction);
                    
            };
        
        }
    
    }

}

#endif
