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

#ifndef __bts_image_expected_buffer_h__
#define __bts_image_expected_buffer_h__

namespace BTS {
    
    namespace Image {
        
        namespace Expected {
            
            class Buffer;
        
        }
    
    }

}

#include "bts/common.h"

#include "bts/image/observed/buffer.h"
#include "bts/image/index.h"
#include "bts/image/container/buffer.h"
#include "bts/image/reference/buffer.h"

#include "bts/fibre/strand.h"
#include "bts/fibre/tractlet.h"

#include "bts/diffusion/model.h"

#define IMAGE_PARAMETERS \
\
  Option ("img_dims", "The number of voxels along each dimension") \
   + Argument ("img_dims", "").type_text ("[3,3,3]"), \
\
  Option ("img_vox_lengths", "The length of the voxels along each dimension") \
   + Argument ("img_vox_lengths", "").type_text ("[0.15,0.15,0.15]"), \
\
  Option ("img_offsets", "The offset of the centre of the image from the origin (0,0,0).") \
   + Argument ("img_offset", "").type_text ("auto") \


//Loads the parameters into variables
#define SET_IMAGE_PARAMETERS \
\
  Triple<size_t>    img_dims          = Image::Expected::Buffer::DIMS_DEFAULT; \
  Triple<double>    img_vox_lengths   = Image::Expected::Buffer::VOX_LENGTHS_DEFAULT; \
  Triple<double>    img_offsets       = Triple<double>::Invalid; \
\
  Options img_opt = get_options("img_dims"); \
  if (img_opt.size()) \
    img_dims = parse_triple<size_t>(std::string(img_opt[0][0]));  \
\
  img_opt = get_options("img_vox_lengths"); \
  if (img_opt.size()) \
    img_vox_lengths = parse_triple<double>(std::string(img_opt[0][0]));  \
\
  img_opt = get_options("img_offsets"); \
  if (img_opt.size()) \
    img_offsets = parse_triple<double>(std::string(img_opt[0][0]));

//  Triple<size_t>  img_dims (img_dim);
//  Triple<double> img_vox_lengths  (img_vox_length);

//Adds the  parameters to the properties to be saved with the data.
#define ADD_IMAGE_PROPERTIES(properties) \
  properties["img_dims"]            = str(img_dims); \
  properties["img_vox_lengths"]     = str(img_vox_lengths); \
  properties["img_offsets"]         = str(img_offsets); \


//Defines the parameters required to initialise a object.
#define EXPECTED_IMAGE_PARAMETERS \
\
  Option ("exp_num_length_sections", "The number of samples along each strand path when generating the corresponding image set.") \
    + Argument ("exp_num_length_sections", "").type_integer (1, Image::Expected::Buffer::NUM_LENGTH_SECTIONS_DEFAULT, LARGE_INT), \
\
  Option ("exp_num_width_sections", "The number of samples across the tractlet width axes when generating the expected image.") \
    + Argument ("exp_num_width_sections", "").type_integer (1, Image::Expected::Buffer::NUM_WIDTH_SECTIONS_DEFAULT, LARGE_INT), \
\
  Option ("exp_type", "Form of the kernel that will be used to interpolate the signal from the tractlets") \
   + Argument ("exp_type", "").type_text (Image::Expected::Buffer::TYPE_DEFAULT).type_text(), \
\
  Option ("exp_interp_extent", "Extent of the kernel that will be used to interpolate the signal from the tractlets") \
   + Argument ("exp_interp_extent", "").type_float (SMALL_FLOAT, Image::Expected::Buffer::INTERP_EXTENT_DEFAULT, LARGE_FLOAT), \
\
  Option ("exp_enforce_bounds", "Only record the signal (and include in likelihood calculations) that falls withinness the bounds of the image."), \
\
  Option ("exp_half_width", "'Variance' of the kernel that will be used to interpolate the signal from the tractlets when using the Gaussian type kernel.") \
   + Argument ("exp_half_width", "").type_float (SMALL_FLOAT, Image::Expected::Buffer::HALF_WIDTH_DEFAULT, LARGE_FLOAT), \
\
  Option ("exp_base_intensity", "The reference b0 for a \"full\" voxel with unity density. This is used to set the base intensity of the strands. If set to zero (the default) the existing base_intensity of the strands will be used instead.") \
   + Argument ("exp_base_intensity", "").type_float (0.0, 0.0, LARGE_FLOAT), \
\
  Option ("exp_untie_width_intensity", "When not set, intensity will be coupled to the average cross-sectional area of the tract.")

//Loads the parameters into variables
#define SET_EXPECTED_IMAGE_PARAMETERS \
\
  size_t        exp_num_length_sections = Image::Expected::Buffer::NUM_LENGTH_SECTIONS_DEFAULT; \
  size_t        exp_num_width_sections  = Image::Expected::Buffer::NUM_WIDTH_SECTIONS_DEFAULT; \
  std::string   exp_type                = Image::Expected::Buffer::TYPE_DEFAULT; \
  double        exp_interp_extent       = Image::Expected::Buffer::INTERP_EXTENT_DEFAULT; \
  bool          exp_enforce_bounds      = Image::Expected::Buffer::ENFORCE_BOUNDS_DEFAULT; \
  double        exp_half_width          = Image::Expected::Buffer::HALF_WIDTH_DEFAULT; \
  double        exp_base_intensity      = 0.0; \
\
  Options exp_opt = get_options("exp_num_length_sections"); \
  if (exp_opt.size()) \
    exp_num_length_sections = exp_opt[0][0]; \
\
exp_opt = get_options("exp_num_width_sections"); \
  if (exp_opt.size()) \
    exp_num_width_sections = exp_opt[0][0]; \
\
  exp_opt = get_options("exp_type"); \
  if (exp_opt.size()) \
    exp_type = exp_opt[0][0].c_str(); \
\
  exp_opt = get_options("exp_interp_extent"); \
  if (exp_opt.size()) \
    exp_interp_extent = exp_opt[0][0]; \
\
  exp_opt = get_options("exp_enforce_bounds"); \
  if (exp_opt.size()) \
    exp_enforce_bounds = true; \
\
  exp_opt = get_options("exp_half_width"); \
  if (exp_opt.size()) \
    exp_half_width = exp_opt[0][0]; \
\
  exp_opt = get_options("exp_base_intensity"); \
  if (exp_opt.size()) \
    exp_base_intensity = exp_opt[0][0]; \
\

//Adds the parameters to the properties to be saved with the data.
#define ADD_EXPECTED_IMAGE_PROPERTIES(properties) \
  properties["exp_num_length_sections"]    = str(exp_num_length_sections); \
  properties["exp_num_width_sections"]     = str(exp_num_width_sections); \
  properties["exp_interp_extent"]          = str(exp_interp_extent); \
  properties["exp_enforce_bounds"]         = str(exp_enforce_bounds); \
  properties["exp_type"]                   = exp_type; \
  properties["exp_base_intensity"]         = str(exp_base_intensity); \
  if (exp_type == "gaussian") { \
    properties["exp_half_width"]       = exp_half_width; \
  } \


#include "math/matrix.h"

#include "bts/common.h"

#include "bts/image/index.h"
#include "bts/image/buffer.h"
#include "bts/fibre/strand.h"
#include "bts/fibre/tractlet.h"

#include "bts/diffusion/model.h"

#include "bts/image/observed/buffer.h"
#include "bts/image/container/buffer.h"

#include "bts/image/expected/voxel.h"

namespace BTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Gaussian {
                class Buffer;
                class Voxel;
            }
            namespace Quartic {
                class Buffer;
                class Voxel;
            }
            namespace Trilinear {
                class Buffer;
                class Voxel;
            }
            namespace TopHat {
                class Buffer;
                class Voxel;
            }
            namespace Sinc {
                class Buffer;
                class Voxel;
            }
            namespace ReverseSqrt {
                class Buffer;
                class Voxel;
            }
            namespace Realistic {
                class Buffer;
                class Voxel;
            }
            
            class Buffer {
                    
                public:
                    
                    typedef std::map<Index, Voxel>::iterator iterator;
                    typedef std::map<Index, Voxel>::const_iterator const_iterator;

                public:
                    
                    const static size_t NUM_LENGTH_SECTIONS_DEFAULT;
                    const static size_t NUM_WIDTH_SECTIONS_DEFAULT;

                    const static Triple<size_t> DIMS_DEFAULT;
                    const static Triple<double> VOX_LENGTHS_DEFAULT;
                    const static Triple<double> OFFSETS_DEFAULT;
                    const static std::string FILE_EXTENSION;

                    const static double INTERP_EXTENT_DEFAULT;
                    const static bool ENFORCE_BOUNDS_DEFAULT;
                    const static double HALF_WIDTH_DEFAULT;
                    const static char* TYPE_DEFAULT;

                    const static std::string STRAND_BASE_INTENSITY_REFERENCE;
                    const static std::string TRACTLET_BASE_INTENSITY_REFERENCE;

                    //Public static methods.
                public:
                    
                    static Buffer* factory(const std::string& type, const Triple<size_t>& dims,
                                           const Triple<double>& vox_lengths,
                                           const Diffusion::Model& diffusion_model,
                                           size_t num_length_sections, size_t num_width_sections,
                                           double interp_extent, const Triple<double>& offsets,
                                           bool enforce_bounds, double gaussian_half_width);

                    static Buffer* factory(const std::string& type,
                                           const Observed::Buffer& obs_image,
                                           const Diffusion::Model& diffusion_model,
                                           size_t num_length_sections, size_t num_width_sections,
                                           double interp_extent, bool enforce_bounds,
                                           double gaussian_half_width)

                                           {
                        return factory(type, obs_image.dims(), obs_image.vox_lengths(),
                                diffusion_model, num_length_sections, num_width_sections,
                                interp_extent, obs_image.offsets(), enforce_bounds,
                                gaussian_half_width);
                    }
                    
                    //Used for pretty printing in gdb. Is set in the constructor of derived classes.
                protected:
                    
                    char type[50];

                public:
                    
                    virtual ~Buffer() {}

                    virtual Voxel& operator()(int x, int y, int z) = 0;

                    virtual const Voxel& operator()(int x, int y, int z) const = 0;

                    virtual Voxel& operator()(Index c) = 0;

                    virtual const Voxel& operator()(Index c) const = 0;

                    virtual const Triple<size_t>& dims() const = 0;

                    virtual const Triple<double>& vox_lengths() const = 0;

                    //! Offset of the lowest corner of the lowest voxel coordinate ([0,0,0]) of the image.
                    virtual const Triple<double>& offsets() const = 0;

                    virtual size_t dim(size_t dim_index) const = 0;

                    virtual double vox_length(size_t dim_index) const = 0;

                    //! Offset of the lowest corner of the lowest voxel coordinate ([0,0,0]) of the image.
                    virtual double offset(size_t dim_index) const = 0;

                    virtual size_t num_encodings() const = 0;

                    virtual const Diffusion::Encoding& encoding(size_t index) const = 0;

                    virtual size_t num_length_sections() const = 0;

                    virtual size_t num_width_sections() const = 0;

                    virtual void save(const std::string& location) const = 0;

                    virtual Buffer& expected_image(const Fibre::Strand::Set& strands) = 0;

                    virtual Buffer& part_image(const Fibre::Strand& strand) = 0;

                    virtual Buffer
                    & expected_image(const Fibre::Tractlet::Set& tractlets) = 0;

                    virtual Reference::Buffer<Fibre::Strand::Section>::Set
                    & expected_image_with_references(const Fibre::Strand::Set& fibres) = 0;

                    virtual Reference::Buffer<Fibre::Tractlet::Section>::Set
                    & expected_image_with_references(const Fibre::Tractlet::Set& fibres) = 0;

                    virtual Buffer& expected_image(
                            const Fibre::Strand::Set& strands,
                            Container::Buffer<Fibre::Strand>::Set& gradients) = 0;

                    virtual Buffer& expected_image(
                            const Fibre::Tractlet::Set& tractlets,
                            Container::Buffer<Fibre::Tractlet>::Set& gradients) = 0;

                    virtual Buffer& expected_image(
                            const Fibre::Strand::Set& strands,
                            Container::Buffer<Fibre::Strand>::Set& gradients,
                            Container::Buffer<Fibre::Strand::Tensor>::Set& hessians) = 0;

                    virtual Buffer& expected_image(
                            const Fibre::Tractlet::Set& tractlets,
                            Container::Buffer<Fibre::Tractlet>::Set& gradients,
                            Container::Buffer<Fibre::Tractlet::Tensor>::Set& hessians) = 0;

                    virtual void precalculate_section_weighting_gradients(
                            const Fibre::Strand& dummy) = 0;

                    virtual void precalculate_section_weighting_gradients_and_hessians(
                            const Fibre::Strand& dummy) = 0;

                    virtual void precalculate_section_weighting_gradients(
                            const Fibre::Tractlet& dummy) = 0;

                    virtual void precalculate_section_weighting_gradients_and_hessians(
                            const Fibre::Tractlet& dummy) = 0;

                    virtual Buffer* clone() const = 0;

                    virtual const Diffusion::Model& get_diffusion_model() const = 0;

                    virtual std::ostream& to_stream(std::ostream& stream) const = 0;

                    Image::Buffer clean_buffer();

                    virtual double get_base_intensity(double ref_b0) = 0;

                    bool dims_match(const Observed::Buffer& reference);

                    virtual Properties& properties() = 0;

                    virtual const Properties& properties() const = 0;

                    virtual double rms(bool include_b0s) const = 0;

                    virtual double max_b0() const = 0;

                    virtual Buffer& zero() = 0;

                    virtual Buffer& clear() = 0;

                    virtual Buffer& negate() = 0;

                    virtual std::set<Index> non_empty() const = 0;

                    virtual std::set<Index> non_empty_or_inbounds() const = 0;

                    virtual std::set<Index> empty_inbounds() const = 0;

                    virtual bool bounds_are_enforced() const = 0;

                    virtual void clear_and_enforce_bounds() = 0;

                    virtual void relax_bounds() = 0;

                    virtual Buffer& operator+=(const Buffer& buffer) = 0;

                    virtual Buffer& operator-=(const Buffer& buffer) = 0;

                    virtual Buffer& operator-=(
                            const Image::Buffer_tpl<Observed::Voxel>& buffer) = 0;

                    virtual iterator begin();
                    
                    virtual iterator end();
            };
            
#define EXPECTED_BUFFER_FUNCTIONS \
 \
            void                        save(const std::string& location) const \
              { this->Observed::Buffer_tpl<Voxel>::save(location); } \
 \
            Buffer&                     expected_image(const Fibre::Strand::Set& strands) \
              { this->Buffer_tpl<Voxel>::expected_image<Fibre::Strand>(strands); return *this; } \
\
            Buffer&                     part_image(const Fibre::Strand& strand) \
              { this->Buffer_tpl<Voxel>::part_image<Fibre::Strand>(strand); return *this; } \
 \
            Buffer&                     expected_image(const Fibre::Tractlet::Set& tractlets) \
              { this->Buffer_tpl<Voxel>::expected_image<Fibre::Tractlet>(tractlets); return *this; } \
 \
            Reference::Buffer<Fibre::Strand::Section>::Set&   expected_image_with_references(const Fibre::Strand::Set& fibres) \
              { return this->Buffer_tpl<Voxel>::expected_image_with_references<Fibre::Strand>(fibres); } \
 \
            Reference::Buffer<Fibre::Tractlet::Section>::Set& expected_image_with_references(const Fibre::Tractlet::Set& fibres) \
              { return this->Buffer_tpl<Voxel>::expected_image_with_references<Fibre::Tractlet>(fibres); } \
 \
            Buffer&                     expected_image(const Fibre::Strand::Set& strands, Container::Buffer<Fibre::Strand>::Set& gradients) \
              { this->Buffer_tpl<Voxel>::expected_image<Fibre::Strand>(strands, gradients); return *this; } \
 \
            Buffer&                     expected_image(const Fibre::Tractlet::Set& tractlets, Container::Buffer<Fibre::Tractlet>::Set& gradients) \
              { this->Buffer_tpl<Voxel>::expected_image<Fibre::Tractlet>(tractlets, gradients); return *this; } \
 \
            Buffer&                     expected_image(const Fibre::Strand::Set& strands, Container::Buffer<Fibre::Strand>::Set& gradients, Container::Buffer<Fibre::Strand::Tensor>::Set& hessians) \
              { this->Buffer_tpl<Voxel>::expected_image<Fibre::Strand>(strands, gradients, hessians); return *this; } \
 \
            Buffer&                     expected_image(const Fibre::Tractlet::Set& tractlets, Container::Buffer<Fibre::Tractlet>::Set& gradients, Container::Buffer<Fibre::Tractlet::Tensor>::Set& hessians) \
              { this->Buffer_tpl<Voxel>::expected_image<Fibre::Tractlet>(tractlets, gradients, hessians); return *this; } \
 \
            Buffer&                     part_image (const Fibre::Strand& strand, Container::Buffer<Fibre::Strand>& gradients, Container::Buffer<Fibre::Strand::Tensor>& hessians) \
              { this->Buffer_tpl<Voxel>::part_image(strand, gradients, hessians); return *this; } \
 \
            Buffer&                     part_image (const Fibre::Tractlet& tractlet, Container::Buffer<Fibre::Tractlet>& gradients, Container::Buffer<Fibre::Tractlet::Tensor>& hessians) \
              { this->Buffer_tpl<Voxel>::part_image(tractlet, gradients, hessians); return *this; } \
 \
            Buffer&                     part_image(const Fibre::Strand& strand, Container::Buffer<Fibre::Strand>& gradients) \
              { this->Buffer_tpl<Voxel>::part_image(strand, gradients); return *this; } \
 \
            Buffer&                     part_image(const Fibre::Tractlet& tractlet, Container::Buffer<Fibre::Tractlet>& gradients) \
              { this->Buffer_tpl<Voxel>::part_image(tractlet, gradients); return *this; } \
 \
            void                        precalculate_section_weighting_gradients(const Fibre::Strand& dummy) \
              { Expected::Buffer_tpl<Voxel>::precalculate_section_weighting_gradients<Fibre::Strand>(); } \
 \
            void                        precalculate_section_weighting_gradients_and_hessians(const Fibre::Strand& dummy) \
              { Expected::Buffer_tpl<Voxel>::precalculate_section_weighting_gradients_and_hessians<Fibre::Strand>(); } \
 \
            void                        precalculate_section_weighting_gradients(const Fibre::Tractlet& dummy) \
              { Expected::Buffer_tpl<Voxel>::precalculate_section_weighting_gradients<Fibre::Tractlet>(); } \
 \
            void                        precalculate_section_weighting_gradients_and_hessians(const Fibre::Tractlet& dummy) \
              { Expected::Buffer_tpl<Voxel>::precalculate_section_weighting_gradients_and_hessians<Fibre::Tractlet>(); } \
 \
            Buffer*                     clone() const \
              { return new Buffer(*this); } \
 \
            const Diffusion::Model&       get_diffusion_model() const \
              { return this->Buffer_tpl<Voxel>::get_diffusion_model(); } \
 \
            double                      get_base_intensity(double ref_b0) \
             { return this->Buffer_tpl<Voxel>::get_base_intensity(ref_b0); } \
            size_t                          num_length_sections() const \
              { return this->Buffer_tpl<Voxel>::num_length_sections(); } \
 \
            size_t                          num_width_sections() const \
              { return this->Buffer_tpl<Voxel>::num_width_sections(); } \
 \
            Voxel&                        operator() (Index coord) \
              { return Image::Buffer_tpl<Voxel>::operator()(coord); } \
 \
            const Voxel&                  operator() (Index coord) const \
              { return Image::Buffer_tpl<Voxel>::operator()(coord); } \
 \
            Voxel&                        operator() (int x, int y, int z) \
              { return Image::Buffer_tpl<Voxel>::operator()(x,y,z); } \
 \
            const Voxel&                  operator() (int x, int y, int z) const \
              { return Image::Buffer_tpl<Voxel>::operator()(x,y,z); } \
 \
            size_t                          num_encodings() const \
              { return this->Buffer_tpl<Voxel>::num_encodings(); } \
 \
            const Diffusion::Encoding&    encoding(size_t index) const \
              { return this->Buffer_tpl<Voxel>::encoding(index); } \
 \
            size_t                          dim(size_t dim_index) const \
              { return this->Image::Buffer_tpl<Voxel>::dim(dim_index); } \
 \
            double                        vox_length(size_t dim_index) const \
              { return this->Buffer_tpl<Voxel>::vox_length(dim_index); } \
 \
            double                        offset(size_t dim_index) const \
              { return this->Buffer_tpl<Voxel>::offset(dim_index); } \
 \
            const Triple<size_t>&           dims() const \
              { return this->Image::Buffer_tpl<Voxel>::dims(); } \
           \
            const Triple<double>&         vox_lengths() const \
              { return this->Buffer_tpl<Voxel>::vox_lengths(); } \
 \
            const Triple<double>&         offsets() const \
              { return this->Buffer_tpl<Voxel>::offsets(); } \
 \
            Buffer&                       zero () \
              { this->Buffer_tpl<Voxel>::zero(); return *this; } \
 \
            Buffer&                       clear () \
              { this->Buffer_tpl<Voxel>::clear(); return *this; } \
\
            Buffer&                       negate () \
             { this->Buffer_tpl<Voxel>::negate(); return *this; } \
 \
            Buffer&                       operator+= (const Buffer& buff) \
              { this->Buffer_tpl<Voxel>::operator+= (buff); return *this; } \
 \
            Buffer&                       operator-= (const Buffer& buff) \
              { this->Buffer_tpl<Voxel>::operator-= (buff); return *this; } \
\
            Buffer&                       operator-= (const Image::Buffer_tpl<Observed::Voxel>& buff) \
              { this->Buffer_tpl<Voxel>::operator-= (buff); return *this; } \
 \
            Buffer&                       operator*= (double M) \
              { this->Buffer_tpl<Voxel>::operator*= (M); return *this; } \
 \
            Buffer&                       operator/= (double M) \
              { this->Buffer_tpl<Voxel>::operator/= (M); return *this; } \
 \
            template <typename T> Buffer& operator+= (const Image::Buffer_tpl<T>& buff) \
              { this->Buffer_tpl<Voxel>::operator+= (buff); return *this; } \
 \
            template <typename T> Buffer& operator-= (const Image::Buffer_tpl<T>& buff) \
              { this->Buffer_tpl<Voxel>::operator-= (buff); return *this; } \
\
           Buffer&                       operator+= (const Expected::Buffer& buff) \
             { this->Buffer_tpl<Voxel>::operator+= (buff); return *this; } \
\
           Buffer&                       operator-= (const Expected::Buffer& buff) \
             { this->Buffer_tpl<Voxel>::operator-= (buff); return *this; } \
\
            std::ostream&                 to_stream (std::ostream& stream) const \
              { return this->Buffer_tpl<Voxel>::to_stream(stream); } \
 \
            Properties&                   properties() \
              { return this->Buffer_tpl<Voxel>::properties(); } \
 \
            const Properties&             properties() const \
              { return this->Buffer_tpl<Voxel>::properties(); } \
 \
            double                        rms(bool include_b0s) const \
              { return Buffer_tpl<Voxel>::rms(include_b0s); } \
 \
            double                        max_b0() const \
              { return Buffer_tpl<Voxel>::max_b0(); } \
 \
            std::set<Index>               non_empty() const \
              { return Image::Buffer_tpl<Voxel>::non_empty(); } \
 \
            std::set<Index>               non_empty_or_inbounds() const \
              { return Image::Buffer_tpl<Voxel>::non_empty_or_inbounds(); } \
\
            std::set<Index>               empty_inbounds() const \
             { return Image::Buffer_tpl<Voxel>::empty_inbounds(); } \
\
            bool                          bounds_are_enforced() const \
              { return Image::Buffer_tpl<Voxel>::bounds_are_enforced(); } \
 \
            void                          clear_and_enforce_bounds() \
              { return Image::Buffer_tpl<Voxel>::clear_and_enforce_bounds(); } \
\
            void                          relax_bounds() \
              { return Image::Buffer_tpl<Voxel>::relax_bounds(); } \
\
          protected: \
\
            virtual Voxel       new_voxel(const Index& coord) \
              { return Voxel(*this,coord); } \
\
            void                          name_init() \
              { strncpy (this->type, SHORT_NAME.c_str(), 50); } \
\
          friend class Voxel
            
            inline std::ostream& operator<<(std::ostream& stream, const Buffer& image) {
                
                return image.to_stream(stream);
                
            }
        
        }
    
    }

}

#endif
