/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Bayesian Tractlet Sampling (BTS).

 BTS is free software: you can reimageibute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 BTS is imageibuted in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.

 */


#ifndef __bts_analysis_imagehessiantester_h__
#define __bts_analysis_imagehessiantester_h__

namespace BTS {

  namespace Analysis {

    template <typename ImageClass, typename State> class ImageHessianTester {


    public:

      //typedef the pointer to the function to test as 'Function'.
      typedef ImageClass& (ImageClass::*Function)(const State&, Image::Container::Buffer<State>&, Image::Container::Buffer<typename State::Tensor>&); //, typename Image::Container::Buffer<typename State::Tensor>&

    protected:

      ImageClass                      image;

    public:

      ImageHessianTester(const ImageClass& image)
        : image(image) {  this->image.clear_and_enforce_bounds(); }

      ~ImageHessianTester() {}

      void                                    test(Function function, State& state, double step_size, Image::Container::Buffer<typename State::Tensor>& analytic_hessian, Image::Container::Buffer<typename State::Tensor>& numeric_hessian);


  };

 }

}


#include "progressbar.h"

#include "bts/image/container/buffer.h"

#include "bts/analysis/image_hessian_tester.h"

namespace BTS {

  namespace Analysis {

    template <typename ImageClass, typename State> void     ImageHessianTester<ImageClass, State>::test(Function function, State& state, double step_size, Image::Container::Buffer<typename State::Tensor>& analytic_hessian, Image::Container::Buffer<typename State::Tensor>& numeric_hessian) {

      Triple<size_t> dims = analytic_hessian.dims();
      size_t num_encodings = analytic_hessian.num_encodings();

      analytic_hessian.clear();
      numeric_hessian.clear();

      Image::Container::Buffer<typename State::Tensor> dummy_hessian (analytic_hessian);
      Image::Container::Buffer<State> gradient(dims,num_encodings), perturbed_gradient(dims, num_encodings);

      //Calculate analytic hessian and unperturbed gradient

      gradient.zero();

      (image.*function)(state, gradient, analytic_hessian);

      //Initialize


      MR::Math::Vector<double>& state_vector = state;

      typename State::Tensor state_tensor(state);
      MR::Math::Matrix<double> blank_tensor_matrix(state_vector.size(), state_vector.size());
      blank_tensor_matrix = 0.0;

      Image::Container::Buffer< MR::Math::Matrix<double> > numeric_hessian_matrices (dims, num_encodings);

      for (size_t z = 0; z < image.dim(Z); z++)
        for (size_t y = 0; y < image.dim(Y); y++)
          for (size_t x = 0; x < image.dim(X); x++)
            for (size_t encode_i = 0; encode_i < num_encodings; encode_i++)
              numeric_hessian_matrices(x,y,z)[encode_i] = blank_tensor_matrix;


      MR::ProgressBar progress_bar ("Testing hessian calculations...", state_vector.size());

      State perturbed_state = state;

      for (size_t row_i = 0; row_i < state_vector.size(); ++row_i) {

        state_vector[row_i] += step_size;

        perturbed_gradient.zero();

        (image.*function)(state, perturbed_gradient, dummy_hessian);

        for (typename Image::Container::Buffer<State>::iterator vox_it = perturbed_gradient.begin(); vox_it != perturbed_gradient.end(); ++vox_it)
          for (size_t encode_i = 0; encode_i < num_encodings; encode_i++)
            numeric_hessian_matrices(vox_it->first)[encode_i].row(row_i) = ((vox_it->second[encode_i] - gradient(vox_it->first)[encode_i]) / step_size);


        state_vector[row_i] -= step_size;

        ++progress_bar;

      }


      for (size_t z = 0; z < image.dim(Z); z++)
        for (size_t y = 0; y < image.dim(Y); y++)
          for (size_t x = 0; x < image.dim(X); x++) {

            for (size_t encode_i = 0; encode_i < num_encodings; encode_i++) {
              numeric_hessian(x,y,z)[encode_i] = state_tensor;
              numeric_hessian(x,y,z)[encode_i] = numeric_hessian_matrices(x,y,z)[encode_i];
            }

            if (analytic_hessian.is_empty(x,y,z))
              for (size_t encode_i = 0; encode_i < image.num_encodings(); encode_i++)
                analytic_hessian(x,y,z)[encode_i] = typename State::Tensor(state_tensor).zero();

          }

    }


  }
}

#undef LOOP

#endif
