/*
    Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

    Written by Tom G. Close, 4/03/10.

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


#include "bts/cmd.h"

#include "bts/common.h"

#include "math/matrix.h"

#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood/gaussian.h"
#include "bts/prob/likelihood/one_sided_gaussian.h"
#include "bts/prob/test/gaussian.h"
#include "bts/prob/test/landscape.h"
#include "bts/prob/test/landscape/peak.h"
#include "bts/diffusion/encoding.h"

#include "bts/analysis/gradient_tester.h"
#include "bts/analysis/image_gradient_tester.h"
#include "bts/analysis/hessian_tester.h"
#include "bts/analysis/image_hessian_tester.h"
#include "bts/analysis/rank3_hessian_tester.h"
#include "bts/analysis/fisher_gradient_tester.h"
#include "bts/analysis/fisher_gradient_matrix_tester.h"

#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"

#include "bts/diffusion/model.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/quartic/buffer.h"

#include "bts/image/observed/buffer.h"

#include "bts/prob/uniform.h"
#include "bts/prob/prior.h"





#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood/one_sided_gaussian.h"
#include "bts/prob/likelihood/gaussian.h"

#include "bts/mcmc/state.h"
#include "bts/utilities/reader.h"
#include "bts/mcmc/state/tensor.h"
#include "bts/mcmc/state/tensor/writer.h"

#include "bts/fibre/strand/basic_section/tensor.h"
#include "bts/fibre/tractlet/section/tensor.h"
#include "bts/math/common.h"

#include "bts/fibre/base/tensor_writer.cpp.h"


#include "bts/inline_functions.h"

#define TEST_GRADIENT(State, Object, function, object_instance) \
  State state; \
\
  State::Reader reader(state_location); \
  State::Writer analytic_writer (analytic_output_location, reader, properties), numeric_writer (numeric_output_location, reader, properties); \
\
  Analysis::GradientTester<Object, State> grad_tester (object_instance); \
  Analysis::GradientTester<Object, State>::Function test_function; \
\
  test_function = &Object::function; \
\
  while (reader.next(state)) { \
    State analytic_gradient(state), numeric_gradient(state); \
\
    grad_tester.test(test_function, state, step_size, analytic_gradient, numeric_gradient); \
\
    analytic_writer.append(analytic_gradient); \
    numeric_writer.append(numeric_gradient); \
\
  } \
\
  reader.close(); \
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_ELEMENT_GRADIENT(State, Object, function, object_instance) \
  State state; \
\
  State::Reader reader(state_location); \
  State::Writer analytic_writer (analytic_output_location, reader, properties), numeric_writer (numeric_output_location, reader, properties); \
\
  Analysis::GradientTester<Object, State::Element> grad_tester (object_instance); \
  Analysis::GradientTester<Object, State::Element>::Function test_function; \
\
  test_function = &Object::function; \
\
  while (reader.next(state)) { \
    State analytic_gradient(state.size(), numeric_gradient(state.size()); \
\
\   for (size_t elem_i = 0; elem_i < state.size(); elem_i++) { \
      grad_tester.test(test_function, state[elem_i], step_size, analytic_gradient[elem_i], numeric_gradient[elem_i]); \
    } \
\
    analytic_writer.append(analytic_gradient); \
    numeric_writer.append(numeric_gradient); \
\
  } \
\
  reader.close(); \
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_SINGLE_GRADIENT(State, Object, function, object_instance) \
\
  State state (state_location); \
  State analytic_gradient(state), numeric_gradient(state); \
\
  Analysis::GradientTester<Object, State> grad_tester (object_instance); \
  Analysis::GradientTester<Object, State>::Function test_function; \
\
  test_function = &Object::function; \
\
  grad_tester.test(test_function, state, step_size, analytic_gradient, numeric_gradient); \
\
  analytic_gradient.save(analytic_output_location, properties); \
  numeric_gradient.save(numeric_output_location, properties); \


#define TEST_IMAGE_GRADIENT(State, Object, function, object_instance) \
  State state; \
\
  State::Reader reader(state_location); \
  State::Writer analytic_writer (analytic_output_location, reader, properties), numeric_writer (numeric_output_location, reader, properties); \
\
  Analysis::ImageGradientTester<Object, State> grad_tester (object_instance); \
  Analysis::ImageGradientTester<Object, State>::Function test_function; \
\
  test_function = &Object::function; \
\
  while (reader.next(state)) { \
\
    object_instance.zero(); \
\
    Image::Container::Buffer<State> analytic_gradient(img_dims, object_instance.num_encodings()); \
    Image::Container::Buffer<State> numeric_gradient(img_dims, object_instance.num_encodings()); \
\
    grad_tester.test(test_function, state, step_size, analytic_gradient, numeric_gradient); \
\
    for (size_t z = 0; z < object_instance.dim(Z); z++) { \
      for (size_t y = 0; y < object_instance.dim(Y); y++) { \
        for (size_t x = 0; x < object_instance.dim(X); x++) { \
          for (size_t encode_i = 0; encode_i < object_instance.num_encodings(); encode_i++) { \
            analytic_writer.append(analytic_gradient(x,y,z)[encode_i]); \
            numeric_writer.append(numeric_gradient(x,y,z)[encode_i]); \
          } \
        } \
      } \
    } \
  } \
\
  reader.close(); \
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_IMAGE_HESSIAN(State, Object, function, object_instance) \
  State state; \
  State::Reader reader (state_location); \
  if (!reader.next(state)) \
    throw Exception ("Did not find any states in file '" + state_location + "'."); \
  reader.close(); \
\
  State::Tensor::Writer analytic_writer (analytic_output_location + ".tnr", state), numeric_writer (numeric_output_location + ".tnr", state); \
  reader.open(state_location); \
\
  Analysis::ImageHessianTester<Object, State> hess_tester (object_instance); \
  Analysis::ImageHessianTester<Object, State>::Function test_function; \
\
  test_function = &Object::function; \
\
  while (reader.next(state)) { \
\
    object_instance.zero(); \
\
    Image::Container::Buffer<State::Tensor> analytic_hessian(img_dims, object_instance.num_encodings()); \
    Image::Container::Buffer<State::Tensor> numeric_hessian(img_dims, object_instance.num_encodings()); \
\
    hess_tester.test(test_function, state, step_size, analytic_hessian, numeric_hessian); \
\
    for (size_t z = 0; z < object_instance.dim(Z); z++) { \
      for (size_t y = 0; y < object_instance.dim(Y); y++) { \
        for (size_t x = 0; x < object_instance.dim(X); x++) { \
          for (size_t encode_i = 0; encode_i < object_instance.num_encodings(); encode_i++) { \
            analytic_writer.append(analytic_hessian(x,y,z)[encode_i]); \
            numeric_writer.append(numeric_hessian(x,y,z)[encode_i]); \
          } \
        } \
      } \
    } \
  } \
\
  reader.close(); \
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_OBJECT_ELEMENT_GRADIENT(State, Object, function, object_instance) \
  State state; \
\
  State::Writer analytic_writer (analytic_output_location, properties), numeric_writer (numeric_output_location, properties); \
  State::Reader reader; \
\
  for (size_t i = 0; i < object_instance.size(); ++i) { \
\
    reader.open(state_location); \
\
    Analysis::GradientTester<Object, State> grad_tester (object_instance[i]); \
    Analysis::GradientTester<Object, State>::Function test_function; \
\
    test_function = &Object::function; \
\
    while (reader.next(state)) { \
      State analytic_gradient(state), numeric_gradient(state); \
\
      grad_tester.test(test_function, state, step_size, analytic_gradient, numeric_gradient); \
\
      analytic_writer.append(analytic_gradient); \
      numeric_writer.append(numeric_gradient); \
\
    } \
\
    reader.close(); \
\
  } \
\
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_OBJECT_ELEMENT_HESSIAN(State, Object, function, object_instance) \
  State state; \
  State::Reader reader (state_location); \
  if (!reader.next(state)) \
    throw Exception ("Did not find any states in file '" + state_location + "'."); \
  reader.close(); \
\
  State::Tensor::Writer analytic_writer (analytic_output_location + ".tnr", state), numeric_writer (numeric_output_location + ".tnr", state); \
\
  for (size_t i = 0; i < object_instance.size(); ++i) { \
\
    reader.open(state_location); \
\
    Analysis::HessianTester<Object, State> grad_tester (object_instance[i]); \
    Analysis::HessianTester<Object, State>::Function test_function; \
\
    test_function = &Object::function; \
\
    while (reader.next(state)) { \
      State::Tensor analytic_hessian(state), numeric_hessian(state); \
\
      grad_tester.test(test_function, state, step_size, analytic_hessian, numeric_hessian); \
\
      analytic_writer.append(analytic_hessian); \
      numeric_writer.append(numeric_hessian); \
  \
    } \
\
    reader.close(); \
\
  } \
\
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_WRITABLE_OBJECT_ELEMENT_GRADIENT(State, Object, function, object_instance) \
  State::Writable state; \
\
  State::Writable::Writer analytic_writer (analytic_output_location, properties), numeric_writer (numeric_output_location, properties); \
  State::Writable::Reader reader(state_location); \
\
  for (size_t i = 0; i < object_instance.size(); ++i) { \
\
    Analysis::GradientTester<Object, State> grad_tester (object_instance[i]); \
    Analysis::GradientTester<Object, State>::Function test_function; \
\
    test_function = &Object::function; \
\
    while (reader.next(state)) { \
      State analytic_gradient(state.size(), numeric_gradient(state.size()); \
\
      grad_tester.test(test_function, state, step_size, analytic_gradient, numeric_gradient); \
\
      analytic_writer.append(analytic_gradient); \
      numeric_writer.append(numeric_gradient); \
\
    } \
  } \
\
  reader.close(); \
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_VOXEL_GRADIENT(State, Object, function, image_instance) \
  State state; \
\
  State::Writer analytic_writer (analytic_output_location, properties), numeric_writer (numeric_output_location, properties); \
\
  for (size_t z = 0; z < image_instance.dim(Z); ++z) { \
    for (size_t y = 0; y < image_instance.dim(Y); ++y) { \
      for (size_t x = 0; x < image_instance.dim(X); ++x) { \
\
        State::Reader reader(state_location); \
\
        Analysis::GradientTester<Object, State> grad_tester (image_instance(x,y,z)); \
        Analysis::GradientTester<Object, State>::Function test_function; \
\
        test_function = &Object::function; \
\
        while (reader.next(state)) { \
          State analytic_gradient(state), numeric_gradient(state); \
\
          grad_tester.test(test_function, state, step_size, analytic_gradient, numeric_gradient); \
\
          analytic_writer.append(analytic_gradient); \
          numeric_writer.append(numeric_gradient); \
\
        } \
\
        reader.close(); \
\
      } \
    } \
  } \
\
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_VOXEL_HESSIAN(State, Object, function, image_instance) \
  State state; \
  State::Reader reader (state_location); \
  if (!reader.next(state)) \
    throw Exception ("Did not find any states in file '" + state_location + "'."); \
  reader.close(); \
\
  State::Tensor::Writer analytic_writer (analytic_output_location + ".tnr", state), numeric_writer (numeric_output_location + ".tnr", state); \
\
  for (size_t z = 0; z < image_instance.dim(Z); ++z) { \
    for (size_t y = 0; y < image_instance.dim(Y); ++y) { \
      for (size_t x = 0; x < image_instance.dim(X); ++x) { \
\
        State::Reader reader(state_location); \
\
        Analysis::HessianTester<Object, State> grad_tester (image_instance(x,y,z)); \
        Analysis::HessianTester<Object, State>::Function test_function; \
\
        test_function = &Object::function; \
\
        while (reader.next(state)) { \
          State::Tensor analytic_hessian(state), numeric_hessian(state); \
\
          grad_tester.test(test_function, state, step_size, analytic_hessian, numeric_hessian); \
\
          analytic_writer.append(analytic_hessian); \
          numeric_writer.append(numeric_hessian); \
\
        } \
\
        reader.close(); \
\
      } \
    } \
  } \
\
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_WRITABLE_VOXEL_GRADIENT(State, Object, function, image_instance) \
  State::Writable state; \
\
  State::Writable::Writer analytic_writer (analytic_output_location, properties), numeric_writer (numeric_output_location, properties); \
\
  for (size_t z = 0; z < image_instance.dim(Z); ++z) { \
    for (size_t y = 0; y < image_instance.dim(Y); ++y) { \
      for (size_t x = 0; x < image_instance.dim(X); ++x) { \
\
        State::Writable::Reader reader(state_location); \
\
        Analysis::GradientTester<Object, State> grad_tester (image_instance(x,y,z)); \
        Analysis::GradientTester<Object, State>::Function test_function; \
\
        test_function = &Object::function; \
\
        while (reader.next(state)) { \
          State analytic_gradient(state.size(), numeric_gradient(state.size()); \
\
          grad_tester.test(test_function, state, step_size, analytic_gradient, numeric_gradient); \
\
          analytic_writer.append(analytic_gradient); \
          numeric_writer.append(numeric_gradient); \
\
        } \
\
        reader.close(); \
\
      } \
    } \
  } \
\
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_DIRECTION_GRADIENT(State, function, image_instance) \
  State state; \
  std::map<std::string,std::string> read_properties; \
  std::vector<std::string> row_properties_list; \
\
  State::Writer analytic_writer (analytic_output_location, properties), numeric_writer (numeric_output_location, properties); \
\
  for (size_t z = 0; z < image_instance.dim(Z); ++z) { \
    for (size_t y = 0; y < image_instance.dim(Y); ++y) { \
      for (size_t x = 0; x < image_instance.dim(X); ++x) { \
        for (size_t encode_i = 0; encode_i < image_instance.num_encodings(); encode_i++) { \
\
          State::Reader reader(state_location); \
\
          Analysis::GradientTester<Image::Expected::Direction, State> grad_tester (image_instance(x,y,z).direction(encode_i)); \
          Analysis::GradientTester<Image::Expected::Direction, State>::Function test_function; \
\
          test_function = &Image::Expected::Direction::function; \
\
          while (reader.next(state)) { \
            State analytic_gradient(state.size(), numeric_gradient(state.size()); \
\
            grad_tester.test(test_function, state, step_size, analytic_gradient, numeric_gradient); \
\
            analytic_writer.append(analytic_gradient); \
            numeric_writer.append(numeric_gradient); \
          } \
\
          reader.close(); \
\
        } \
      } \
    } \
  } \
\
  analytic_writer.close(); \
  numeric_writer.close();


#define TEST_DIRECTION_HESSIAN(State, function, image_instance) \
  State state; \
  State::Reader reader (state_location); \
  if (!reader.next(state)) \
    throw Exception ("Did not find any states in file '" + state_location + "'."); \
  reader.close(); \
\
  State::Tensor::Writer analytic_writer (analytic_output_location + ".tnr", state), numeric_writer (numeric_output_location + ".tnr", state); \
\
  for (size_t z = 0; z < image_instance.dim(Z); ++z) { \
    for (size_t y = 0; y < image_instance.dim(Y); ++y) { \
      for (size_t x = 0; x < image_instance.dim(X); ++x) { \
        for (size_t encode_i = 0; encode_i < image_instance.num_encodings(); encode_i++) { \
\
          State::Reader reader(state_location); \
\
          Analysis::HessianTester<Image::Expected::Direction, State> hess_tester (image_instance(x,y,z).direction(encode_i)); \
          Analysis::HessianTester<Image::Expected::Direction, State>::Function test_function; \
\
          test_function = &Image::Expected::Direction::function; \
\
          while (reader.next(state)) { \
            State::Tensor analytic_hessian(state.size(), numeric_hessian(state.size()); \
\
            hess_tester.test(test_function, state, step_size, analytic_hessian, numeric_hessian); \
\
            analytic_writer.append(analytic_hessian); \
            numeric_writer.append(numeric_hessian); \
          } \
\
          reader.close(); \
\
        } \
      } \
    } \
  } \
\
  analytic_writer.close(); \
  numeric_writer.close();


#define TEST_IMAGE_POINTER_GRADIENT(State, Object, function, object_instance, state) \
  State* analytic_gradient = state.clone(); \
  State* numeric_gradient  = state.clone(); \
\
  Analysis::PointerGradientTester<Object, State> grad_tester (object_instance); \
  Analysis::PointerGradientTester<Object, State>::Function test_function; \
\
  test_function = &Object::function; \
\
  grad_tester.test(test_function, state, step_size, *analytic_gradient, *numeric_gradient); \
\
  analytic_gradient->save(analytic_output_location); \
  numeric_gradient->save(numeric_output_location); \
\
  delete analytic_gradient; \
  delete numeric_gradient; \


#define TEST_HESSIAN(State, Object, function, object_instance) \
  State state; \
  State::Reader reader (state_location); \
  if (!reader.next(state)) \
    throw Exception ("Did not find any states in file '" + state_location + "'."); \
\
  State::Tensor::Writer analytic_writer (analytic_output_location + ".tnr", state), numeric_writer (numeric_output_location + ".tnr", state); \
\
  reader.close(); \
  reader.open(state_location); \
\
  Analysis::HessianTester<Object, State> grad_tester (object_instance); \
  Analysis::HessianTester<Object, State>::Function test_function; \
\
  test_function = &Object::function; \
\
  while (reader.next(state)) { \
    State::Tensor analytic_hessian(state), numeric_hessian(state); \
\
    grad_tester.test(test_function, state, step_size, analytic_hessian, numeric_hessian); \
\
    analytic_writer.append(analytic_hessian); \
    numeric_writer.append(numeric_hessian); \
\
  } \
\
  reader.close(); \
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_RANK3_HESSIAN(State, Object, function, object_instance) \
  State state; \
  State::Reader reader (state_location); \
  if (!reader.next(state)) \
    throw Exception ("Did not find any states in file '" + state_location + "'."); \
\
  State::Tensor::Writer analytic_writer (analytic_output_location + ".tnr3", state), numeric_writer (numeric_output_location + ".tnr3", state); \
\
  reader.close(); \
  reader.open(state_location); \
\
  Analysis::Rank3HessianTester<Object, State> grad_tester (object_instance); \
  Analysis::Rank3HessianTester<Object, State>::Function test_function; \
\
  test_function = &Object::function; \
\
  while (reader.next(state)) { \
    std::vector<State::Tensor> analytic_rank3_hessian, numeric_rank3_hessian; \
\
    grad_tester.test(test_function, state, step_size, analytic_rank3_hessian, numeric_rank3_hessian); \
\
    for (size_t dim_i = 0; dim_i < state.size(); dim_i++) { \
      analytic_writer.append(analytic_rank3_hessian[dim_i]); \
      numeric_writer.append(numeric_rank3_hessian[dim_i]); \
    } \
\
  } \
\
  reader.close(); \
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_FISHER_GRADIENT(State, Object, function, object_instance) \
  State state; \
  State::Reader reader (state_location); \
  if (!reader.next(state)) \
    throw Exception ("Did not find any states in file '" + state_location + "'."); \
\
  State::Tensor::Writer analytic_writer (analytic_output_location + ".tnr", state), numeric_writer (numeric_output_location + ".tnr", state); \
\
  reader.close(); \
  reader.open(state_location); \
\
  Analysis::FisherGradientTester<Object, State> grad_tester (object_instance); \
  Analysis::FisherGradientTester<Object, State>::Function test_function; \
\
  test_function = &Object::function; \
\
  while (reader.next(state)) { \
    std::vector<State::Tensor> analytic_fisher_gradient, numeric_fisher_gradient; \
\
    grad_tester.test(test_function, state, step_size, analytic_fisher_gradient, numeric_fisher_gradient); \
\
    if (analytic_fisher_gradient.size() != numeric_fisher_gradient.size()) \
      throw Exception ("Sizes of Fisher gradients do not match"); \
\
    for (size_t i = 0; i < analytic_fisher_gradient.size(); i++) { \
      analytic_writer.append(analytic_fisher_gradient[i]); \
      numeric_writer.append(numeric_fisher_gradient[i]); \
    } \
\
  } \
\
  reader.close(); \
  analytic_writer.close(); \
  numeric_writer.close(); \


#define TEST_FISHER_GRADIENT_MATRIX(State, Object, function, object_instance) \
  State state; \
  State::Reader reader (state_location); \
  if (!reader.next(state)) \
    throw Exception ("Did not find any states in file '" + state_location + "'."); \
\
  State::Tensor::Writer analytic_writer (analytic_output_location + ".tnr", state), numeric_writer (numeric_output_location + ".tnr", state); \
\
  reader.close(); \
  reader.open(state_location); \
\
  Analysis::FisherGradientMatrixTester<Object, State> grad_tester (object_instance); \
  Analysis::FisherGradientMatrixTester<Object, State>::Function test_function; \
\
  test_function = &Object::function; \
\
  while (reader.next(state)) { \
    std::vector<State::Tensor> analytic_fisher_gradient, numeric_fisher_gradient; \
\
    grad_tester.test(test_function, state, step_size, analytic_fisher_gradient, numeric_fisher_gradient); \
\
    if (analytic_fisher_gradient.size() != numeric_fisher_gradient.size()) \
      throw Exception ("Sizes of Fisher gradients do not match"); \
\
    for (size_t i = 0; i < analytic_fisher_gradient.size(); i++) { \
      analytic_writer.append(analytic_fisher_gradient[i]); \
      numeric_writer.append(numeric_fisher_gradient[i]); \
    } \
\
  } \
\
  reader.close(); \
  analytic_writer.close(); \
  numeric_writer.close(); \


#define ONLY_ACTIVE_FUNCTIONS

using namespace BTS; 

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Compare analytically calculated gradients against numerical approximations.",
  NULL
};

ARGUMENTS = {
  Argument ("state", "State about which the gradient will be tested.").type_file (),
  Argument ("output", "Output analysis file").type_file (),
  Argument()
};

OPTIONS = {


  Option ("object_type", "The object_type that the function belongs to.")
   + Argument ("object_type", "").type_text(),

  Option ("state_type", "The state to pass to the test_function.")
   + Argument ("state_type", "").type_text(),
  
  Option ("function_name", "The name of the function to test.")
   + Argument ("function_name", "").type_text(),
  
  Option ("hessian", "Calculate Hessian instead of gradient."),

//  Option ("fisher_info", "Calculate gradient of the Fisher information gradient."),

  Option ("rank3_hessian", "Calculate 3rd order Hessian instead of gradient."),

  Option ("step_size", "The size of the steps used to build the numerical approximation.")
   + Argument ("step_size", "").type_float (SMALL_FLOAT, 1e-4, LARGE_FLOAT),

  Option ("encoding_orientation", "Test encoding orientation.")
   + Argument ("encoding_orientation", "").type_text(),

  Option ("axis_scales_location", "Location of the file containing the scales of the Gaussian axes")
   + Argument ("axis_scales_location", "").type_text (Prob::Test::Gaussian::AXIS_SCALES_LOCATION_DEFAULT),

//  Option ("diff_encode_index", "Index of the diffusion encoding to use")
//   + Argument ("diff_encode_index", "").type_integer (0, 1000, 1),

  Option ("element_index", "Index of the 'Object' element to use.")
   + Argument("element_index", "").type_integer(0,0, LARGE_INT),

  DIFFUSION_PARAMETERS,
  
  IMAGE_PARAMETERS,

  EXPECTED_IMAGE_PARAMETERS,

  LIKELIHOOD_PARAMETERS,

  PRIOR_PARAMETERS,
  
  TEST_LANDSCAPE_PARAMETERS,

  COMMON_PARAMETERS,  

  Option ("obs_image", "The location of the reference image that is to be set")
   + Argument ("obs_image", "").type_image_in(),

  Option()
};




EXECUTE {




  //-----------------//
  //  Load Arguments //
  //-----------------//

  std::string state_location = argument[0];
  std::string output_location = argument[1];
  
  if (File::extension(state_location) != File::extension(output_location))
    throw Exception ("Extension of output, '" +  File::extension(output_location) + ", does not match extension of input, '" + File::extension(state_location) + "'.");

  std::string analytic_output_location = File::strip_extension(output_location) + str(".analytic.") + File::extension(output_location);
  std::string numeric_output_location = File::strip_extension(output_location) + str(".numeric.") + File::extension(output_location);

//
//
//
//
//  //----------------------------------//
//  //  Get and Set Optional Parameters //
//  //----------------------------------//
//
//  double step_size = 1e-4;
//  Triple<double> encoding_orientation (1.0, 0.0, 0.0);
//  std::string axis_scales_location = Prob::Test::Gaussian::AXIS_SCALES_LOCATION_DEFAULT;
//  std::string object_type  = "Prob::Prior::Magnitude";
//  std::string function_name = "log_prob";
//  bool  calculate_hessian = false;
//  bool calculate_rank3_hessian = false;
////  bool calculate_fisher_gradient = false;
////  size_t diff_encode_index = 1;
//  std::string prior_b_intens_gauss_mean_type;
//  std::string obs_image_name;
//  size_t element_index = 0;
//
//  std::string state_type;
//  if (File::has_extension<Fibre::Strand::Set>(state_location))
//    state_type = "Fibre::Strand::Set";
//
//  else if (File::has_extension<Fibre::Tractlet::Set>(state_location))
//    state_type = "Fibre::Tractlet::Set";
//
//  else if (File::has_extension<Fibre::Strand>(state_location))
//    state_type = "Fibre::Strand";
//
//  else if (File::has_extension<Fibre::Tractlet>(state_location))
//    state_type = "Fibre::Tractlet";
//
//  else if (File::has_extension< Triple<double> >(state_location))
//    state_type = "Triple";
//
//  else if (File::has_extension<Fibre::Strand::Section>(state_location))
//    state_type = "Fibre::Strand::Section";
////
////  else if (File::has_extension<Image::Expected::Buffer>(state_location))
////    state_type = "Image::Observed";
//
//  else if (File::has_extension<MCMC::State>(state_location))
//    state_type = "MCMC::State";
//
//  else
//    throw Exception ("Unrecognised extension for state, '" + File::extension(state_location) + "'.");
//
//
//
//  Options opt;
//
//  opt = get_options("object_type");
//  if (opt.size())
//    object_type = (std::string)opt[0][0];
//
//  opt = get_options("state_type");
//  if (opt.size())
//    state_type = (std::string)opt[0][0];
//
//  opt = get_options("function_name");
//  if (opt.size())
//    function_name = (std::string)opt[0][0];
//
//  opt = get_options("hessian");
//  if (opt.size())
//    calculate_hessian = true;
//
////  opt = get_options("fisher_info");
////  if (opt.size())
////    calculate_fisher_gradient = true;
//
//  opt = get_options("rank3_hessian");
//  if (opt.size())
//    calculate_rank3_hessian = true;
//
//  opt = get_options("step_size");
//  if (opt.size())
//    step_size = opt[0][0];
//
//  opt = get_options("encoding_orientation");
//  if (opt.size())
//    encoding_orientation = parse_triple<double>(opt[0][0]);
//
//  opt = get_options("axis_scales_location");
//  if (opt.size())
//    axis_scales_location = (std::string)opt[0][0];
//
////  opt = get_options("diff_encode_index");
////  if (opt.size())
////    diff_encode_index = opt[0][0];
//
//  opt = get_options("element_index");
//  if (opt.size())
//    element_index = opt[0][0];
//
//  // Loads parameters to construct Diffusion::Model ('diff_' prefix)
//  SET_DIFFUSION_PARAMETERS;
//
//  // Loads parameters to construct Image::Expected::*::Buffer that are inherited from Image::Observed::Buffer ('img_' prefix)
//  SET_IMAGE_PARAMETERS;
//
//  // Loads parameters to construct Image::Expected::*::Buffer ('img_' prefix)
//  SET_EXPECTED_IMAGE_PARAMETERS;
//
//  // Loads parameters to construct Prob::Likelihood('like_' prefix)
//  SET_LIKELIHOOD_PARAMETERS;
//
//  // Loads parameters to construct Prob::Prior('prior_' prefix)
//  SET_PRIOR_PARAMETERS;
//
//  // Loads parameters to construct Prob::Test::Landscape ('test_gen_' prefix)
//  SET_TEST_LANDSCAPE_PARAMETERS;
//
//  // Loads parameters that are common to all commands.
//  SET_COMMON_PARAMETERS;
//
//
////----------------------//
////  Load observed image //
////----------------------//
//
//  Image::Observed::Buffer obs_image;
//
//  bool obs_image_provided = false;
//
//  opt = get_options("obs_image");
//  if (opt.size()) {
//
//    obs_image_name = (std::string)opt[0][0];
//    obs_image.load(obs_image_name);
//
//    if ((obs_image.properties().count("diff_response_SH")) && (Math::matlab_str(diff_response_SH) != obs_image.properties()["diff_response_SH"]))
//      std::cout << std::endl << "Warning! Diffusion response function harmonics (" << Math::matlab_str(diff_response_SH)  << ") do not match reference image (" << obs_image.properties()["diff_response_SH"] + ")!" << std::endl;
//
//    obs_image_provided = true;
//
//  } else {
//
//    if (object_type.substr(0,18) == "Prob::Likelihood::")
//      throw Exception("Observed image location (-obs_image option) required for 'Prob::Likelihood::*' object types (" + object_type + ").");
//
//    obs_image.reset(img_dims, img_vox_lengths, img_offsets, Diffusion::Encoding::Set(diff_encodings));
//
//  }
//
//
//  //----------------------------//
//  //  Initialize Expected Image //
//  //----------------------------//
//
//
//  Diffusion::Model diffusion_model = Diffusion::Model::factory (diff_encodings,
//                                                                diff_response_SH,
//                                                                diff_adc,
//                                                                diff_fa,
//                                                                diff_isotropic);
//
//  Image::Expected::Buffer* exp_image = Image::Expected::Buffer::factory(exp_type,
//                                                       img_dims,
//                                                       img_vox_lengths,
//                                                       diffusion_model,
//                                                       exp_num_length_sections,
//                                                       exp_num_width_sections,
//                                                       exp_interp_extent,
//                                                       img_offsets,
//                                                       exp_enforce_bounds,
//                                                       exp_half_width);
//
//
//  //----------------------------------------------------------------//
//  // Auto-generate base intensity initial value/b_intens_gauss_mean //
//  //----------------------------------------------------------------//
//
//  if ((exp_base_intensity < 0) && obs_image_provided)
//    exp_base_intensity = exp_image->base_intensity_default(obs_image, state_location);
//
//  //-----------------------//
//  // Initialize Likelihood //
//  //-----------------------//
//
//  Prob::Likelihood* likelihood = 0;
//
//  if (obs_image_provided)
//    likelihood = Prob::Likelihood::factory(like_type, obs_image, exp_image, like_snr, like_b0_include, like_outside_scale, like_ref_b0, like_ref_signal);
//
//
//  //--------------------//
//  //---     Prior    ---//
//  //--------------------//
//
//  Prob::Prior prior (prior_scale,
//                    prior_mag_scale,
//                    prior_mag_aux_scale,
//                    prior_hook_scale,
//                    prior_hook_num_points,
//                    prior_density_high_scale,
//                    prior_density_low_scale,
//                    prior_density_num_points,
//                    prior_acs_scale,
//                    prior_acs_mean);
//
//
//
//  MCMC::State axis_scale (axis_scales_location);
//  Prob::Test::Gaussian test_gaussian(axis_scale);
//
//
//  Prob::Test::Landscape test_landscape;
//
//  if (lnd_location.size())
//    test_landscape = Prob::Test::Landscape(lnd_location, lnd_roi_radius, lnd_barrier_rate);
//
//
//
////-----------------//
//// Save properties //
////-----------------//
//
//  std::map<std::string,std::string> properties;
//
//  properties["method"]              = "test_gradient";
//  properties["object_type"]         = object_type;
//  properties["state_type"]          = state_type;
//  properties["state_location"]      = state_location;
//  properties["function_name"]       = function_name;
//  properties["step_size"]           = str(step_size);
//
//  if (state_type == "MCMC::State")
//    properties["axis_scales_location"] = str(axis_scales_location);
//
//  if (object_type == "Diffusion::Response")
//    properties["diff_encode_index"]   = str(obs_image_name);
//
//  if (obs_image_name.size())
//    properties["obs_image"]           = obs_image_name;
//
//  if (object_type.substr(0,9) == "Diffusion") {
//    ADD_DIFFUSION_PROPERTIES(properties);
//
//  } else if (object_type.substr(0,17) == "Prob::Likelihood") {
//    ADD_LIKELIHOOD_PROPERTIES(properties);
//
//    ADD_EXPECTED_IMAGE_PROPERTIES(properties);
//
//    ADD_DIFFUSION_PROPERTIES(properties);
//
////    ADD_DENSITY_IMAGE_PROPERTIES(properties);
//
//
//  } else if (object_type.substr(0,5) == "Image") {
//    ADD_IMAGE_PROPERTIES(properties);
//
//    if (object_type.substr(7,8) == "Expected") {
//
//      ADD_EXPECTED_IMAGE_PROPERTIES(properties);
//
//      ADD_DIFFUSION_PROPERTIES(properties);
//
//    }
////    } else if (object_type.substr(7,7) == "Density") {
////
////      ADD_DENSITY_IMAGE_PROPERTIES(properties);
////
////    }
//
//  } else if (object_type.substr(0,11) == "Prob::Prior") {
//    ADD_PRIOR_PROPERTIES(properties);
//
//  } else if (object_type == "Prob::Test::Landscape")
//    ADD_TEST_LANDSCAPE_PROPERTIES(properties);
//
//
//
//
//
//
//
////----------------------------------------------------------------------------------------------------------------//
////---------------------------------------- Select State, Object and Function -----------------------------------//
////----------------------------------------------------------------------------------------------------------------//
//
//
//
//
//
////----------------------------//
////------ Fibre::Strand -------//
////----------------------------//
//
//
//  if (state_type == "Fibre::Strand") {
//
//  //--- Prob::Prior ---//
//
//    if (object_type == "Prob::Prior") {
//
//      if (function_name == "log_prob") {
//
//        if (calculate_hessian) {
//          TEST_HESSIAN(Fibre::Strand, Prob::Prior, log_prob, prior);
//        } else {
//          TEST_GRADIENT(Fibre::Strand, Prob::Prior, log_prob, prior);
//        }
//
//      } else if (function_name == "log_prob_and_fisher") {
//
//        TEST_FISHER_GRADIENT(Fibre::Strand, Prob::Prior, log_prob_and_fisher, prior);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//
////--- Prob::Prior::Magnitude ---//
//
//    } else if (object_type == "Prob::PriorComponent::Magnitude") {
//
//      Prob::PriorComponent::Magnitude prior_curvature_magnitude (prior_mag_scale);
//
//      if (function_name == "log_prob") {
//
//        if (calculate_hessian) {
//          TEST_HESSIAN(Fibre::Strand, Prob::PriorComponent::Magnitude, log_prob, prior_curvature_magnitude);
//        } else {
//          TEST_GRADIENT(Fibre::Strand, Prob::PriorComponent::Magnitude, log_prob, prior_curvature_magnitude);
//        }
//
//      } else if (function_name == "log_prob_and_fisher") {
//
//        TEST_FISHER_GRADIENT(Fibre::Strand, Prob::PriorComponent::Magnitude, log_prob_and_fisher, prior_curvature_magnitude);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
////--- Prob::Prior::Hook ---//
//
//    } else if (object_type == "Prob::Prior::Hook") {
//
//      Prob::Prior::Hook prior_midpoint_insidecube (prior_hook_scale, prior_hook_num_points, prior_mid_cube_sd);
//
//      if (function_name == "log_prob") {
//
//        TEST_GRADIENT(Fibre::Strand, Prob::Prior::Hook, log_prob, prior_midpoint_insidecube);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//    }
////--- Image::Expected::Quartic::Buffer ---//
//
//    } else if (object_type == "Image::Expected::Quartic::Buffer") {
//
//      Image::Expected::Quartic::Buffer exp_image_quartic ( img_dims,
//                                                                  img_vox_lengths,
//                                                                  diffusion_model,
//                                                                  exp_num_length_sections,
//                                                                  exp_num_width_sections,
//                                                                  exp_interp_extent,
//                                                                  img_offsets,
//                                                                  exp_enforce_bounds);
//
//      if (function_name == "part_image") {
//
//        if (calculate_hessian) {
//          TEST_IMAGE_HESSIAN(Fibre::Strand, Image::Expected::Quartic::Buffer, part_image, exp_image_quartic);
//        } else {
//          TEST_IMAGE_GRADIENT(Fibre::Strand, Image::Expected::Quartic::Buffer, part_image, exp_image_quartic);
//        }
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//
//    } else
//      throw Exception ("Unrecognised object_type, '" + object_type + "', for state type '" + state_type + "'.");
//
//
////----------------------------//
////------ Fibre::Tractlet -------//
////----------------------------//
//
//
//  } else if (state_type == "Fibre::Tractlet") {
//
//    if (object_type == "Prob::Prior::Density") {
//
//        Prob::Prior::Density prior_density (prior_density_high_scale, prior_density_low_scale, prior_density_num_points);
//
//        if (function_name == "log_prob") {
//
//          if (calculate_hessian) {
//
//            TEST_HESSIAN(Fibre::Tractlet, Prob::Prior::Density, log_prob, prior_density);
//
//          } else {
//
//            TEST_GRADIENT(Fibre::Tractlet, Prob::Prior::Density, log_prob, prior_density);
//
//          }
//
//        } else if (function_name == "log_prob_and_fisher") {
//
//          TEST_FISHER_GRADIENT(Fibre::Tractlet, Prob::Prior::Density, log_prob_and_fisher, prior_density);
//
//        } else
//          throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//    } else if (object_type == "Image::Expected::Quartic::Buffer") {
//
//      Image::Expected::Quartic::Buffer exp_image_quartic ( img_dims,
//                                                                  img_vox_lengths,
//                                                                  diffusion_model,
//                                                                  exp_num_length_sections,
//                                                                  exp_num_width_sections,
//                                                                  exp_interp_extent,
//                                                                  img_offsets,
//                                                                  exp_enforce_bounds);
//
//      if (function_name == "part_image") {
//        if (calculate_hessian) {
//
//          Fibre::Tractlet state;
//          Fibre::Tractlet::Reader reader (state_location);
//          if (!reader.next(state))
//            throw Exception ("Did not find any states in file '" + state_location + "'.");
//          reader.close();
//
//          Fibre::Tractlet::Tensor::Writer analytic_writer (analytic_output_location + ".tnr", state), numeric_writer (numeric_output_location + ".tnr", state);
//          reader.open(state_location);
//
//          Analysis::ImageHessianTester<Image::Expected::Quartic::Buffer, Fibre::Tractlet> hess_tester (exp_image_quartic);
//          Analysis::ImageHessianTester<Image::Expected::Quartic::Buffer, Fibre::Tractlet>::Function test_function;
//
//          test_function = &Image::Expected::Quartic::Buffer::part_image;
//
//          while (reader.next(state)) {
//
//            exp_image_quartic.zero();
//
//            Image::Container::Buffer<Fibre::Tractlet::Tensor> analytic_hessian(img_dims, exp_image_quartic.num_encodings());
//            Image::Container::Buffer<Fibre::Tractlet::Tensor> numeric_hessian(img_dims, exp_image_quartic.num_encodings());
//
//            hess_tester.test(test_function, state, step_size, analytic_hessian, numeric_hessian);
//
//            for (size_t z = 0; z < exp_image_quartic.dim(Z); z++) {
//              for (size_t y = 0; y < exp_image_quartic.dim(Y); y++) {
//                for (size_t x = 0; x < exp_image_quartic.dim(X); x++) {
//                  for (size_t encode_i = 0; encode_i < exp_image_quartic.num_encodings(); encode_i++) {
//                    analytic_writer.append(analytic_hessian(x,y,z)[encode_i]);
//                    numeric_writer.append(numeric_hessian(x,y,z)[encode_i]);
//                  }
//                }
//              }
//            }
//          }
//
//          reader.close();
//          analytic_writer.close();
//          numeric_writer.close();
//
////          TEST_IMAGE_HESSIAN(Fibre::Tractlet, Image::Expected::Quartic::Buffer, part_image, exp_image_quartic);
//        } else {
//          TEST_IMAGE_GRADIENT(Fibre::Tractlet, Image::Expected::Quartic::Buffer, part_image, exp_image_quartic);
//        }
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//
//    } else
//      throw Exception ("Unrecognised object_type, '" + object_type + "', for state type '" + state_type + "'.");
//
//
////----------------------------//
////----- Fibre::Strand::Set----//
////----------------------------//
//
//  } else if (state_type == "Fibre::Strand::Set") {
//
////--- Prob::Likelihood::Gaussian ---//
//
//    if (object_type == "Prob::Likelihood::Gaussian") {
//
//      Prob::Likelihood::Gaussian likelihood_imagediff_gaussian (obs_image, exp_image, 1);
//
//      if (function_name == "log_prob") {
//        if (calculate_hessian) {
//          TEST_HESSIAN(Fibre::Strand::Set, Prob::Likelihood::Gaussian, log_prob, likelihood_imagediff_gaussian);
//        } else {
//          TEST_GRADIENT(Fibre::Strand::Set, Prob::Likelihood::Gaussian, log_prob, likelihood_imagediff_gaussian);
//        }
//      } else if (function_name == "log_prob_and_fisher") {
//
//        if (calculate_hessian) {
//
//          std::cout << "WARNING!! Requires that observed image is equal to expected image for given strand set." << std::endl;
//
//          TEST_HESSIAN(Fibre::Strand::Set, Prob::Likelihood::Gaussian, log_prob_and_fisher, likelihood_imagediff_gaussian);
//        } else {
//          TEST_FISHER_GRADIENT_MATRIX(Fibre::Strand::Set, Prob::Likelihood::Gaussian, log_prob_and_fisher, likelihood_imagediff_gaussian);
//        }
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//#ifndef ONLY_ACTIVE_FUNCTIONS
////--- Prob::Likelihood::OneSidedGaussian ---//
//    } else if (object_type == "Prob::Likelihood::OneSidedGaussian") {
//
//      Prob::Likelihood::OneSidedGaussian likelihood_weak (&obs_image, exp_image, 1);
//
//      if (function_name == "log_prob") {
//
//        TEST_GRADIENT(Fibre::Strand::Set, Prob::Likelihood::OneSidedGaussian, log_prob, likelihood_weak);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//#endif
//
////--- Prob::Likelihood::OneSidedGaussian ---//
//    } else if (object_type == "Prob::Prior") {
//
//      Prob::Prior prior_strand = *dynamic_cast<Prob::Prior*>(prior);
//
//      if (function_name == "log_prob") {
//        if (calculate_hessian) {
//          TEST_HESSIAN(Fibre::Strand::Set, Prob::Prior, log_prob, prior_strand);
//        } else {
//          TEST_GRADIENT(Fibre::Strand::Set, Prob::Prior, log_prob, prior_strand);
//        }
//
//      } else if (function_name == "log_prob_and_fisher") {
//
//        TEST_FISHER_GRADIENT_MATRIX(Fibre::Strand::Set, Prob::Prior, log_prob_and_fisher, prior_strand);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//    } else
//      throw Exception ("Unrecognised object_type, '" + object_type + "', for state type '" + state_type + "'.");
//
//
//
////----------------------------//
////----- Fibre::Tractlet::Set ----//
////----------------------------//
//
//  } else if (state_type == "Fibre::Tractlet::Set") {
//
//
////--- Prob::Likelihood::Gaussian ---//
//    if (object_type == "Prob::Likelihood::Gaussian") {
//
//      Prob::Likelihood::Gaussian likelihood_imagediff_gaussian (obs_image, exp_image, 1);
//
//      if (function_name == "log_prob") {
//
//        if (calculate_hessian) {
//          TEST_HESSIAN(Fibre::Tractlet::Set, Prob::Likelihood::Gaussian, log_prob, likelihood_imagediff_gaussian);
//        } else {
//          TEST_GRADIENT(Fibre::Tractlet::Set, Prob::Likelihood::Gaussian, log_prob, likelihood_imagediff_gaussian);
//        }
//
//      } else if (function_name == "log_prob_and_fisher") {
//
//        TEST_FISHER_GRADIENT_MATRIX(Fibre::Tractlet::Set, Prob::Likelihood::Gaussian, log_prob_and_fisher, likelihood_imagediff_gaussian);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//#ifndef ONLY_ACTIVE_FUNCTIONS
////--- Prob::Likelihood::OneSidedGaussian ---//
//
//    } else if (object_type == "Prob::Likelihood::OneSidedGaussian") {
//
//      Prob::Likelihood::OneSidedGaussian likelihood_weak (&obs_image, exp_image, 1);
//
////      if (function_name == "log_prob") {
////
////        TEST_GRADIENT(Fibre::Tractlet::Set, Prob::Likelihood::OneSidedGaussian, log_prob, likelihood_weak);
////
////      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//#endif
//
//    } else
//      throw Exception ("Unrecognised object_type, '" + object_type + "', for state type '" + state_type + "'.");
//
//
//
////-----------------//
////----- Triple ----//
////-----------------//
//
//  } else if (state_type == "Triple") {
//
//    if (object_type == "Diffusion::Response") {
//
//      if (function_name == "weighting")  {
//
//        if (calculate_hessian) {
//          TEST_OBJECT_ELEMENT_HESSIAN(Coord, Diffusion::Response, weighting, diffusion_model);
//        } else {
//          TEST_OBJECT_ELEMENT_GRADIENT(Coord, Diffusion::Response, weighting, diffusion_model);
//        }
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//#ifndef ONLY_ACTIVE_FUNCTIONS
//    } else if (object_type == "Image::Expected::Trilinear::Voxel") {
//
//
//      Image::Expected::Trilinear::Buffer exp_image_trilinear ( img_dims,
//                                                                  img_vox_lengths,
//                                                                  diffusion_model,
//                                                                  exp_num_length_sections,
//                                                                  exp_num_width_sections,
//                                                                  exp_interp_extent,
//                                                                  img_offsets,
//                                                                  exp_enforce_bounds);
//
//      if (function_name == "interpolate") {
//
//        TEST_VOXEL_GRADIENT(Triple<double>, Image::Expected::Trilinear::Voxel, interpolate, exp_image_trilinear);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//#endif
//#ifndef ONLY_ACTIVE_FUNCTIONS
//    } else if (object_type == "Image::Expected::Gaussian::Voxel") {
//
//
//      Image::Expected::Gaussian::Buffer exp_image_gaussian ( img_dims,
//                                                                img_vox_lengths,
//                                                                diffusion_model,
//                                                                exp_num_length_sections,
//                                                                exp_num_width_sections,
//                                                                exp_interp_extent,
//                                                                exp_half_width,
//                                                                img_offsets,
//                                                                exp_enforce_bounds);
//
//      if (function_name == "interpolate") {
//
//        TEST_VOXEL_GRADIENT(Triple<double>, Image::Expected::Gaussian::Voxel, interpolate, exp_image_gaussian);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//#endif
//    } else if (object_type == "Image::Expected::Quartic::Voxel") {
//
//
//      Image::Expected::Quartic::Buffer exp_image_quartic ( img_dims,
//                                                              img_vox_lengths,
//                                                              diffusion_model,
//                                                              exp_num_length_sections,
//                                                              exp_num_width_sections,
//                                                              exp_interp_extent,
//                                                              img_offsets,
//                                                              exp_enforce_bounds);
//
//      if (function_name == "interpolate") {
//
//        if (calculate_hessian) {
//          TEST_VOXEL_HESSIAN(Coord, Image::Expected::Quartic::Voxel, interpolate, exp_image_quartic);
//        } else {
//          TEST_VOXEL_GRADIENT(Coord, Image::Expected::Quartic::Voxel, interpolate, exp_image_quartic);
//        }
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//    } else
//          throw Exception ("Unrecognised -object_type '" + object_type + "' for state type '" + state_type + "'.");
//
//
////---------------------------------------//
////----- Fibre::Strand::Section ----//
////---------------------------------------//
//
//  } else if (state_type == "Fibre::Strand::Section") {
//
//#ifdef CODE_TESTING
//
//    if (object_type == "Image::Expected::Quartic::Voxel") {
//
//
//        Image::Expected::Quartic::Buffer exp_image_quartic ( img_dims,
//                                                                img_vox_lengths,
//                                                                diffusion_model,
//                                                                exp_num_length_sections,
//                                                                exp_num_width_sections,
//                                                                exp_interp_extent,
//                                                                img_offsets,
//                                                                exp_enforce_bounds);
//
//        if (function_name == "interpolate") {
//
//          if (calculate_hessian) {
//            TEST_VOXEL_HESSIAN(Fibre::Strand::Section, Image::Expected::Quartic::Voxel, interpolate, exp_image_quartic);
//          } else {
//            TEST_VOXEL_GRADIENT(Fibre::Strand::Section, Image::Expected::Quartic::Voxel, interpolate, exp_image_quartic);
//          }
//
//        } else
//          throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//
//#ifndef ONLY_ACTIVE_FUNCTIONS
//    } else if (object_type == "Image::Expected::Trilinear::Direction") {
//
//
//      Image::Expected::Trilinear::Buffer exp_image_trilinear ( img_dims,
//                                                                  img_vox_lengths,
//                                                                  diffusion_model,
//                                                                  exp_num_length_sections,
//                                                                  exp_num_width_sections,
//                                                                  exp_interp_extent,
//                                                                  img_offsets,
//                                                                  exp_enforce_bounds);
//
//      if (function_name == "signal") {
//
//        TEST_DIRECTION_GRADIENT(Fibre::Strand::Section, signal, exp_image_trilinear);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//#endif
//#ifndef ONLY_ACTIVE_FUNCTIONS
//    } else if (object_type == "Image::Expected::Gaussian::Direction") {
//
//
//      Image::Expected::Gaussian::Buffer exp_image_gaussian ( img_dims,
//                                                                img_vox_lengths,
//                                                                diffusion_model,
//                                                                exp_num_length_sections,
//                                                                exp_num_width_sections,
//                                                                exp_interp_extent,
//                                                                exp_half_width,
//                                                                img_offsets,
//                                                                exp_enforce_bounds);
//
//       if (function_name == "signal") {
//
//        TEST_DIRECTION_GRADIENT(Fibre::Strand::Section, signal, exp_image_gaussian);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//#endif
//    } else if (object_type == "Image::Expected::Quartic::Direction") {
//
//
//      Image::Expected::Quartic::Buffer exp_image_quartic ( img_dims,
//                                                              img_vox_lengths,
//                                                              diffusion_model,
//                                                              exp_num_length_sections,
//                                                              exp_num_width_sections,
//                                                              exp_interp_extent,
//                                                              img_offsets,
//                                                              exp_enforce_bounds);
//
//      if (function_name == "signal") {
//        if (calculate_hessian) {
//          TEST_DIRECTION_HESSIAN(Fibre::Strand::Section, signal, exp_image_quartic);
//        } else {
//          TEST_DIRECTION_GRADIENT(Fibre::Strand::Section, signal, exp_image_quartic);
//        }
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//
//
//    } else
//      throw Exception ("Unrecognised -object_type '" + object_type + "' for state type '" + state_type + "'.");
//
//#else
//
//  throw Exception ("test_gradient must be compiled with CODE_TESTING macro to test Fibre::Strand::Sections.");
//
//#endif /* CODE_TESTING */
//
//  } else if (state_type == "Fibre::Tractlet::Section") {
//
//#ifdef CODE_TESTING
//
//    if (object_type == "Image::Expected::Direction::Quartic") {
//
//      Image::Expected::Quartic::Buffer exp_image_quartic ( img_dims,
//                                                              img_vox_lengths,
//                                                              diffusion_model,
//                                                              exp_num_length_sections,
//                                                              exp_num_width_sections,
//                                                              exp_interp_extent,
//                                                              img_offsets,
//                                                              exp_enforce_bounds);
//
//      if (function_name == "signal") {
//
//        TEST_DIRECTION_GRADIENT(Fibre::Tractlet::Section, signal, exp_image_quartic);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//#ifndef ONLY_ACTIVE_FUNCTIONS
//    } else if (object_type == "Image::Expected::Direction::Trilinear") {
//
//
//      Image::Expected::Trilinear::Buffer exp_image_trilinear ( img_dims,
//                                                                  img_vox_lengths,
//                                                                  diffusion_model,
//                                                                  exp_num_length_sections,
//                                                                  exp_num_width_sections,
//                                                                  exp_interp_extent,
//                                                                  img_offsets,
//                                                                  exp_enforce_bounds);
//
//      if (function_name == "signal") {
//
//        TEST_DIRECTION_GRADIENT(Fibre::Tractlet::Section, signal, exp_image_trilinear);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//#endif
//#ifndef ONLY_ACTIVE_FUNCTIONS
//    } else if (object_type == "Image::Expected::Direction::Gaussian") {
//
//
//      Image::Expected::Gaussian::Buffer exp_image_gaussian ( img_dims,
//                                                                img_vox_lengths,
//                                                                diffusion_model,
//                                                                exp_num_length_sections,
//                                                                exp_num_width_sections,
//                                                                exp_interp_extent,
//                                                                exp_half_width,
//                                                                img_offsets,
//                                                                exp_enforce_bounds);
//
//       if (function_name == "signal") {
//
//        TEST_DIRECTION_GRADIENT(Fibre::Tractlet::Section, signal, exp_image_gaussian);
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//#endif
//    } else
//      throw Exception ("Unrecognised -object_type '" + object_type + "' for state type '" + state_type + "'.");
//
//
//
//#else
//
//    throw Exception ("test_gradient must be compiled with CODE_TESTING macro to test Fibre::Tractlet::Sections.");
//
//#endif /* CODE_TESTING */
//
//  //---------------------------//
//  //----- Image::Observed ----//
//  //---------------------------//
//
//  } else if (state_type == "Image::Observed") {
//
//
//    Image::Observed::Buffer image (state_location); //, Diffusion::Encoding::Set(diff_encodings));
//
//
////    if (object_type == "Image::Observed") {
////
////      if (function_name == "sum_strong_difference") {
////
////        TEST_IMAGE_POINTER_GRADIENT(Image::Observed, Image::Observed, sum_strong_difference, obs_image, image);
////
////      } else if (function_name == "sum_weak_difference") {
////
////        TEST_IMAGE_POINTER_GRADIENT(Image::Observed, Image::Observed, sum_weak_difference, obs_image, image);
////
////      } else
////        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
////    } else
////      throw Exception ("Unrecognised -object_type '" + object_type + "' for state type '" + state_type + "'.");
////
////
//    throw Exception ("deprecated.");
//
////----------------------//
////----- MCMC::State ----//
////----------------------//
//
//  } else if (state_type == "MCMC::State") {
//
//    if (object_type == "Prob::Test::Gaussian") {
//
//      if (function_name == "log_prob") {
//
//        if (calculate_hessian) {
//          TEST_HESSIAN(MCMC::State, Prob::Test::Gaussian, log_prob, test_gaussian);
//        } else {
//          TEST_GRADIENT(MCMC::State, Prob::Test::Gaussian, log_prob, test_gaussian);
//        }
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//    } else if (object_type == "Prob::Test::Landscape") {
//
//      if (function_name == "log_prob") {
//
//        if (calculate_rank3_hessian) {
//          TEST_RANK3_HESSIAN(MCMC::State, Prob::Test::Landscape, log_prob, test_landscape);
//        } else if (calculate_hessian) {
//          TEST_HESSIAN(MCMC::State, Prob::Test::Landscape, log_prob, test_landscape);
//        } else {
//          TEST_GRADIENT(MCMC::State, Prob::Test::Landscape, log_prob, test_landscape);
//        }
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//
//    }  else if (object_type == "Prob::Test::Landscape::Peak") {
//
//      if (function_name == "gaussian_log_prob") {
//
//
//        if (calculate_rank3_hessian) {
//          TEST_RANK3_HESSIAN(MCMC::State, Prob::Test::Landscape::Peak, gaussian_log_prob, test_landscape[element_index]);
//        } else if (calculate_hessian) {
//          TEST_HESSIAN(MCMC::State, Prob::Test::Landscape::Peak, gaussian_log_prob, test_landscape[element_index]);
//        } else {
//          TEST_GRADIENT(MCMC::State, Prob::Test::Landscape::Peak, gaussian_log_prob, test_landscape[element_index]);
//        }
//
//      } else
//        throw Exception ("Unrecognised function '" + object_type + "::" + function_name + "(const " + state_type + "&, " + state_type + "&)'.");
//
//
//    } else
//      throw Exception ("Unrecognised -object_type '" + object_type + "' for state type '" + state_type + "'.");
//
//
//  } else
//    throw Exception ("Unrecognised state_type, '" + state_type + "'.");
//
//
//  delete likelihood;
//  delete exp_image;
//
//  std::cout << std::endl << "Generated analytic and numerical gradients: " << std::endl;
//  std::cout << "  State: " << state_type << std::endl;
//  std::cout << "  Object: " << object_type << std::endl;
//  std::cout << "  Function: " + function_name << "\n" << std::endl;
//
//  std::cout << "\nMATLAB Plot command:\n\nclose all; ";
//
//  if (function_name == "log_prob_and_fisher")
//    std::cout << "plot_fisher_gradient " + output_location << ".tnr\n\n" << std::endl;
//  else if (calculate_rank3_hessian)
//    std::cout << "plot_rank3_hessian " + output_location << ".tnr3\n\n" << std::endl;
//  else if (calculate_hessian)
//    std::cout << "plot_hessian " + output_location << ".tnr\n\n" << std::endl;
//  else
//    std::cout << "plot_gradient " + output_location << "\n\n" << std::endl;
////

}


