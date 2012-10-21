/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Aug 2, 2010.

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

#include "bts/mcmc/proposal/momentum.h"
#include "bts/mcmc/proposal/walker.h"


#ifndef __bts_mcmc_state_h__
#define __bts_mcmc_state_h__

#include <map>

namespace BTS {

  namespace MCMC {

    class State;

  }
}


#include "math/vector.h"
#include "bts/utilities/reader.h"
#include "bts/utilities/writer.h"


namespace BTS {

	namespace MCMC {

		class State : public MR::Math::Vector<double> {

			//Public static variables, nested classes and typedefs
			public:

        typedef Proposal::Walker Walker;
        typedef Proposal::Momentum Momentum;
        typedef std::map<std::string, std::string> Properties;
        typedef Utilities::Reader<MCMC::State> Reader;
        typedef Utilities::Writer<MCMC::State> Writer;

        template <typename T> class Tensor_tpl;
        class Tensor;

			public:

        Properties properties;


			public:

        static const std::string FILE_EXTENSION;

			public:


        //Do nothing, no characteristic properties in MCMC::State class.
        static std::vector<std::string>&  append_characteristic_keys(std::vector<std::string>& header) { return header; }

			//Public member functions
			public:

        //! construct empty state
        State () throw () {}

        //! copy constructor
        State (const State& state)
          : MR::Math::Vector<double> (state), properties(state.properties) {}

        State(const MR::Math::Vector<double>& vector)
          : MR::Math::Vector<double>(vector) {}

        //! construct state of size \a nelements
        /** \note the elements of the state are left uninitialised. */
        State (size_t nelements)
        : MR::Math::Vector<double>(nelements) {}

        //! construct from existing data array
        State (double* state_data, size_t nelements, size_t skip = 1) throw ()
          : MR::Math::Vector<double>(state_data,nelements,skip) {}

        //! construct a state by reading from the text file \a filename
        State (const std::string& file)
          : MR::Math::Vector<double>(file) {}


        State&  operator= (const State& state)
          { MR::Math::Vector<double>::operator=(state); properties = state.properties; return *this; }


        //! destructor
        ~State () {}



        double                         norm2() const
          { return MR::Math::norm2(*this); }

        double                         norm() const
          { return MR::Math::norm(*this); }

        double                         inner_product(const State& state) const
          { return MR::Math::dot(*this, state); }

        State&                         zero()
          { for (size_t i = 0; i < size(); i++) operator[](i) = 0.0; return *this;}

        State&                         invalidate()
          { for (size_t i = 0; i < size(); i++) operator[](i) = NAN; return *this; }

        State                          operator- ()    { operator*=(-1.0); return *this; }


        void                           set_extend_prop(const std::string& prop, const std::string& value) {

          properties[prop] = value;

        }




        State                          operator+(const State& state) const { State answer = *this; answer += state; return answer; }
        State                          operator-(const State& state) const { State answer = *this; answer -= state; return answer; }

        State                          operator* (const State& state) const { State answer = *this; answer *= state; return answer; }
        State                          operator/ (const State& state) const { State answer = *this; answer /= state; return answer; }

        State                          operator* (double scalar) const { State answer = *this; answer *= scalar; return answer; }
        State                          operator/ (double scalar) const { State answer = *this; answer /= scalar; return answer; }


        MR::Math::Vector<double>      vector() const
          { return *this; }

        void                          from_vector(const MR::Math::Vector<double>& vec)
          { MR::Math::Vector<double>::operator=(vec); }

        size_t                          vsize() const
          { return size(); }


        void                          set_characteristics() {}

		};


		inline State                      operator*(double scalar, const State& state)
      { return state * scalar; }

	}

}

#endif /* __bts_mcmc_state_h__ */
