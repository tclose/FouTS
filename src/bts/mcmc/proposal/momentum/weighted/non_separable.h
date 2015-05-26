/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Sep 27, 2010.

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

#ifndef __bts_mcmc_proposal_momentum_weighted_nonseparable_h__
#define __bts_mcmc_proposal_momentum_weighted_nonseparable_h__

#include "bts/math/common.h"
#include "bts/mcmc/naninf_exception.h"
#include "bts/mcmc/proposal/momentum/weighted.h"

namespace FTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            class Momentum::Weighted::NonSeparable: public Momentum::Weighted {
                    
                    //Public static variables, nested classes and typedefs
                public:
                    
                    //Protected member variables
                protected:
                    
                    //Working matrices and vectors.
                    MR::Math::Vector<double> tmp_momen;
                    MR::Math::Vector<double> working1;
                    MR::Math::Vector<double> working2;
                    std::vector<MR::Math::Matrix<double> > Winv_dWdx;
                    std::vector<double> trace_Winv_dWdx;

                    //Public static functions
                public:
                    
                    template<typename T> static NonSeparable factory(
                            const T& state, double step_scale, const std::string& step_location,
                            MCMC::Proposal::Distribution* const proposal_distribution,
                            size_t num_newton_steps);

                    //Public member functions
                public:
                    
                    NonSeparable(Distribution* const prop_distr = 0)
                            : Weighted(prop_distr) {
                    }
                    
                    NonSeparable(Distribution* const proposal_distribution,
                                 const MR::Math::Vector<double>& step_sizes,
                                 size_t num_newton_steps)
                            : Weighted(proposal_distribution, step_sizes), tmp_momen(
                                      step_sizes.size()), working1(step_sizes.size()), working2(
                                      step_sizes.size()), Winv_dWdx(step_sizes.size(),
                                      MR::Math::Matrix<double>(step_sizes.size(),
                                              step_sizes.size())), trace_Winv_dWdx(
                                      step_sizes.size()) {
                    }
                    
                    virtual ~NonSeparable() {
                    }
                    
                    NonSeparable(const NonSeparable& NS)
                            : Weighted(NS), tmp_momen(NS.tmp_momen), working1(NS.working1), working2(
                                      NS.working2), Winv_dWdx(NS.Winv_dWdx), trace_Winv_dWdx(
                                      NS.trace_Winv_dWdx) {
                    }
                    
                    NonSeparable& operator=(const NonSeparable& NS) {
                        
                        Weighted::operator=(NS);
                        tmp_momen = NS.tmp_momen;
                        working1 = NS.working1;
                        working2 = NS.working2;
                        Winv_dWdx = NS.Winv_dWdx;
                        trace_Winv_dWdx = NS.trace_Winv_dWdx;
                        
                        return *this;
                    }
                    
                    void randomize(const MR::Math::Matrix<double>& weights_chol) {
                        Momentum::Weighted::randomize(weights_chol);
                    }
                    
//          double                                  log_kinetic_energy(const MR::Math::Matrix<double>& weights_chol) const;
                    
                    template<typename T> void half_update_momentum(
                            const MR::Math::Vector<double>& gradient,
                            const MR::Math::Matrix<double>& weights,
                            std::vector<typename T::Tensor>& weights_gradient,
                            double time_direction, size_t num_newton_steps = 1);

                    double predicted_change(const MR::Math::Vector<double>& gradient,
                                            const MR::Math::Matrix<double>& fisher_chol,
                                            double time_direction) const;

                    template<typename T> double predicted_change(
                            const T& gradient, const MR::Math::Matrix<double>& fisher_chol,
                            double time_direction) const {
                        return predicted_change(gradient, fisher_chol, time_direction);
                    }
                    
            };
            
            template<typename T> void Momentum::Weighted::NonSeparable::half_update_momentum(
                    const MR::Math::Vector<double>& gradient,
                    const MR::Math::Matrix<double>& weights_chol,
                    std::vector<typename T::Tensor>& weights_gradient, double time_direction,
                    size_t num_newton_steps) {
                
                assert((time_direction == -1.0) || (time_direction == 1.0));
                
                for (size_t elem_i = 0; elem_i < size(); ++elem_i) {
                    
                    //More numerically sound than calculating the inverse and multiplying it with the gradient.
                    for (size_t col_i = 0; col_i < size(); ++col_i) {
                        MR::Math::Vector<double> Winv_dWdx_col = Winv_dWdx[elem_i].column(col_i);
                        MR::Math::Cholesky::solve(Winv_dWdx_col, weights_chol,
                                weights_gradient[elem_i].column(col_i));
                    }
                    
                    trace_Winv_dWdx[elem_i] = Math::trace(Winv_dWdx[elem_i]) / 2.0;
                }
                
                tmp_momen = momen;
                
                MR::Math::Vector<double> weights_momen(size());
                
                for (size_t newton_i = 0; newton_i < num_newton_steps; ++newton_i) {
                    
                    MR::Math::Cholesky::solve(weights_momen, weights_chol, tmp_momen);
                    
                    for (size_t elem_i = 0; elem_i < size(); ++elem_i) {
                        MR::Math::mult(working1, Winv_dWdx[elem_i], weights_momen);
                        working2[elem_i] = MR::Math::dot(working1, tmp_momen) / 2.0;
                        
                    }
                    
                    for (size_t elem_i = 0; elem_i < size(); ++elem_i) {
                        
                        tmp_momen[elem_i] = momen[elem_i]
                                + 0.5 * time_direction * step[elem_i]
                                  * (gradient[elem_i] - trace_Winv_dWdx[elem_i] + working2[elem_i]);
                        
                        if (isnan(tmp_momen[elem_i]) || isinf(tmp_momen[elem_i]))
                            throw NanInfException();
//              throw Exception ("Momentum has NaN or inf values in " + str(newton_i) + " newton step.\n\n" + str(tmp_momen) + "\n\n");
                    }
                    
                }
                
                momen = tmp_momen;
                
            }
        
        }
    
    }

}

#endif /* __bts_mcmc_proposal_momentum_weighted_non_separable_h__ */
