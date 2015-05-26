/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 01/10/2010.

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

#include "bts/prob/test/bayes_log_regression.h"
#include "bts/mcmc/state/tensor.h"

namespace FTS {
    
    namespace Prob {
        
        namespace Test {
            
            const char* BayesLogRegression::DATA_LOCATION_DEFAULT =
                    "/home/tclose/data/prob/test/bayes_log_regression/australian.txt";
            const size_t BayesLogRegression::POLY_ORDER_DEFAULT = 1;
            const float BayesLogRegression::PRIOR_VARIANCE_DEFAULT = 100.0;
            
            BayesLogRegression::BayesLogRegression(const std::string& location, size_t poly_order,
                                                   double alpha)
                    : poly_order(poly_order), alpha(alpha) {
                
                MR::Math::Matrix<double> loaded_data(location);
                MR::Math::Matrix<double> X = loaded_data.sub(0, loaded_data.rows(), 0,
                        loaded_data.columns() - 1);
                t = loaded_data.column(loaded_data.columns() - 1);
                
                XX.resize(X.rows(), X.columns() * poly_order + 1);
                XX.column(0) = 1.0;
                
                for (size_t poly_i = 0; poly_i < poly_order; ++poly_i)
                    for (size_t col_i = 0; col_i < X.columns(); ++col_i)
                        for (size_t row_i = 0; row_i < X.rows(); ++row_i)
                            XX(row_i, poly_i * X.columns() + col_i + 1) = MR::Math::pow(
                                    X(row_i, col_i), poly_i + 1);
                
                N = XX.rows();
                D = XX.columns();
                
                f.resize(N);
                working.resize(D, N);
                v.resize(N);
                p.resize(N);
                V.resize(D, N);
                
                Z2.resize(N, D);
                
            }
            
            double BayesLogRegression::log_prob(const MR::Math::Vector<double>& w) {
                
                double lprior = 0.0;
                
                for (size_t d = 0; d < D; ++d)
                    lprior -= 0.5 * MR::Math::log(2.0 * M_PI * alpha)
                            - MR::Math::pow2(w[d]) / (2.0 * alpha);
                
                MR::Math::mult(f, XX, w);
                
                double llike = MR::Math::dot(f, t);
                
                for (size_t n = 0; n < N; ++n)
                    llike -= MR::Math::log(1.0 + MR::Math::exp(f[n]));
                
                return lprior + llike;
                
            }
            
            double BayesLogRegression::log_prob(const MR::Math::Vector<double>& w,
                                                MR::Math::Vector<double>& d_w) {
                
                MR::Math::mult(f, XX, w);
                
                for (size_t n = 0; n < N; ++n)
                    v[n] = t[n] - 1.0 / (1.0 + MR::Math::exp(-f[n]));
                
                MR::Math::mult(d_w, 0.0, 1.0, CblasTrans, XX, v);
                
                for (size_t d = 0; d < D; ++d)
                    d_w[d] -= 1.0 / alpha * w[d];
                
                return log_prob(w);
                
            }
            
            double BayesLogRegression::log_prob_and_fisher(const MR::Math::Vector<double>& w,
                                                           MR::Math::Vector<double>& d_w,
                                                           MR::Math::Matrix<double>& G) {
                
                assert(w.size() == D);
                
                G.resize(D, D);
                d_w.resize(D);
                
                MR::Math::mult(f, XX, w);
                
                for (size_t n = 0; n < N; ++n) {
                    p[n] = 1 / (1 + MR::Math::exp(-f[n]));
                    v[n] = p[n] * (1 - p[n]);
                    for (size_t d = 0; d < D; ++d)
                        working(d, n) = v[n] * XX(n, d);
                }
                
                MR::Math::mult(G, working, XX);
                
                for (size_t d = 0; d < D; ++d)
                    G(d, d) += 1.0 / alpha;
                
                return log_prob(w, d_w);
                
            }
            
            double BayesLogRegression::log_prob_and_fisher(const MR::Math::Vector<double>& w,
                                                           MR::Math::Vector<double>& d_w,
                                                           MR::Math::Matrix<double>& G,
                                                           std::vector<MCMC::State::Tensor>& d_G) {
                
                d_G.clear();
                d_G.resize(D, MCMC::State::Tensor(w));
                
                Z = 0.0;
                
                MR::Math::mult(f, XX, w);
                
                for (size_t n = 0; n < N; ++n) {
                    p[n] = 1 / (1 + MR::Math::exp(-f[n]));
                    v[n] = p[n] * (1 - p[n]);
                    for (size_t d = 0; d < D; ++d)
                        V(d, n) = v[n];
                }
                
                for (size_t d = 0; d < D; ++d) {
                    
                    for (size_t n = 0; n < N; ++n) {
                        double z1 = (1.0 - 2.0 * p[n]) * XX(n, d) * v[n];
                        for (size_t d2 = 0; d2 < D; ++d2)
                            Z2(n, d2) = XX(n, d2) * z1;
                    }
                    
                    MR::Math::Matrix<double> d_G_d;
                    
                    MR::Math::mult(d_G_d, 1.0, CblasTrans, XX, CblasNoTrans, Z2);
                    
                    d_G[d] = d_G_d;
                    
                }
                
                return log_prob_and_fisher(w, d_w, G);
                
            }
        
        }
    
    }

}
