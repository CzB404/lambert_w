
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_LAMBERT_W_RAW_HPP_INCLUDED
#define BOOST_MATH_LAMBERT_W_RAW_HPP_INCLUDED

#include <cmath>
#include <complex>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/policies/policy.hpp>
#include "lambert_w_tools.hpp"
#include "lambert_w_tay.hpp"
#include "lambert_w_laur.hpp"
#include "lambert_w_sing.hpp"
#include "lambert_w_asy.hpp"
#include "lambert_w_mid.hpp"
#include "lambert_w_iter.hpp"

#define BOOST_MATH_LAMBERT_W_NEWTON_ITERATION
#define BOOST_MATH_LAMBERT_W_FRITSCH_ITERATION

namespace boost
{
namespace math
{
namespace lambw
{
     //Unfinished dynamic approximation selector
     /*
     template<class ArgumentType, std::size_t L>
     ArgumentType _best_approx(const ArgumentType &z, const boost::array<ArgumentType,L> &approxes)
     {
          int prec;
          int prec_max = _test_score(z,approxes[0]);
          std::size_t max_place = 0;

          for(std::size_t k = 1; k < approxes.size(); ++k)
          {
               prec = _test_score(z,approxes[k]);
               if(prec_max < prec)
               {
                    prec_max = prec;
                    max_place = k;
               }
          }
          return approxes[max_place];
     }
     */

     //This function determines which approximation to use as input for iterations, then applies those iterations.
     //It makes a few Newton's method iterations to refine the input for Fritsch's iteration, which runs until it hits machine precision.
     //Policies are accepted as arguments, but because this function template accepts both real and complex values,
     //policy-based errors have not been introduced in here. Instead, there are separate checking function following this one.
     template<class ArgumentType, class CoeffType, class IndexType, class Policy>
     ArgumentType _raw(const ArgumentType &z, IndexType k, const Policy& pol)
     {
          using std::abs;
          using std::exp;

          ArgumentType w;

          CoeffType z_imag = _complex_imag(z);

          if(k==-1)
          {
               if(abs(z+(CoeffType)0.04)<=0.14 && z_imag>=0)
               {
                    w = _N1<ArgumentType,CoeffType>(z);
               }
               else
               {
                    if(abs(z+rec_e<CoeffType>())<0.3 && z_imag>=0)
                    {
                         w = _sing<ArgumentType,CoeffType,IndexType>(z,k);
                    }
                    else
                    {
                         w = _asy<ArgumentType,CoeffType,IndexType>(z,-1);
                    }
               }
          }
          else
          {
               if(k==1)
               {
                    if(abs(z+(CoeffType)0.04)<=0.14 && z_imag<0)
                    {
                         w = _N1<ArgumentType,CoeffType>(z);
                    }
                    else
                    {
                         if(abs(z+rec_e<CoeffType>())<0.3 && z_imag<0)
                         {
                              w = _sing<ArgumentType,CoeffType,IndexType>(z,k);
                         }
                         else
                         {
                              w = _asy<ArgumentType,CoeffType,IndexType>(z,1);
                         }
                    }
               }
               else
               {
                    if(k!=0)
                    {
                         w = _asy<ArgumentType,CoeffType,IndexType>(z,k);
                    }
                    else
                    {
                         if(abs(z-(CoeffType)0.02)<0.28)
                         {
                              w = _laur<ArgumentType,CoeffType>(z);
                         }
                         else
                         {
                              if(abs(z+rec_e<CoeffType>())<0.3)
                              {
                                   w = _sing<ArgumentType,CoeffType,IndexType>(z,k);
                              }
                              else
                              {
                                   if(abs(z-(CoeffType)8.)<13)
                                   {
                                        w = _mid<ArgumentType,CoeffType>(z);
                                   }
                                   else
                                   {
                                        w = _asy<ArgumentType,CoeffType,IndexType>(z,0);
                                   }
                              }
                         }
                    }
               }
          }

          ///Unfinished dynamic approximation selector
          /*
          if(k == 0)
          {
               boost::array<ArgumentType,4> approxes = {_laur<ArgumentType,CoeffType>(z),
                                                            _sing<ArgumentType,CoeffType,IndexType>(z,0),
                                                            _mid<ArgumentType,CoeffType>(z),
                                                            _asy<ArgumentType,CoeffType,IndexType>(z,0)};
               w = _best_approx(z,approxes);
          }
          else if(k == -1)
          {
               if(_complex_imag(z) >= 0)
               {
                    boost::array<ArgumentType,3> approxes = {_sing<ArgumentType,CoeffType,IndexType>(z,-1),
                                                                 _N1<ArgumentType,CoeffType>(z),
                                                                 _asy<ArgumentType,CoeffType,IndexType>(z,-1)};
                    w = _best_approx(z,approxes);
               }
               else
               {
                    w = _asy<ArgumentType,CoeffType,IndexType>(z,-1);
               }
          }
          else if(k == 1)
          {
               if(_complex_imag(z) < 0)
               {
                    boost::array<ArgumentType,3> approxes = {_sing<ArgumentType,CoeffType,IndexType>(z,1),
                                                                 _N1<ArgumentType,CoeffType>(z),
                                                                 _asy<ArgumentType,CoeffType,IndexType>(z,1)};
                    w = _best_approx(z,approxes);
               }
               else
               {
                    w = _asy<ArgumentType,CoeffType,IndexType>(z,1);
               }
          }
          else
          {
               w = _asy<ArgumentType,CoeffType,IndexType>(z,k);
          }
          */

          if( (boost::math::isnan(_complex_real(w))
               ||boost::math::isnan(_complex_imag(w)))
               || (_complex_real(z)==0.
               && _complex_imag(z)==0.) )
          {
               return w;
          }

          CoeffType err = abs(_test_s<ArgumentType>(z,w,pol));

          #ifdef BOOST_MATH_LAMBERT_W_NEWTON_ITERATION
          for(std::size_t cntr = 0; cntr < 3; ++cntr)
          {
               _iter_newt<ArgumentType,CoeffType,IndexType>(z,k,w);
          }
          #endif

          ArgumentType w_prev = w;

          CoeffType err_prev;

          #ifdef BOOST_MATH_LAMBERT_W_FRITSCH_ITERATION
          do
          {
               err_prev = err;
               w_prev = w;
               _iter_frit<ArgumentType,CoeffType,IndexType>(z,k,w);
               err = abs(_test_s<ArgumentType>(z,w,pol));
          }
          while(err < err_prev);

          _iter_newt<ArgumentType,CoeffType,IndexType>(z,k,w_prev);
          #endif

          return w_prev;
     }

     //Checks for errors in input arguments and the return value using the policy pol for real valued arguments.
     template<class ArgumentType, class IndexType, class Policy>
     inline ArgumentType _with_real_range_checks(const ArgumentType& z, const IndexType& k, const Policy& pol)
     {
          ArgumentType w = (ArgumentType)0.;

          if(boost::math::isnan(z))
          {
               w = boost::math::policies::raise_evaluation_error<ArgumentType>("lambert_w(%1%,IndexType)", "Argument value %1% is Not-a-Number!",z,pol);
          }
          else if( ( k == 0 || k == -1 ) )
          {
               if( z < -lambw::rec_e<ArgumentType>())
               {
                    w = boost::math::policies::raise_domain_error<ArgumentType>("lambert_w(%1%,IndexType)", "Argument value %1% out of range for real arguments.", z, pol);
               }
               else if( (k == 0) && (z == std::numeric_limits<ArgumentType>::infinity()) )
               {
                    w = boost::math::policies::raise_overflow_error<ArgumentType>("lambert_w(%1%,IndexType)", "Result value overflow for an infinite argument at 0 index.", pol);
               }
               else if( k == -1 )
               {
                    if(z == -0.)
                    {
                         w = -boost::math::policies::raise_overflow_error<ArgumentType>("lambert_w(%1%,IndexType)", "Result value overflow for 0 argument at -1 index.", pol);
                    }
                    else if( z > 0. )
                    {
                         w = boost::math::policies::raise_domain_error<ArgumentType>("lambert_w(%1%,IndexType)", "Argument value %1% out of range for real arguments.", z, pol);
                    }
                    else
                    {
                         w = lambw::_raw<ArgumentType,ArgumentType,IndexType>(z, boost::math::round(k,pol),pol);
                    }
               }
               else
               {
                    w = lambw::_raw<ArgumentType,ArgumentType,IndexType>(z, boost::math::round(k,pol),pol);
               }
          }
          else
          {
               w = boost::math::policies::raise_domain_error<ArgumentType>("lambert_w(%1%,IndexType)", "Invalid index value of %1% for real arguments. Valid values are 0 and -1.", k, pol);
          }

          if(boost::math::isnan(w))
          {
               return boost::math::policies::raise_evaluation_error<ArgumentType>("lambert_w(%1%,IndexType)", "Function return value %1% is Not-a-Number", w, pol);
          }
          else if( !(k == 0 && z == 0.) && w == 0 )
          {
               return boost::math::policies::raise_underflow_error<ArgumentType>("lambert_w(%1%,IndexType)", "Function return value has underflowed to zero.", pol);
          }
          else if( !boost::math::isnormal(w) )
          {
               return boost::math::policies::raise_denorm_error<ArgumentType>("lambert_w(%1%,IndexType)", "Function return value %1% is denormalised.", w, pol);
          }
          else if( boost::math::isinf(w) )
          {
               return boost::math::policies::raise_overflow_error<ArgumentType>("lambert_w(%1%,IndexType)", "Result value overflow.", pol);
          }
          else
          {
               return w;
          }
     }

     //Checks for errors in input arguments and the return value using the policy pol for complex valued arguments.
     template<class ArgumentType, class IndexType, class Policy>
     inline std::complex<ArgumentType> _with_complex_range_checks(const std::complex<ArgumentType>& z, const IndexType& k, const Policy& pol)
     {
          std::complex<ArgumentType> w;

          if(boost::math::isnan(_complex_real(z)) || boost::math::isnan(_complex_imag(z)))
          {
               w = std::complex<ArgumentType>(
                    boost::math::policies::raise_evaluation_error<ArgumentType>(
                         "lambert_w(%1%,IndexType)",
                         "Argument with real value %1% is Not-a-Number!",
                         _complex_real(z),pol),
                    std::numeric_limits<ArgumentType>::quiet_NaN());
          }
          else if(boost::math::isinf(_complex_real(z)) || boost::math::isinf(_complex_imag(z)))
          {
               boost::math::policies::raise_overflow_error<ArgumentType>(
                                "lambert_w(%1%,IndexType)",
                                "Argument value overflow.",pol);
               return std::complex<ArgumentType>(std::numeric_limits<ArgumentType>::quiet_NaN(), std::numeric_limits<ArgumentType>::quiet_NaN());
          }
          else if(k != 0 && _complex_real(z) == 0. && _complex_imag(z) == 0.)
          {
               w = std::complex<ArgumentType>(
                    boost::math::policies::raise_pole_error<ArgumentType>(
                         "lambert_w(std::complex<%1%>,IndexType)",
                         "Function has a pole at %1% for indices other than 0.",
                         _complex_real(z),pol),
                    std::numeric_limits<ArgumentType>::quiet_NaN());
          }
          else
          {
               w = lambw::_raw<std::complex<ArgumentType>,ArgumentType,IndexType,Policy>(z,boost::math::round(k,pol),pol);
          }

          if(boost::math::isnan(_complex_real(w)) || boost::math::isnan(_complex_imag(w)))
          {
               boost::math::policies::raise_evaluation_error<ArgumentType>(
                                "lambert_w(%1%,IndexType)",
                                "Return value with real value %1% is Not-a-Number!",
                                _complex_real(w),pol);
               return w;
          }
          else if(boost::math::isinf(_complex_real(w)) || boost::math::isinf(_complex_imag(w)))
          {
               boost::math::policies::raise_overflow_error<ArgumentType>(
                                "lambert_w(%1%,IndexType)",
                                "Return value overflow.",pol);
               return w;
          }
          else if(!boost::math::isnormal(_complex_real(w)) || !boost::math::isnormal(_complex_imag(w)))
          {
               boost::math::policies::raise_denorm_error<ArgumentType>(
                                "lambert_w(%1%,IndexType)",
                                "Return value with real value %1% is denormalised.",_complex_real(w),pol);
               return w;
          }
          else
          {
               return w;
          }
     }
} // namespace lambw
} // namespace math
} // namespace boost

#endif // BOOST_MATH_LAMBERT_W_RAW_HPP_INCLUDED
