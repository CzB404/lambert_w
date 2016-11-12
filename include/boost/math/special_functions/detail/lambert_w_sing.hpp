
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_LAMBERT_W_SING_HPP_INCLUDED
#define BOOST_MATH_LAMBERT_W_SING_HPP_INCLUDED

#include <boost/array.hpp>
#include "lambert_w_tools.hpp"

namespace boost
{
namespace math
{
namespace lambw
{
     //Constant that controls series evaluation length
     const std::size_t _sing_series_L = 9;

     //The following functions calculate the coefficients for the series expansion
     //of W(z) and W_{+-1}(z) for abs(z-exp(-1))<exp(-1), or in other words: around the singularity at z=-exp(-1)
     template<class CoeffType>
     /*BOOST_CONSTEXPR*/ CoeffType _sing_coeffs_mu(std::size_t k);

     template<class CoeffType>
     /*BOOST_CONSTEXPR*/ CoeffType _sing_coeffs_alpha(std::size_t k)
     {
          if(k == 0)
          {
               return 2.;
          }
          else if(k == 1)
          {
               return -1.;
          }
          else
          {
               CoeffType sum = 0;
               for(std::size_t j = 2; j <= (k-1); ++j)
               {
                    sum+=_sing_coeffs_mu<CoeffType>(j)*_sing_coeffs_mu<CoeffType>(k+1-j);
               }
               return sum;
          }
     }

     template<class CoeffType>
     /*BOOST_CONSTEXPR*/ CoeffType _sing_coeffs_mu(std::size_t k)
     {
          if(k == 0)
          {
               return -1.;
          }
          else if(k == 1)
          {
               return 1.;
          }
          else
          {
               return (_sing_coeffs_mu<CoeffType>(k-2)/2
                         +_sing_coeffs_alpha<CoeffType>(k-2)/4)*(k-1)/(k+1)
                         -_sing_coeffs_alpha<CoeffType>(k)/2
                         -_sing_coeffs_mu<CoeffType>(k-1)/(k+1);
          }
     }

     template<class CoeffType>
     /*BOOST_CONSTEXPR*/ CoeffType _sing_coeffs(std::size_t k)
     {
          CoeffType ans = _sing_coeffs_mu<CoeffType>(k+1)/_sing_coeffs_mu<CoeffType>(k);
          return ans;
     }

     //Array to store the coefficients for the singularity series expansion
     template<class CoeffType>
     inline const boost::array<CoeffType,_sing_series_L>& _lw_singc()
     {
          static const boost::array<CoeffType,_sing_series_L> ans = _coeff_array<CoeffType,_sing_series_L>(_sing_coeffs<CoeffType>);
          return ans;
     }

     //Calculates an approximation for W(z) around the singularity at z=-exp(-1) in an exp(-1) radius.
     template<class ArgumentType, class CoeffType, class IndexType>
     ArgumentType _sing(const ArgumentType &z, IndexType k)
     {
          using std::pow;
          using std::sqrt;
          using std::exp;
          using std::log;

          ArgumentType p;

          p = ((CoeffType)pow(-1.,k%2))*sqrt((CoeffType)2.*(boost::math::constants::e<CoeffType>()*z+(CoeffType)1.));

          /*
          if(_complex_real(z) >= 0)
          {
               p = ((CoeffType)pow(-1.,k%2))*sqrt((CoeffType)2.*(exp((CoeffType)1.+log(z))+(CoeffType)1.));
          }
          else
          {
               p = ((CoeffType)pow(-1.,k%2))*sqrt((CoeffType)2.*(-exp((CoeffType)1.+log(-z))+(CoeffType)1.));
          }
          */

          return -_series_sum(p,_lw_singc<CoeffType>());
     }
} //namespace lambw
} //namespace math
} //namespace boost

#endif // BOOST_MATH_LAMBERT_W_SING_HPP_INCLUDED
