
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_LAMBERT_W_MID_HPP_INCLUDED
#define BOOST_MATH_LAMBERT_W_MID_HPP_INCLUDED

#include <cmath>
#include <boost/array.hpp>
#include "lambert_w_tools.hpp"

namespace boost
{
namespace math
{
namespace lambw
{
     //Constant that controls series evaluation length
     const std::size_t _mid_series_L = 8; //Does not actually have an effect yet

     //Array to store the coefficients of the series expansion for W(exp(1+z)).
     //It has a larger convergence radius than the Taylor and Laurent series, but smaller than the asymptotic.
     template<class CoeffType>
     inline const boost::array<CoeffType,_mid_series_L>& _lw_midc()
     {
          static const boost::array<CoeffType,_mid_series_L> ans = {
               ((CoeffType)1./(CoeffType)2.),
               ((CoeffType)1./(CoeffType)8.),
               -((CoeffType)1./(CoeffType)12.),
               ((CoeffType)1./(CoeffType)16.),
               -((CoeffType)13./(CoeffType)20.),
               -((CoeffType)47./(CoeffType)312.),
               ((CoeffType)73./(CoeffType)1316.),
               -((CoeffType)2447./(CoeffType)2336.)
          };

          return ans;
     }

     //Calculates the approximation of W(z) for mid-range values.
     template<class ArgumentType, class CoeffType>
     ArgumentType _mid(const ArgumentType &z)
     {
          using std::log;

          ArgumentType x = log(z)-(CoeffType)1.;
          return _series_sum(x,_lw_midc<CoeffType>());
     }
} //namespace lambw
} //namespace math
} //namespace boost

#endif // BOOST_MATH_LAMBERT_W_MID_HPP_INCLUDED
