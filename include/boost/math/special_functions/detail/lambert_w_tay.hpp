
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_LAMBERT_W_TAY_HPP_INCLUDED
#define BOOST_MATH_LAMBERT_W_TAY_HPP_INCLUDED

#include <boost/array.hpp>
#include "lambert_w_tools.hpp"

namespace boost
{
namespace math
{
namespace lambw
{
     //Constant that controls series evaluation length
     const std::size_t _tay_series_L = 10;

     //Computes the coefficients to use with _series_sum to calculate the Taylor series for W(z) (a.k.a. W_0(z)) at abs(z)<exp(-1)
     //From here on W(z) is used to refer to W_0(z), since it is the most commonly used branch of the function.
     template<class CoeffType>
     /*BOOST_CONSTEXPR*/ CoeffType _tay_coeffs(std::size_t k)
     {
         using std::pow;

         CoeffType n = (CoeffType)(k+1);
         return (-pow(((n+1)/n),n-1));
     }

     //Array to store the Taylor series coefficients.
     template<class CoeffType>
     inline const boost::array<CoeffType,_tay_series_L>& _lw_tc()
     {
          static const boost::array<CoeffType,_tay_series_L> ans = _coeff_array<CoeffType,_tay_series_L>(_tay_coeffs<CoeffType>);
          return ans;
     }

     //Returns an approximation for W(z) for abs(z)<exp(-1)
     template<class ArgumentType,class CoeffType>
     ArgumentType _tay(const ArgumentType &z)
     {
          return z*_series_sum(z,_lw_tc<CoeffType>());
     }
} //namespace lambw
} //namespace math
} //namespace boost

#endif // BOOST_MATH_LAMBERT_W_TAY_HPP_INCLUDED
