
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_LAMBERT_W_LAUR_HPP_INCLUDED
#define BOOST_MATH_LAMBERT_W_LAUR_HPP_INCLUDED

#include <boost/array.hpp>
#include "lambert_w_tools.hpp"

namespace boost
{
namespace math
{
namespace lambw
{
     //Constant that controls series evaluation length
     const std::size_t _laur_series_L = 10;

     //Computes the coefficients to use with _series_sum to calculate the Taylor series for f(z)=z/W(z)
     //which was derived from the Laurent series for W(z)^r, hence the name.
     template<class CoeffType>
     /*BOOST_CONSTEXPR*/ CoeffType _laur_coeffs(std::size_t k)
     {
          using std::pow;

          CoeffType n = (CoeffType)k;
          switch(k)
          {
          case 0:
               return 1.;
               break;
          case 1:
               return -0.5;
               break;
          default:
               return pow(n/(n-(CoeffType)1.),n)*((CoeffType)1.-n)/(n+(CoeffType)1.);
               break;
          }
     }

     //Array to store the Laurent series coefficients.
     template<class CoeffType>
     inline const boost::array<CoeffType,_laur_series_L>& _lw_lc()
     {
          static const boost::array<CoeffType,_laur_series_L> ans = _coeff_array<CoeffType,_laur_series_L>(_laur_coeffs<CoeffType>);
          return ans;
     }

     //Returns an approximation for W(z) for abs(z)<exp(-1)
     template<class ArgumentType,class CoeffType>
     ArgumentType _laur(const ArgumentType &z)
     {
          return z/_series_sum(z,_lw_lc<CoeffType>());
     }
} //namespace lambw
} //namespace math
} //namespace boost

#endif // BOOST_MATH_LAMBERT_W_LAUR_HPP_INCLUDED
