
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_LAMBERT_W_TOOLS_HPP_INCLUDED
#define BOOST_MATH_LAMBERT_W_TOOLS_HPP_INCLUDED

#include <cmath>
#include <complex>
#include <limits>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_complex.hpp>
//#include <boost/typeof/typeof.hpp>
#include <boost/array.hpp>

namespace boost
{
namespace math
{
namespace lambw
{
     //Function overloads returning the imaginary unit for complex types, NaN otherwise.
     template<class ArgumentType>
     inline const ArgumentType& _imaginary_unit(const ArgumentType &z)
     {
          static const ArgumentType i = std::numeric_limits<ArgumentType>::quiet_NaN();

          return i;
     }

     template<class ArgumentType>
     inline const std::complex<ArgumentType>& _imaginary_unit(const std::complex<ArgumentType> &z)
     {
          static const std::complex<ArgumentType> i((ArgumentType)0.,(ArgumentType)1.);

          return i;
     }

     //Function returning exp(-1)
     template<class CoeffType>
     inline const CoeffType& rec_e()
     {
          using std::exp;

          static const CoeffType ans = exp((CoeffType)(-1.));

          return ans;
     }

     //Function overloads to return the real part of any number
     template<class CoeffType>
     inline CoeffType _complex_real(const CoeffType &z)
     {
          return z;
     }

     template<class CoeffType>
     inline CoeffType _complex_real(const std::complex<CoeffType> &z)
     {
          return z.real();
     }

     //Function overloads to return the imaginary part of any number
     template<class CoeffType>
     inline CoeffType _complex_imag(const CoeffType &z)
     {
          static const CoeffType zero = (CoeffType)0.;
          return zero;
     }

     template<class CoeffType>
     inline CoeffType _complex_imag(const std::complex<CoeffType> &z)
     {
          return z.imag();
     }

     //Function overloads to return 1. For some reason when using boost::multiprecision types _series_sum does not compile for sum = (ArgumentType)1.
     template<class CoeffType>
     inline const CoeffType& _complex_one(const CoeffType &z)
     {
          static const CoeffType ans = (CoeffType)1.;
          return ans;
     }

     template<class CoeffType>
     inline const std::complex<CoeffType>& _complex_one(const std::complex<CoeffType> &z)
     {
          static const std::complex<CoeffType> ans((CoeffType)1.,(CoeffType)0.);
          return ans;
     }

     //Evaluates a power series, but it requires the coefficients to be processed as follows: c'_k = c_{k+1}/c_k
     template<class ArgumentType, class CoeffType, std::size_t L>
     ArgumentType _series_sum(const ArgumentType &z, const boost::array<CoeffType,L> &c)
     {
          ArgumentType sum = _complex_one(z);

          for(std::size_t k = 0; k < L; ++k)
          {
               sum = _complex_one(z)+c[L-k-1]*z*sum;
          }

          return sum;
     }

     //Creates an array of the coefficients defined by the function argument coeffs.
     template<class CoeffType, std::size_t L>
     boost::array<CoeffType,L> _coeff_array(CoeffType (&coeffs)(std::size_t) )
     {
          boost::array<CoeffType,L> ans;
          for(std::size_t k = 0; k < L; ++k)
          {
               ans[k] = coeffs(k);
          }
          return ans;
     }

     //constexpr-s in these require C++14. BOOST_CONSTEXPR inserts constexpr-s from C++11, without checking for C++14 compatibility.
     //Cuts up the complex plane on the imaginary axis in 2pi intervals, effectively repeating z-2pi*k for each integer k
     template<class ArgumentType, class Policy>
     /*BOOST_CONSTEXPR*/ ArgumentType _linstrips(const ArgumentType &z, const Policy &pol)
     {
          return z;
     }

     template<class ArgumentType, class Policy>
     /*BOOST_CONSTEXPR*/ std::complex<ArgumentType> _linstrips(const std::complex<ArgumentType> &z, const Policy &pol)
     {
          return z-boost::math::constants::two_pi<ArgumentType>()
               *_imaginary_unit(z)*boost::math::round(_complex_imag(z)/boost::math::constants::two_pi<ArgumentType>(),pol);
     }

     //Tests if w satisfies w*exp(w)=z. In order to avoid overflows it
     //actually checks the equivalent, but numerically less precise log(w)+w=log(z), which is the logarithm of the previous.
     //The relative errors (x1-x2)/x1 and (x1-x2)/x2 can be expressed as 1-x2/x1 and x1/x2-1, which in turn using logarithms:
     //1-exp(log(x2)-log(x1)) and exp(log(x1)-log(x2))-1, or -expm1(log(x2)-log(x1)) and expm1(log(x1)-log(x2)).
     //Since expm1(x) ~ x for small arguments both can equivalently be approximated with log(x1)-log(x2).
     //On that basis this function tests for relative errors.
     template<class ArgumentType, class Policy>
     ArgumentType _test_s(const ArgumentType &z, const ArgumentType &w, const Policy &pol)
     {
          using std::log;
          using std::abs;

          if(_complex_real(z) < 0.)
          {
               return _linstrips(log(-w)+w-log(-z),pol);
          }
          else
          {
               return _linstrips(log(w)+w-log(z),pol);
          }
     }

     template<class ArgumentType, class Policy>
     ArgumentType _test(const ArgumentType &z, const ArgumentType &w, const Policy &pol)
     {
          using std::exp;
          using std::abs;

          ArgumentType t = w*exp(w);

          return std::max(abs((z-t)/z),abs((z-t)/t));
     }

     template<class ArgumentType, class Policy>
     ArgumentType _test(const std::complex<ArgumentType> &z, const std::complex<ArgumentType> &w, const Policy &pol)
     {
          using std::exp;
          using std::abs;

          std::complex<ArgumentType> t = w*exp(w);

          return std::max(abs((z-t)/z),abs((z-t)/t));
     }

     //Abandoned exact evaluation of relative error, using expm1.
     //expm1 is not defined for complex arguments, which is why this was abandoned.
     /*
     template<class ArgumentType>
     ArgumentType _test(const ArgumentType &z, const ArgumentType &w)
     {
          return expm1(_test_s(z,w));
     }
     */
     //Abandoned scoring system
     /*
     template<class ArgumentType, class CoeffType>
     CoeffType _test_score(const ArgumentType &z, const ArgumentType &w)
     {
          using std::abs;
          using std::log10;

          return -log10(abs(_test_s(z,w)));
     }
     */

     //Returns the k-th branch of the complex logarithm: log(z)+i*2pi*k
     template<class ArgumentType, class IndexType>
     ArgumentType _log_k(const ArgumentType &z, const IndexType& k)
     {
          using std::log;

          return (k==0)?(log(z)):(std::numeric_limits<ArgumentType>::quiet_NaN());
     }

     template<class ArgumentType, class IndexType>
     std::complex<ArgumentType> _log_k(const std::complex<ArgumentType> &z, const IndexType& k)
     {
          using std::log;

          return log(z) + (boost::math::constants::two_pi<ArgumentType>()*k)*_imaginary_unit(z);
     }
} //namespace lambw
} //namespace math
} //namespace boost

#endif // BOOST_MATH_LAMBERT_W_TOOLS_HPP_INCLUDED
