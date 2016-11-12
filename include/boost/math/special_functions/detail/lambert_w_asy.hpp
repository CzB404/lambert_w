
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_LAMBERT_W_ASY_HPP_INCLUDED
#define BOOST_MATH_LAMBERT_W_ASY_HPP_INCLUDED

//#include <boost/math/special_functions/factorials.hpp>
#include "lambert_w_tools.hpp"

namespace boost
{
namespace math
{
namespace lambw
{
     //Constant that controls series evaluation length
     const std::size_t _asy_series_L = 4; //Does not actually have an effect yet

     //Calculates Stirling numbers of the first kind.
     //To be used for the asymptotic approximation.
     /*
     template<class IntegralType>
     BOOST_CONSTEXPR IntegralType _stirling(IntegralType n, IntegralType k)
     {
          if( n == 0 || k == 0)
          {
               if( n == 0 && k == 0)
               {
                    return 1;
               }
               else
               {
                    return 0;
               }
          }
          else
          {
               return (n-1)*_stirling((n-1),k)+_stirling((n-1),k-1);
          }
     }
     */

     //Calculates the asymptotic approximation
     //Depending on L1 and L2 it either calculates W_k(z) for large abs(z) or W_{+-1}(z) for small abs(z).
     //W_{+-1}(z), as in the W_1 and W_{-1} have a direct branch cut between them on the real axis at [-exp(-1):0).
     //This expansion works for complex numbers close to that branch cut. From here on W_{+-1}(z) refers to expansion working in similar fashion.
     template<class ArgumentType, class CoeffType>
     ArgumentType _asy_body(const ArgumentType &L1, const ArgumentType &L2)
     {
          ArgumentType ans = L1-L2;
          /*
          for(std::size_t l = 0; l < _asy_series_L; ++l)
          {
               for(std::size_t m = 1; m <= l+1; ++m)
               {
                    ans+=std::pow(-1.,l%2)*_stirling(l+m,l+1)/boost::math::factorial(m)*std::pow(L1,-(CoeffType)l-(CoeffType)m)*std::pow(L2,(CoeffType)m);
               }
          }
          */
          ans += L2/L1 + L2*(((CoeffType)(-2.))+L2)/(CoeffType)2./L1/L1 + L2*((CoeffType)6.-(CoeffType)9.*L2+(CoeffType)2.*L2*L2)/(CoeffType)6./L1/L1/L1
               + L2*((CoeffType)(-12.)+(CoeffType)36.*L2-(CoeffType)22.*L2*L2+(CoeffType)3.*L2*L2*L2)/(CoeffType)12./L1/L1/L1/L1;


          return ans;
     }

     //Calculates the asymptotic approximation for W_k(z) at large abs(z)
     template<class ArgumentType, class CoeffType, class IndexType>
     ArgumentType _asy(const ArgumentType &z, IndexType k = 0)
     {
          using std::log;

          ArgumentType L1 = _log_k(z,k);
          ArgumentType L2 = log(L1);

          return _asy_body<ArgumentType,CoeffType>(L1,L2);
     }

     //Calculates the asymptotic approximation for W_{+-1}(z) at small abs(z)
     template<class ArgumentType, class CoeffType>
     ArgumentType _N1(const ArgumentType &z)
     {
          using std::log;

          ArgumentType L1 = log(-z);
          ArgumentType L2 = log(-L1);

          return _asy_body<ArgumentType,CoeffType>(L1,L2);
     }
} //namespace lambw
} //namespace math
} //namespace boost

#endif // BOOST_MATH_LAMBERT_W_ASY_HPP_INCLUDED
