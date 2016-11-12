
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_LAMBERT_W_ITER_HPP_INCLUDED
#define BOOST_MATH_LAMBERT_W_ITER_HPP_INCLUDED

#include <cmath>
#include "lambert_w_tools.hpp"

namespace boost
{
namespace math
{
namespace lambw
{
     //Fritsch's iteration for the Lambert W function with index k. It is a fast approximation that works for complex arguments,
     //but preliminary tests using gnuplot revealed numeric weaknesses for poor starting approximations.
     template<class ArgumentType, class CoeffType, class IndexType>
     void _iter_frit(const ArgumentType &x, const IndexType& k, ArgumentType &Wn)
     {
          using std::log;
          using std::abs;
          using std::exp;

          if((k!=0 && k!=-1) || (x!=(CoeffType)0. && x!=(ArgumentType)(-rec_e<CoeffType>())))
          {
               ArgumentType z_n;

               if(k==0)
               {
                    if(_complex_imag(x)==0
                         &&(_complex_real(x)<=0)
                         &&(_complex_real(x)>=-rec_e<CoeffType>()) )
                    {
                         z_n = log(-x)-log(-Wn)-Wn;
                    }
                    else
                    {
                         z_n = log(x)-log(Wn)-Wn;
                    }
               }
               else
               {
                    if((k==-1)&&(_complex_imag(x)==0)
                         &&(_complex_real(x)<=0)
                         &&(_complex_real(x)>=-rec_e<CoeffType>()))
                    {
                         z_n = log(-x)-log(-Wn)-Wn;
                    }
                    else
                    {
                         z_n = _log_k(x,k)+log((CoeffType)1./Wn)-Wn;
                    }
               }

               ArgumentType q_n=(CoeffType)2.*((CoeffType)1.+Wn)*((CoeffType)1.+Wn+(CoeffType)2./(CoeffType)3.*z_n);

               ArgumentType e_n=z_n/((CoeffType)1.+Wn)*(q_n-z_n)/(q_n-(CoeffType)2.*z_n);

               Wn *= ((CoeffType)1.+e_n);
          }
     }

     //Calculates Newton's method without logarithms for W_k(z).
     //Newton's method was chosen instead of Halley's because Halley's method tended to overflow for large arguments.
     template<class ArgumentType, class CoeffType, class IndexType>
     void _iter_newt_nolog(const ArgumentType &x, ArgumentType &fn)
     {
          using std::exp;

          fn += -fn/((CoeffType)1.+fn)+x/((CoeffType)1.+fn)*exp(-fn);
     }

     //Calculates Newton's method with logarithms. This version works for large arguments.
     template<class ArgumentType, class CoeffType, class IndexType>
     void _iter_newt_wlog(const ArgumentType &x, ArgumentType &fn)
     {
          using std::exp;
          using std::log;

          if(_complex_real(x)<0)
          {
               if(_complex_real(fn)<0)
               {
                    fn += -fn/((CoeffType)1.+fn)+exp(log(-x)-log(-(CoeffType)1.-fn)-(fn));
               }
               else
               {
                    fn += -fn/((CoeffType)1.+fn)-exp(log(-x)-log((CoeffType)1.+fn)-(fn));
               }
          }
          else
          {
               if(_complex_real(fn)<0)
               {
                    //return fn-(CoeffType)1./((CoeffType)1.+(CoeffType)1./fn)-exp(log(x)-log(-(CoeffType)1.-fn)-(fn));
                    fn += -fn/((CoeffType)1.+fn)-exp(log(x)-log(-(CoeffType)1.-fn)-(fn));
               }
               else
               {
                    fn += -fn/((CoeffType)1.+fn)+exp(log(x)-log((CoeffType)1.+fn)-(fn));
               }
          }
     }

     //Chooses which version of Newton's method to use and returns the value.
     template<class ArgumentType, class CoeffType, class IndexType>
     void _iter_newt(const ArgumentType &x, const IndexType& k, ArgumentType &fn)
     {
          using std::abs;
          using std::exp;

          if(((((k==0)||(k==-1))&&(_complex_imag(x)>=0))
               ||((k==1)&&(_complex_imag(x)<0)))
               &&(abs(x+rec_e<CoeffType>())<7e-4))
          {

          }
          else
          {
               if(((k==0)&&(abs(x)<rec_e<CoeffType>()))
                    ||(k==-1&&_complex_imag(x)>=0&&abs(x+rec_e<CoeffType>())<0.3)
                    ||(k==1&&_complex_imag(x)<0&&abs(x+rec_e<CoeffType>())<0.3))
               {
                    _iter_newt_nolog<ArgumentType,CoeffType,IndexType>(x,fn);
               }
               else
               {
                    _iter_newt_wlog<ArgumentType,CoeffType,IndexType>(x,fn);
               }
          }
     }
} //namespace lambw
} //namespace math
} //namespace boost

#endif // BOOST_MATH_LAMBERT_W_ITER_HPP_INCLUDED
