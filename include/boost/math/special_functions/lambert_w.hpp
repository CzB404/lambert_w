
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//Sources:
//
//Lambert W function - Wikipedia, the free encyclopedia
//
//Lagrange inversion theorem - Wikipedia, the free encyclopedia
//
//Numerical Evaluation of the Lambert W Function
//and Application to Generation of Generalized
//Gaussian Noise With Exponent 1/2
// - François Chapeau-Blondeau, Member, IEEE, and Abdelilah Monir
//(p2160 IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 50, NO. 9, SEPTEMBER 2002)
//
//Having Fun with Lambert W(x) Function
// - Darko Veberic a,b,c
//   a University of Nova Gorica, Slovenia
//   b IK, Forschungszentrum Karlsruhe, Germany
//   c J. Stefan Institute, Ljubljana, Slovenia
//(arXiv:1003.1628v1 [cs.MS] 8 Mar 2010)

#ifndef BOOST_MATH_LAMBW_HPP_INCLUDED
#define BOOST_MATH_LAMBW_HPP_INCLUDED

#include <cmath>
#include <complex>
#include <limits>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/config.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/array.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/nonfinite_num_facets.hpp>
//#include <boost/math/special_functions/factorials.hpp>

namespace boost
{
namespace math
{
namespace lambw
{
     ///Constants controlling series evaluation lengths
     const std::size_t _tay_series_L = 10;
     const std::size_t _laur_series_L = 10;
     const std::size_t _asy_series_L = 4; //Does not actually have an effect yet
     const std::size_t _mid_series_L = 8; //Does not actually have an effect yet
     const std::size_t _sing_series_L = 9;

     ///Function returning the imaginary unit for complex types, NaN otherwise
     template<typename ArgumentType>
     inline const ArgumentType& _imaginary_unit()
     {
          using std::sqrt;
          using std::abs;

          static ArgumentType z;
          static const ArgumentType i = sqrt(ArgumentType((BOOST_TYPEOF(abs(z)))(-1.)));

          return i;
     }

     ///Function returning exp(-1)
     template<typename CoeffType>
     inline const CoeffType& rec_e()
     {
          using std::exp;

          static const CoeffType ans = exp((CoeffType)(-1.));

          return ans;
     }

     ///Function overloads to return the real part of any number
     template<typename CoeffType>
     inline CoeffType _complex_real(const CoeffType &z)
     {
          return z;
     }

     template<typename CoeffType>
     inline CoeffType _complex_real(const std::complex<CoeffType> &z)
     {
          return z.real();
     }

     ///Function overloads to return the imaginary part of any number
     template<typename CoeffType>
     inline CoeffType _complex_imag(const CoeffType &z)
     {
          return (CoeffType)0.;
     }

     template<typename CoeffType>
     inline CoeffType _complex_imag(const std::complex<CoeffType> &z)
     {
          return z.imag();
     }

     ///Function overloads to return 1. For some reason when using boost::multiprecision types _series_sum does not compile for sum = 1.
     template<typename CoeffType>
     inline const CoeffType& _complex_one(const CoeffType &z)
     {
          static const CoeffType ans = (CoeffType)1.;
          return ans;
     }

     template<typename CoeffType>
     inline const std::complex<CoeffType>& _complex_one(const std::complex<CoeffType> &z)
     {
          static const std::complex<CoeffType> ans((CoeffType)1.,(CoeffType)0.);
          return ans;
     }

     ///Evaluates a power series, but it requires the coefficients to be preprocessed: c'_k = c_{k+1}/c_k
     template<typename ArgumentType, typename CoeffType, std::size_t L>
     ArgumentType _series_sum(const ArgumentType &z, const boost::array<CoeffType,L> &c)
     {
          ArgumentType sum = _complex_one(z);

          for(std::size_t k = 0; k < L; ++k)
          {
               sum = _complex_one(z)+c[L-k-1]*z*sum;
          }

          return sum;
     }

     ///Creates an array of the coefficients defined by the function argument coeffs.
     template<typename CoeffType, std::size_t L>
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
     ///Cuts up the complex plane on the imaginary axis in 2pi intervals, effectively repeating z-2pi*k for each integer k
     template<typename ArgumentType>
     /*BOOST_CONSTEXPR*/ ArgumentType _linstrips(const ArgumentType &z)
     {
          using std::abs;

          if(!(boost::is_complex<ArgumentType>::value) )
          {
               return z;
          }
          else
          {
               return z-boost::math::constants::two_pi<BOOST_TYPEOF(abs(z))>()
               *_imaginary_unit<ArgumentType>()*boost::math::round(_complex_imag(z)/boost::math::constants::two_pi<BOOST_TYPEOF(abs(z))>());
          }
     }

     ///Tests if w satisfies w*exp(w)=z. In order to avoid overflows it
     ///actually checks the equivalent, but numerically less precise log(w)+w=log(z), which is the logarithm of the previous.
     ///The relative errors (x1-x2)/x1 and (x1-x2)/x2 can be expressed as 1-x2/x1 and x1/x2-1, which in turn using logarithms:
     ///1-exp(log(x2)-log(x1)) and exp(log(x1)-log(x2))-1, or -expm1(log(x2)-log(x1)) and expm1(log(x1)-log(x2)).
     ///Since expm1(x) ~ x for small arguments both can equivalently be approximated with log(x1)-log(x2).
     ///On that basis this function tests for relative errors.
     template<typename ArgumentType>
     ArgumentType _test_s(const ArgumentType &z, const ArgumentType &w)
     {
          using std::log;
          using std::abs;

          if(_complex_real(z) < 0.)
          {
               return _linstrips(log(-w)+w-log(-z));
          }
          else
          {
               return _linstrips(log(w)+w-log(z));
          }
     }

     ///Abandoned exact evaluation of relative error, using expm1.
     ///expm1 is not defined for complex arguments, which is why this was abandoned.
     /*
     template<typename ArgumentType>
     ArgumentType _test(const ArgumentType &z, const ArgumentType &w)
     {
          return expm1(_test_s(z,w));
     }
     */
     ///Abandoned scoring system
     /*
     template<typename ArgumentType, typename CoeffType>
     CoeffType _test_score(const ArgumentType &z, const ArgumentType &w)
     {
          using std::abs;
          using std::log10;

          return -log10(abs(_test_s(z,w)));
     }
     */
     ///Computes the coefficients to use with _series_sum to calculate the Taylor series for W(z) (a.k.a. W_0(z)) at abs(z)<exp(-1)
     ///From here on W(z) is used to refer to W_0(z), since it is the most commonly used branch of the function.
     template<typename CoeffType>
     /*BOOST_CONSTEXPR*/ CoeffType _tay_coeffs(std::size_t k)
     {
         using std::pow;

         CoeffType n = (CoeffType)(k+1);
         return (-pow(((n+1)/n),n-1));
     }

     ///Computes the coefficients to use with _series_sum to calculate the Taylor series for f(z)=z/W(z)
     ///which was derived from the Laurent series for W(z)^r, hence the name.
     template<typename CoeffType>
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

     ///Array to store the Taylor series coefficients.
     template<typename CoeffType>
     inline const boost::array<CoeffType,_tay_series_L>& _lw_tc()
     {
          static const boost::array<CoeffType,_tay_series_L> ans = _coeff_array<CoeffType,_tay_series_L>(_tay_coeffs<CoeffType>);
          return ans;
     }

     ///Array to store the Laurent series coefficients.
     template<typename CoeffType>
     inline const boost::array<CoeffType,_laur_series_L>& _lw_lc()
     {
          static const boost::array<CoeffType,_laur_series_L> ans = _coeff_array<CoeffType,_laur_series_L>(_laur_coeffs<CoeffType>);
          return ans;
     }

     ///Returns an approximation for W(z) for abs(z)<exp(-1)
     template<typename ArgumentType,typename CoeffType>
     ArgumentType _tay(const ArgumentType &z)
     {
          return z*_series_sum(z,_lw_tc<CoeffType>());
     }

     ///Returns an approximation for W(z) for abs(z)<exp(-1)
     template<typename ArgumentType,typename CoeffType>
     ArgumentType _laur(const ArgumentType &z)
     {
          return z/_series_sum(z,_lw_lc<CoeffType>());
     }

     ///Calculates Stirling numbers of the first kind.
     ///To be used for the asymptotic approximation.
     template<typename IntegralType>
     /*BOOST_CONSTEXPR*/ IntegralType _stirling(IntegralType n, IntegralType k)
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

     ///Returns the k-th branch of the complex logarithm: log(z)+i*2pi*k
     template<typename ArgumentType, typename IndexType>
     ArgumentType _log_k(const ArgumentType &z, IndexType k)
     {
          using std::log;
          using std::abs;

          if(!(boost::is_complex<ArgumentType>::value) )
          {
               return (k==0)?(log(z)):(std::numeric_limits<ArgumentType>::quiet_NaN());
          }
          else
          {
               return log(z) + (boost::math::constants::two_pi<BOOST_TYPEOF(abs(z))>()*k)*_imaginary_unit<ArgumentType>();
          }
     }

     ///Calculates the asymptotic approximation
     ///Depending on L1 and L2 it either calculates W_k(z) for large abs(z) or W_{+-1}(z) for small abs(z).
     ///W_{+-1}(z), as in the W_1 and W_{-1} have a direct branch cut between them on the real axis at [-exp(-1):0).
     ///This expansion works for complex numbers close to that branch cut. From here on W_{+-1}(z) refers to expansion working in similar fashion.
     template<typename ArgumentType, typename CoeffType>
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

     ///Calculates the asymptotic approximation for W_k(z) at large abs(z)
     template<typename ArgumentType, typename CoeffType, typename IndexType>
     ArgumentType _asy(const ArgumentType &z, IndexType k = 0)
     {
          using std::log;

          ArgumentType L1 = _log_k(z,k);
          ArgumentType L2 = log(L1);

          return _asy_body<ArgumentType,CoeffType>(L1,L2);
     }

     ///Calculates the asymptotic approximation for W_{+-1}(z) at small abs(z)
     template<typename ArgumentType, typename CoeffType>
     ArgumentType _N1(const ArgumentType &z)
     {
          using std::log;

          ArgumentType L1 = log(-z);
          ArgumentType L2 = log(-L1);

          return _asy_body<ArgumentType,CoeffType>(L1,L2);
     }

     ///Array to store the coefficients of the series expansion for W(exp(1+z)).
     ///It has a larger convergence radius than the Taylor and Laurent series, but smaller than the asymptotic.
     template<typename CoeffType>
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

     ///Calculates the approximation of W(z) for mid-range values.
     template<typename ArgumentType, typename CoeffType>
     ArgumentType _mid(const ArgumentType &z)
     {
          using std::log;

          ArgumentType x = log(z)-(CoeffType)1.;
          return _series_sum(x,_lw_midc<CoeffType>());
     }

     ///The following functions calculate the coefficients for the series expansion
     ///of W(z) and W_{+-1}(z) for abs(z-exp(-1))<exp(-1), or in other words: around the singularity at z=-exp(-1)
     template<typename CoeffType>
     /*BOOST_CONSTEXPR*/ CoeffType _sing_coeffs_mu(std::size_t k);

     template<typename CoeffType>
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

     template<typename CoeffType>
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

     template<typename CoeffType>
     /*BOOST_CONSTEXPR*/ CoeffType _sing_coeffs(std::size_t k)
     {
          CoeffType ans = _sing_coeffs_mu<CoeffType>(k+1)/_sing_coeffs_mu<CoeffType>(k);
          return ans;
     }

     ///Array to store the coefficients for the singularity series expansion
     template<typename CoeffType>
     inline const boost::array<CoeffType,_sing_series_L>& _lw_singc()
     {
          static const boost::array<CoeffType,_sing_series_L> ans = _coeff_array<CoeffType,_sing_series_L>(_sing_coeffs<CoeffType>);
          return ans;
     }

     ///Calculates an approximation for W(z) around the singularity at z=-exp(-1) in an exp(-1) radius.
     template<typename ArgumentType, typename CoeffType, typename IndexType>
     ArgumentType _sing(const ArgumentType &z, IndexType k)
     {
          using std::pow;
          using std::sqrt;
          using std::exp;
          using std::log;

          ArgumentType p;
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
          p = ((CoeffType)pow(-1.,k%2))*sqrt((CoeffType)2.*(boost::math::constants::e<CoeffType>()*z+(CoeffType)1.));

          return -_series_sum(p,_lw_singc<CoeffType>());
     }

     ///Fritsch's iteration for the Lambert W function with index k. It is a fast approximation that works for complex arguments,
     ///but preliminary tests using gnuplot revealed numeric weaknesses for poor starting approximations.
     template<typename ArgumentType, typename CoeffType, typename IndexType>
     ArgumentType _iter_frit(const ArgumentType &x, IndexType k, const ArgumentType &Wn)
     {
          using std::log;
          using std::abs;
          using std::exp;

          if((k==0||k==1||k==-1)&&(x==(CoeffType)0.||x==-rec_e<CoeffType>()))
          {
               return Wn;
          }
          else
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

               return Wn*((CoeffType)1.+e_n);
          }
     }

     ///Calculates Newton's method without logarithms for W_k(z).
     ///Newton's method was chosen instead of Halley's because Halley's method tended to overflow for large arguments.
     template<typename ArgumentType, typename CoeffType, typename IndexType>
     ArgumentType _iter_newt_nolog(const ArgumentType &x, const ArgumentType &fn)
     {
          using std::exp;

          return fn*((CoeffType)1.-(CoeffType)1./((CoeffType)1.+fn))+x/((CoeffType)1.+fn)*exp(-fn);
     }

     ///Calculates Newton's method with logarithms. This version works for large arguments.
     template<typename ArgumentType, typename CoeffType, typename IndexType>
     ArgumentType _iter_newt_wlog(const ArgumentType &x, const ArgumentType &fn)
     {
          using std::exp;
          using std::log;

          if(_complex_real(x)<0)
          {
               if(_complex_real(fn)<0)
               {
                    return fn-(CoeffType)1./((CoeffType)1.+(CoeffType)1./fn)+exp(log(-x)-log(-(CoeffType)1.-fn)-(fn));
               }
               else
               {
                    return fn-(CoeffType)1./((CoeffType)1.+(CoeffType)1./fn)-exp(log(-x)-log((CoeffType)1.+fn)-(fn));
               }
          }
          else
          {
               if(_complex_real(fn)<0)
               {
                    return fn-(CoeffType)1./((CoeffType)1.+(CoeffType)1./fn)-exp(log(x)-log(-(CoeffType)1.-fn)-(fn));
               }
               else
               {
                    return fn-(CoeffType)1./((CoeffType)1.+(CoeffType)1./fn)+exp(log(x)-log((CoeffType)1.+fn)-(fn));
               }
          }
     }

     ///Chooses which version of Newton's method to use and returns the value.
     template<typename ArgumentType, typename CoeffType, typename IndexType>
     ArgumentType _iter_newt(const ArgumentType &x, IndexType k, const ArgumentType &fn)
     {
          using std::abs;
          using std::exp;

          if(((((k==0)||(k==-1))&&(_complex_imag(x)>=0))
               ||((k==1)&&(_complex_imag(x)<0)))
               &&(abs(x+rec_e<CoeffType>())<7e-4))
          {
               return fn;
          }
          else
          {
               if(((k==0)&&(abs(x)<rec_e<CoeffType>()))
                    ||(k==-1&&_complex_imag(x)>=0&&abs(x+rec_e<CoeffType>())<0.3)
                    ||(k==1&&_complex_imag(x)<0&&abs(x+rec_e<CoeffType>())<0.3))
               {
                    return _iter_newt_nolog<ArgumentType,CoeffType,IndexType>(x,fn);
               }
               else
               {
                    return _iter_newt_wlog<ArgumentType,CoeffType,IndexType>(x,fn);
               }
          }
     }

     ///Unfinished dynamic approximation selector
     /*
     template<typename ArgumentType, std::size_t L>
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

     ///This function determines which approximation to use as input for iterations, then applies those iterations.
     ///It makes a few Newton's method iterations to refine the input for Fritsch's iteration, which runs until it hits machine precision.
     ///Policies are accepted as arguments, but because this function template accepts both real and complex values,
     ///policy-based errors have not been introduced in here. Instead, there are separate checking function following this one.
     template<typename ArgumentType, typename CoeffType, typename IndexType, class Policy>
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

          CoeffType err = abs(_test_s<ArgumentType>(z,w));

          while(err > 1e-4)
          {
               w = _iter_newt<ArgumentType,CoeffType,IndexType>(z,k,w);
               err = abs(_test_s<ArgumentType>(z,w));
          }

          CoeffType err_prev;
          ArgumentType w_prev;

          do
          {
               err_prev = err;
               w_prev = w;
               w = _iter_frit<ArgumentType,CoeffType,IndexType>(z,k,w);
               err = abs(_test_s<ArgumentType>(z,w));
          }
          while(err < err_prev);

          return w_prev;
     }

     ///Checks for errors in input arguments and the return value using the policy pol for real valued arguments.
     template<typename ArgumentType, typename IndexType, class Policy>
     inline ArgumentType _with_real_range_checks(const ArgumentType& z, const IndexType& k, const Policy& pol)
     {
          ArgumentType w;

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

     ///Checks for errors in input arguments and the return value using the policy pol for complex valued arguments.
     template<typename ArgumentType, typename IndexType, class Policy>
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
          else if(k != 0 && _complex_real(z) == 0. && _complex_imag(z) == 0.)
          {
               w = std::complex<ArgumentType>(boost::math::policies::raise_pole_error<ArgumentType>("lambert_w(std::complex<%1%>,IndexType)", "Function has a pole at %1% for indices other than 0.",_complex_real(z),pol),std::numeric_limits<ArgumentType>::quiet_NaN());
          }
          else
          {
               w = lambw::_raw<std::complex<ArgumentType>,ArgumentType,IndexType,Policy>(z,boost::math::round(k,pol),pol);
          }

          ///Add result checks later!!!
          return w;
     }
} // namespace lambw

     ///Returns W(z) with the default policy for real argument z.
     template<typename ArgumentType>
     inline ArgumentType lambert_w(const ArgumentType& z)
     {
          return lambw::_with_real_range_checks<ArgumentType,int>(z,0,policies::policy<>());
     }

     ///Returns W_k(z) with the default policy for real argument z.
     template<typename ArgumentType, typename IndexType>
     inline ArgumentType lambert_w(const ArgumentType& z, const IndexType& k)
     {
          return lambw::_with_real_range_checks<ArgumentType,IndexType>(z,k,policies::policy<>());
     }

     ///Returns W_k(z) with the policy pol for real argument z.
     template<typename ArgumentType, typename IndexType, class Policy>
     inline ArgumentType lambert_w(const ArgumentType& z, const IndexType& k, const Policy& pol)
     {
          BOOST_FPU_EXCEPTION_GUARD

          //typedef typename tools::promote_args<T, U>::type ArgumentType result_type;
          typedef ArgumentType result_type;

          typedef typename policies::evaluation<result_type, Policy>::type value_type;

          typedef typename policies::normalise<
               Policy,
               policies::promote_float<false>,
               policies::promote_double<false>,
               policies::discrete_quantile<>,
               policies::assert_undefined<> >::type forwarding_policy;

          return policies::checked_narrowing_cast<result_type, forwarding_policy>(
               lambw::_with_real_range_checks<ArgumentType,IndexType,Policy>(
                  static_cast<value_type>(z),
                  static_cast<IndexType>(k),
                  forwarding_policy()),
               "boost::math::lambert_w<%1%,IndexType>(%1%,IndexType)");
     }

     ///Returns W(z) with the default policy for complex argument z.
     template<typename ArgumentType>
     inline std::complex<ArgumentType> lambert_w(const std::complex<ArgumentType> &z)
     {
          return lambw::_with_complex_range_checks<ArgumentType,int>(z, 0,policies::policy<>());
     }

     ///Return W_k(z) with the default policy for complex argument z.
     template<typename ArgumentType, typename IndexType>
     inline std::complex<ArgumentType> lambert_w(const std::complex<ArgumentType> &z, const IndexType& k)
     {
          return lambw::_with_complex_range_checks<ArgumentType,IndexType>(z, k,policies::policy<>());
     }

     ///Return W_k(z) with the policy pol for complex argument z.
     template<typename ArgumentType, typename IndexType, class Policy>
     inline std::complex<ArgumentType> lambert_w(const std::complex<ArgumentType>& z, const IndexType& k, const Policy& pol)
     {
          BOOST_FPU_EXCEPTION_GUARD

          //typedef typename tools::promote_args<T, U>::type ArgumentType result_type;
          typedef ArgumentType result_type;

          typedef typename policies::evaluation<result_type, Policy>::type value_type;

          typedef typename policies::normalise<
               Policy,
               policies::promote_float<false>,
               policies::promote_double<false>,
               policies::discrete_quantile<>,
               policies::assert_undefined<> >::type forwarding_policy;

          return policies::checked_narrowing_cast<result_type, forwarding_policy>(
               lambw::_with_complex_range_checks<ArgumentType,IndexType,Policy>(
                  static_cast<value_type>(z),
                  static_cast<IndexType>(k),
                  forwarding_policy()),
               "boost::math::lambert_w<%1%,IndexType>(%1%,IndexType)");
     }
} // namespace math
} // namespace boost

#endif // BOOST_MATH_LAMBW_HPP_INCLUDED
