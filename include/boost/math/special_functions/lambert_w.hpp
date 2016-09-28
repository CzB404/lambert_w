
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

//#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <boost/math/special_functions/fpclassify.hpp>
//#include <type_traits>
#include <boost/math/special_functions/round.hpp>
#include <boost/config.hpp>
//#include <boost/type_traits/is_arithmetic.hpp>
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
	const std::size_t _lambw_series_L = 10;
	const std::size_t _lambw_asy_series_L = 4;
	const std::size_t _lambw_mid_series_L = 8;
	const std::size_t _lambw_sing_series_L = 10;

	template<typename CoeffType>
	const CoeffType& _lambw_imaginary_unit()
	{
		using std::sqrt;
		using std::abs;

		static CoeffType z;
		static const CoeffType i = sqrt(CoeffType((BOOST_TYPEOF(abs(z)))(-1.)));

		return i;
	}

	template<typename CoeffType>
	CoeffType _lambw_complex_real(const CoeffType &z)
	{
		return z;
	}

	template<typename CoeffType>
	CoeffType _lambw_complex_real(const std::complex<CoeffType> &z)
	{
		return z.real();
	}

	template<typename CoeffType>
	CoeffType _lambw_complex_imag(const CoeffType &z)
	{
		return (CoeffType)0.;
	}

	template<typename CoeffType>
	CoeffType _lambw_complex_imag(const std::complex<CoeffType> &z)
	{
		return z.imag();
	}

	template<typename CoeffType>
	const CoeffType& _lambw_complex_one(const CoeffType &z)
	{
		static const CoeffType ans = (CoeffType)1.;
		return ans;
	}

	template<typename CoeffType>
	const std::complex<CoeffType>& _lambw_complex_one(const std::complex<CoeffType> &z)
	{
		static const std::complex<CoeffType> ans((CoeffType)1.,(CoeffType)0.);
		return ans;
	}

	template<typename ArgumentType, typename CoeffType, std::size_t L>
	ArgumentType _lambw_series_sum(const ArgumentType &z, const boost::array<CoeffType,L> &c)
	{
		ArgumentType sum = _lambw_complex_one(z);

		for(std::size_t k = 0; k < L; ++k)
		{
			sum = ((CoeffType)1.)+c[L-k-1]*z*sum;
		}

		return sum;
	}

	template<typename CoeffType, std::size_t L>
	boost::array<CoeffType,L> _lambw_coeff_array(CoeffType (&coeffs)(std::size_t) )
	{
		boost::array<CoeffType,L> ans;
		for(std::size_t k = 0; k < L; ++k)
		{
			ans[k] = coeffs(k);
		}
		return ans;
	}

	//constexpr-s in these require C++14
    template<typename ArgumentType>
	/*BOOST_CONSTEXPR*/ ArgumentType _lambw_linstrips(const ArgumentType &z)
	{
		using std::abs;

		if(!boost::is_complex<ArgumentType>::value)
		{
			return z;
		}
		else
		{
			return z-boost::math::constants::two_pi<BOOST_TYPEOF(abs(z))>()
			*_lambw_imaginary_unit<ArgumentType>()*boost::math::round(_lambw_complex_imag(z)/boost::math::constants::two_pi<BOOST_TYPEOF(abs(z))>());
		}
	}

	template<typename ArgumentType>
	ArgumentType _lambw_test_s(const ArgumentType &z, const ArgumentType &w)
	{
		using std::log;
		using std::abs;

		if(_lambw_complex_real(z) < 0.)
		{
			return _lambw_linstrips(log(-w)+w-log(-z));
		}
		else
		{
			return _lambw_linstrips(log(w)+w-log(z));
		}
	}

	/*
	template<typename ArgumentType>
	ArgumentType _lambw_test(const ArgumentType &z, const ArgumentType &w)
	{
		return expm1(_lambw_test_s(z,w));
	}
	*/

	template<typename ArgumentType, typename CoeffType>
	CoeffType _lambw_test_score(const ArgumentType &z, const ArgumentType &w)
	{
		using std::abs;
		using std::log10;

		return -log10(abs(_lambw_test_s(z,w)));

		//Abandoned integer scoring
		/*
		if(score == std::numeric_limits<BOOST_TYPEOF(abs(z))>::infinity())
		{
			return std::numeric_limits<int>::max();
		}
		else
		{
			return (int)round(score);
		}
		*/
	}

	template<typename CoeffType>
    /*BOOST_CONSTEXPR*/ CoeffType _lambw_tay_coeffs(std::size_t k)
    {
    	using std::pow;

    	CoeffType n = (CoeffType)(k+1);
    	return (-pow(((n+1)/n),n-1));
    }

    template<typename CoeffType>
    /*BOOST_CONSTEXPR*/ CoeffType _lambw_laur_coeffs(std::size_t k)
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

	template<typename CoeffType>
	const boost::array<CoeffType,_lambw_series_L>& _lw_tc()
	{
		static const boost::array<CoeffType,_lambw_series_L> ans = _lambw_coeff_array<CoeffType,_lambw_series_L>(_lambw_tay_coeffs<CoeffType>);
		return ans;
	}

	template<typename CoeffType>
	const boost::array<CoeffType,_lambw_series_L>& _lw_lc()
	{
		static const boost::array<CoeffType,_lambw_series_L> ans = _lambw_coeff_array<CoeffType,_lambw_series_L>(_lambw_laur_coeffs<CoeffType>);
		return ans;
	}

	template<typename ArgumentType,typename CoeffType>
	ArgumentType _lambw_tay(const ArgumentType &z)
	{
		return z*_lambw_series_sum(z,_lw_tc<CoeffType>());
	}

	template<typename ArgumentType,typename CoeffType>
	ArgumentType _lambw_laur(const ArgumentType &z)
	{
		return z/_lambw_series_sum(z,_lw_lc<CoeffType>());
	}

	template<typename IntegralType>
	/*BOOST_CONSTEXPR*/ IntegralType _lambw_stirling(IntegralType n, IntegralType k)
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

	template<typename ArgumentType, typename IndexType>
	ArgumentType _log_k(const ArgumentType &z, IndexType k)
	{
		using std::log;
		using std::abs;

		if(!boost::is_complex<ArgumentType>::value)
		{
			return (k==0)?(log(z)):(std::numeric_limits<ArgumentType>::quiet_NaN());
		}
		else
		{
			return log(z) + (boost::math::constants::two_pi<BOOST_TYPEOF(abs(z))>()*k)*_lambw_imaginary_unit<ArgumentType>();
		}
	}

	template<typename ArgumentType, typename CoeffType>
	ArgumentType _lambw_asy_body(const ArgumentType &L1, const ArgumentType &L2)
	{
		ArgumentType ans = L1-L2;
		/*
		for(std::size_t l = 0; l < _lambw_asy_series_L; ++l)
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

	template<typename ArgumentType, typename CoeffType, typename IndexType>
	ArgumentType _lambw_asy(const ArgumentType &z, IndexType k = 0)
	{
		using std::log;

		ArgumentType L1 = _log_k(z,k);
		ArgumentType L2 = log(L1);

		return _lambw_asy_body<ArgumentType,CoeffType>(L1,L2);
	}

	template<typename ArgumentType, typename CoeffType>
	ArgumentType _lambw_N1(const ArgumentType &z)
	{
		using std::log;

		ArgumentType L1 = log(-z);
		ArgumentType L2 = log(-L1);

		return _lambw_asy_body<ArgumentType,CoeffType>(L1,L2);
	}

	template<typename CoeffType>
	const boost::array<CoeffType,_lambw_series_L>& _lw_midc()
	{
		static const boost::array<CoeffType,_lambw_series_L> ans = {
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

	template<typename ArgumentType, typename CoeffType>
	ArgumentType _lambw_mid(const ArgumentType &z)
	{
		using std::log;

		ArgumentType x = log(z)-(CoeffType)1.;
		return _lambw_series_sum(x,_lw_midc<CoeffType>());
	}

	template<typename CoeffType>
	/*BOOST_CONSTEXPR*/ CoeffType _lambw_sing_coeffs_mu(std::size_t k);

	template<typename CoeffType>
	/*BOOST_CONSTEXPR*/ CoeffType _lambw_sing_coeffs_alpha(std::size_t k)
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
				sum+=_lambw_sing_coeffs_mu<CoeffType>(j)*_lambw_sing_coeffs_mu<CoeffType>(k+1-j);
			}
			return sum;
		}
	}

	template<typename CoeffType>
	/*BOOST_CONSTEXPR*/ CoeffType _lambw_sing_coeffs_mu(std::size_t k)
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
			return (_lambw_sing_coeffs_mu<CoeffType>(k-2)/2
					+_lambw_sing_coeffs_alpha<CoeffType>(k-2)/4)*(k-1)/(k+1)
					-_lambw_sing_coeffs_alpha<CoeffType>(k)/2
					-_lambw_sing_coeffs_mu<CoeffType>(k-1)/(k+1);
		}
	}

	template<typename CoeffType>
	/*BOOST_CONSTEXPR*/ CoeffType _lambw_sing_coeffs(std::size_t k)
	{
		CoeffType ans = _lambw_sing_coeffs_mu<CoeffType>(k+1)/_lambw_sing_coeffs_mu<CoeffType>(k);
		return ans;
	}

	template<typename CoeffType>
	const boost::array<CoeffType,_lambw_series_L>& _lw_singc()
	{
		static const boost::array<CoeffType,_lambw_series_L> ans = _lambw_coeff_array<CoeffType,_lambw_series_L>(_lambw_sing_coeffs<CoeffType>);
		return ans;
	}

	template<typename ArgumentType, typename CoeffType, typename IndexType>
	ArgumentType _lambw_sing(const ArgumentType &z, IndexType k)
	{
		using std::pow;
		using std::sqrt;
		using std::exp;
		using std::log;

		ArgumentType p;
		if(_lambw_complex_real(z) >= 0)
		{
			p = ((CoeffType)pow(-1.,k%2))*sqrt((CoeffType)2.*(exp((CoeffType)1.+log(z))+(CoeffType)1.));
		}
		else
		{
			p = ((CoeffType)pow(-1.,k%2))*sqrt((CoeffType)2.*(-exp((CoeffType)1.+log(-z))+(CoeffType)1.));
		}
		return -_lambw_series_sum(p,_lw_singc<CoeffType>());
	}

	template<typename ArgumentType, typename CoeffType, typename IndexType>
	ArgumentType _lambw_iter_frit(const ArgumentType &x, IndexType k, const ArgumentType &Wn)
	{
		using std::log;
		using std::abs;
		using std::exp;

		if((k==0||k==1||k==-1)&&(x==(CoeffType)0.||x==-exp((CoeffType)(-1.))))
		{
			return Wn;
		}
		else
		{
			ArgumentType z_n;

			if(k==0)
			{
				if(_lambw_complex_imag(x)==0
					&&(_lambw_complex_real(x)<=0)
					&&(_lambw_complex_real(x)>=-exp((CoeffType)(-1.))) )
				{
					z_n = log(abs(x))-log(abs(Wn))-Wn;
				}
				else
				{
					z_n = log(x)-log(Wn)-Wn;
				}
			}
			else
			{
				if((k==-1)&&(_lambw_complex_imag(x)==0)
					&&(_lambw_complex_real(x)<=0)
					&&(_lambw_complex_real(x)>=-(CoeffType)std::exp(-1.l)))
				{
					z_n = log(abs(x))-log(abs(Wn))-Wn;
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

	template<typename ArgumentType, typename CoeffType, typename IndexType>
	ArgumentType _lambw_iter_newt_nolog(const ArgumentType &x, const ArgumentType &fn)
	{
		using std::exp;

		return fn*((CoeffType)1.-(CoeffType)1./((CoeffType)1.+fn))+x/((CoeffType)1.+fn)*exp(-fn);
	}

	template<typename ArgumentType, typename CoeffType, typename IndexType>
	ArgumentType _lambw_iter_newt_wlog(const ArgumentType &x, const ArgumentType &fn)
	{
		using std::exp;
		using std::log;

		if(_lambw_complex_real(x)<0)
		{
			if(_lambw_complex_real(fn)<0)
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
			if(_lambw_complex_real(fn)<0)
			{
				return fn-(CoeffType)1./((CoeffType)1.+(CoeffType)1./fn)-exp(log(x)-log(-(CoeffType)1.-fn)-(fn));
			}
			else
			{
				return fn-(CoeffType)1./((CoeffType)1.+(CoeffType)1./fn)+exp(log(x)-log((CoeffType)1.+fn)-(fn));
			}
		}
	}

	template<typename ArgumentType, typename CoeffType, typename IndexType>
	ArgumentType _lambw_iter_newt(const ArgumentType &x, IndexType k, const ArgumentType &fn)
	{
		using std::abs;
		using std::exp;

		if(((((k==0)||(k==-1))&&(_lambw_complex_imag(x)>=0))
			||((k==1)&&(_lambw_complex_imag(x)<0)))
			&&(abs(x+exp((CoeffType)(-1.)))<7e-4))
		{
			return fn;
		}
		else
		{
			if(((k==0)&&(abs(x)<exp((CoeffType)(-1.))))
				||(k==-1&&_lambw_complex_imag(x)>=0&&abs(x+exp((CoeffType)(-1.)))<0.3)
				||(k==1&&_lambw_complex_imag(x)<0&&abs(x+exp((CoeffType)(-1.)))<0.3))
			{
				return _lambw_iter_newt_nolog<ArgumentType,CoeffType,IndexType>(x,fn);
			}
			else
			{
				return _lambw_iter_newt_wlog<ArgumentType,CoeffType,IndexType>(x,fn);
			}
		}
	}

	///Unfinished dynamic approximation selector
	/*
	template<typename ArgumentType, std::size_t L>
	ArgumentType _lambw_best_approx(const ArgumentType &z, const boost::array<ArgumentType,L> &approxes)
	{
		int prec;
		int prec_max = _lambw_test_score(z,approxes[0]);
		std::size_t max_place = 0;

		for(std::size_t k = 1; k < approxes.size(); ++k)
		{
			prec = _lambw_test_score(z,approxes[k]);
			if(prec_max < prec)
			{
				prec_max = prec;
				max_place = k;
			}
		}
		return approxes[max_place];
	}
	*/

	template<typename ArgumentType, typename CoeffType, typename IndexType>
	ArgumentType _lambw_raw(const ArgumentType &z, IndexType k = 0)
	{
		using std::abs;
		using std::exp;

		ArgumentType w;

		CoeffType z_imag = _lambw_complex_imag(z);

		if(k==-1)
		{
			if(abs(z+(CoeffType)0.04)<=0.14 && z_imag>=0)
			{
				w = _lambw_N1<ArgumentType,CoeffType>(z);
			}
			else
			{
				if(abs(z+exp((CoeffType)(-1.)))<0.3 && z_imag>=0)
				{
					w = _lambw_sing<ArgumentType,CoeffType,IndexType>(z,k);
				}
				else
				{
					w = _lambw_asy<ArgumentType,CoeffType,IndexType>(z,-1);
				}
			}
		}
		else
		{
			if(k==1)
			{
				if(abs(z+(CoeffType)0.04)<=0.14 && z_imag<0)
				{
					w = _lambw_N1<ArgumentType,CoeffType>(z);
				}
				else
				{
					if(abs(z+exp((CoeffType)(-1.)))<0.3 && z_imag<0)
					{
						w = _lambw_sing<ArgumentType,CoeffType,IndexType>(z,k);
					}
					else
					{
						w = _lambw_asy<ArgumentType,CoeffType,IndexType>(z,1);
					}
				}
			}
			else
			{
				if(k!=0)
				{
					w = _lambw_asy<ArgumentType,CoeffType,IndexType>(z,k);
				}
				else
				{
					if(abs(z-(CoeffType)0.02)<0.28)
					{
						w = _lambw_laur<ArgumentType,CoeffType>(z);
					}
					else
					{
						if(abs(z+exp((CoeffType)(-1.)))<0.3)
						{
							w = _lambw_sing<ArgumentType,CoeffType,IndexType>(z,k);
						}
						else
						{
							if(abs(z-(CoeffType)8.)<13)
							{
								w = _lambw_mid<ArgumentType,CoeffType>(z);
							}
							else
							{
								w = _lambw_asy<ArgumentType,CoeffType,IndexType>(z,0);
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
			boost::array<ArgumentType,4> approxes = {_lambw_laur<ArgumentType,CoeffType>(z),
												_lambw_sing<ArgumentType,CoeffType,IndexType>(z,0),
												_lambw_mid<ArgumentType,CoeffType>(z),
												_lambw_asy<ArgumentType,CoeffType,IndexType>(z,0)};
			w = _lambw_best_approx(z,approxes);
		}
		else if(k == -1)
		{
			if(_lambw_complex_imag(z) >= 0)
			{
				boost::array<ArgumentType,3> approxes = {_lambw_sing<ArgumentType,CoeffType,IndexType>(z,-1),
													_lambw_N1<ArgumentType,CoeffType>(z),
													_lambw_asy<ArgumentType,CoeffType,IndexType>(z,-1)};
				w = _lambw_best_approx(z,approxes);
			}
			else
			{
				w = _lambw_asy<ArgumentType,CoeffType,IndexType>(z,-1);
			}
		}
		else if(k == 1)
		{
			if(_lambw_complex_imag(z) < 0)
			{
				boost::array<ArgumentType,3> approxes = {_lambw_sing<ArgumentType,CoeffType,IndexType>(z,1),
													_lambw_N1<ArgumentType,CoeffType>(z),
													_lambw_asy<ArgumentType,CoeffType,IndexType>(z,1)};
				w = _lambw_best_approx(z,approxes);
			}
			else
			{
				w = _lambw_asy<ArgumentType,CoeffType,IndexType>(z,1);
			}
		}
		else
		{
			w = _lambw_asy<ArgumentType,CoeffType,IndexType>(z,k);
		}
		*/

		if( (boost::math::isnan(_lambw_complex_real(w))
			||boost::math::isnan(_lambw_complex_imag(w)))
			|| (_lambw_complex_real(z)==0
			&& _lambw_complex_imag(z)==0) )
		{
			return w;
		}

		CoeffType prec = _lambw_test_score<ArgumentType,CoeffType>(z,w);

		while(prec <= 4)
		{
			w = _lambw_iter_newt<ArgumentType,CoeffType,IndexType>(z,k,w);
			prec = _lambw_test_score<ArgumentType,CoeffType>(z,w);
		}

		CoeffType prec_prev;
		ArgumentType w_prev;

		do
		{
			prec_prev = prec;
			w_prev = w;
			w = _lambw_iter_frit<ArgumentType,CoeffType,IndexType>(z,k,w);
			prec = _lambw_test_score<ArgumentType,CoeffType>(z,w);
		}
		while(prec > prec_prev);

		return w_prev;
	}

	template<typename ArgumentType>
	ArgumentType lambert_w(ArgumentType z)
	{
		if( (z == std::numeric_limits<ArgumentType>::infinity()) && (z > 0) )
		{
			return std::numeric_limits<ArgumentType>::infinity();
		}
		else
		{
			return _lambw_raw<ArgumentType,ArgumentType,int>(z, 0);
		}
	}

	template<typename ArgumentType, typename IndexType>
	ArgumentType lambert_w(ArgumentType z, IndexType k)
	{
		if( (k == 0) || (k == -1) )
		{
			if( (k == -1) )
			{
				if(z == 0)
				{
					return -std::numeric_limits<ArgumentType>::infinity();
				}
				else if( z < -std::exp(-1.f) || z > 0 )
				{
					return std::numeric_limits<ArgumentType>::quiet_NaN();
				}
				else
				{
					return _lambw_raw<ArgumentType,ArgumentType,IndexType>(z, k);
				}
			}
			else if( (k == 0) && (z == std::numeric_limits<ArgumentType>::infinity()) && (z > 0) )
			{
				return std::numeric_limits<ArgumentType>::infinity();
			}
			else
			{
				return _lambw_raw<ArgumentType,ArgumentType,IndexType>(z, k);
			}
		}
		else
		{
			return std::numeric_limits<ArgumentType>::quiet_NaN();
		}
	}

	template<typename ArgumentType>
	std::complex<ArgumentType> lambert_w(const std::complex<ArgumentType> &z)
	{
		return _lambw_raw<std::complex<ArgumentType>,ArgumentType,int>(z, 0);
	}

	template<typename ArgumentType, typename IndexType>
	std::complex<ArgumentType> lambert_w(const std::complex<ArgumentType> &z, IndexType k)
	{
		return _lambw_raw<std::complex<ArgumentType>,ArgumentType,IndexType>(z, k);
	}

	/*
	template<typename ArgumentType>
	ArgumentType wright_omega(const ArgumentType &z)
	{
		return lambert_w(std::exp(z),(long long)boost::math::round(_lambw_complex_imag(z)/boost::math::constants::two_pi<BOOST_TYPEOF(std::abs(z))>()));
	}
	*/

	template<typename ArgumentType, typename IndexType>
	ArgumentType alt_lambert_w(const ArgumentType &z, IndexType k)
	{
		if(_lambw_complex_imag(z)<0)
		{
			return lambert_w(z,-k);
		}
		else
		{
			return lambert_w(z,k);
		}
	}
} // namespace math
} // namespace boost

#endif // BOOST_MATH_LAMBW_HPP_INCLUDED
