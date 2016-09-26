
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
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/array.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/math/constants/constants.hpp>
//#include <boost/math/special_functions/factorials.hpp>

namespace boost
{
namespace math
{
	const std::size_t _lambw_series_L = 10;
	const std::size_t _lambw_asy_series_L = 4;
	const std::size_t _lambw_mid_series_L = 8;
	const std::size_t _lambw_sing_series_L = 10;

	template<typename ConstantType>
	const ConstantType _lambw_imaginary_unit()
	{
		return std::numeric_limits<ConstantType>::quiet_NaN();
	}

	template<>
	const std::complex<float> _lambw_imaginary_unit<std::complex<float> >()
	{
		return std::complex<float>((float)0.,(float)1.);
	}

	template<>
	const std::complex<double> _lambw_imaginary_unit<std::complex<double> >()
	{
		return std::complex<double>((double)0.,(double)1.);
	}

	template<>
	const std::complex<long double> _lambw_imaginary_unit<std::complex<long double> >()
	{
		return std::complex<long double>((long double)0.,(long double)1.);
	}

	template<typename ArgumentType, typename ConstantType, std::size_t L>
	ArgumentType _lambw_series_sum(const ArgumentType &z, const boost::array<ConstantType,L> &c)
	{
		ArgumentType sum = 1;

		for(std::size_t k = 0; k < L; ++k)
		{
			sum = ((ConstantType)1.)+c[L-k-1]*z*sum;
		}

		return sum;
	}

	template<typename ConstantType, std::size_t L>
	boost::array<ConstantType,L> _lambw_coeff_array(ConstantType (&coeffs)(std::size_t) )
	{
		boost::array<ConstantType,L> ans;
		for(std::size_t k = 0; k < L; ++k)
		{
			ans[k] = coeffs(k);
		}
		return ans;
	}

	///constexpr-s in these require C++14
    template<typename ArgumentType>
	/*BOOST_CONSTEXPR*/ ArgumentType _lambw_linstrips(const ArgumentType &z)
	{
		if(boost::is_arithmetic<ArgumentType>::value)
		{
			return z;
		}
		else
		{
			return z-boost::math::constants::two_pi<BOOST_TYPEOF(std::real(z))>()
			*_lambw_imaginary_unit<ArgumentType>()*boost::math::round(std::imag(z)/boost::math::constants::two_pi<BOOST_TYPEOF(std::real(z))>());
		}
	}

	template<typename ArgumentType>
	ArgumentType _lambw_test_s(const ArgumentType &z, const ArgumentType &w)
	{
		if(boost::is_arithmetic<ArgumentType>::value)
		{
			return _lambw_linstrips(std::log(std::abs(w))+w-std::log(std::abs(z)));
		}
		else
		{
			return _lambw_linstrips(std::log(w)+w-std::log(z));
		}
	}

	/*
	template<typename ArgumentType>
	ArgumentType _lambw_test(const ArgumentType &z, const ArgumentType &w)
	{
		return expm1(_lambw_test_s(z,w));
	}
	*/

	template<typename ArgumentType>
	int _lambw_test_score(const ArgumentType &z, const ArgumentType &w)
	{
		BOOST_TYPEOF(std::real(z)) score = -std::log10(std::abs(_lambw_test_s(z,w)));

		if(score == std::numeric_limits<BOOST_TYPEOF(std::real(z))>::infinity())
		{
			return std::numeric_limits<int>::max();
		}
		else
		{
			return (int)boost::math::round(score);
		}
	}

	template<typename ConstantType>
    /*BOOST_CONSTEXPR*/ ConstantType _lambw_tay_coeffs(std::size_t k)
    {
    	ConstantType n = (ConstantType)(k+1);
    	return (-std::pow(((n+1)/n),n-1));
    }

    template<typename ConstantType>
    /*BOOST_CONSTEXPR*/ ConstantType _lambw_laur_coeffs(std::size_t k)
    {
    	ConstantType n = (ConstantType)k;
    	switch(k)
    	{
		case 0:
			return 1.;
			break;
		case 1:
			return -0.5;
			break;
		default:
			return std::pow(n/(n-(ConstantType)1.),n)*((ConstantType)1.-n)/(n+(ConstantType)1.);
			break;
    	}
    }

	template<typename ConstantType>
	boost::array<ConstantType,_lambw_series_L> _lw_tc();

	const boost::array<float,_lambw_series_L> _lw_tc_float = _lambw_coeff_array<float,_lambw_series_L>(_lambw_tay_coeffs<float>);

	template<>
	boost::array<float,_lambw_series_L> _lw_tc<float>()
	{
		return _lw_tc_float;
	}

	const boost::array<double,_lambw_series_L> _lw_tc_double = _lambw_coeff_array<double,_lambw_series_L>(_lambw_tay_coeffs<double>);

	template<>
	boost::array<double,_lambw_series_L> _lw_tc<double>()
	{
		return _lw_tc_double;
	}

	const boost::array<long double,_lambw_series_L> _lw_tc_long_double = _lambw_coeff_array<long double,_lambw_series_L>(_lambw_tay_coeffs<long double>);

	template<>
	boost::array<long double,_lambw_series_L> _lw_tc<long double>()
	{
		return _lw_tc_long_double;
	}

	template<typename ConstantType>
	boost::array<ConstantType,_lambw_series_L> _lw_lc();

	const boost::array<float,_lambw_series_L> _lw_lc_float = _lambw_coeff_array<float,_lambw_series_L>(_lambw_laur_coeffs<float>);

	template<>
	boost::array<float,_lambw_series_L> _lw_lc<float>()
	{
		return _lw_lc_float;
	}

	const boost::array<double,_lambw_series_L> _lw_lc_double = _lambw_coeff_array<double,_lambw_series_L>(_lambw_laur_coeffs<double>);

	template<>
	boost::array<double,_lambw_series_L> _lw_lc<double>()
	{
		return _lw_lc_double;
	}

	const boost::array<long double,_lambw_series_L> _lw_lc_long_double = _lambw_coeff_array<long double,_lambw_series_L>(_lambw_laur_coeffs<long double>);

	template<>
	boost::array<long double,_lambw_series_L> _lw_lc<long double>()
	{
		return _lw_lc_long_double;
	}

	template<typename ArgumentType,typename ConstantType>
	ArgumentType _lambw_tay(const ArgumentType &z)
	{
		return z*_lambw_series_sum(z,_lw_tc<ConstantType>());
	}

	template<typename ArgumentType,typename ConstantType>
	ArgumentType _lambw_laur(const ArgumentType &z)
	{
		return z/_lambw_series_sum(z,_lw_lc<ConstantType>());
	}

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

	template<typename ArgumentType, typename IndexType>
	ArgumentType _log_k(const ArgumentType &z, IndexType k = 0)
	{
		if(boost::is_arithmetic<ArgumentType>::value)
		{
			return (k==0)?(std::log(z)):(std::numeric_limits<ArgumentType>::quiet_NaN());
		}
		else
		{
			return std::log(z) + (boost::math::constants::two_pi<BOOST_TYPEOF(std::real(z))>()*k)*_lambw_imaginary_unit<ArgumentType>();
		}
	}

	template<typename ArgumentType, typename ConstantType>
	ArgumentType _lambw_asy_body(const ArgumentType &L1, const ArgumentType &L2)
	{
		ArgumentType ans = L1-L2;
		/*
		for(std::size_t l = 0; l < _lambw_asy_series_L; ++l)
		{
			for(std::size_t m = 1; m <= l+1; ++m)
			{
				ans+=std::pow(-1.,l%2)*_stirling(l+m,l+1)/boost::math::factorial(m)*std::pow(L1,-(ConstantType)l-(ConstantType)m)*std::pow(L2,(ConstantType)m);
			}
		}
		*/
		ans += L2/L1 + L2*(((ConstantType)(-2.))+L2)/(ConstantType)2./L1/L1 + L2*((ConstantType)6.-(ConstantType)9.*L2+(ConstantType)2.*L2*L2)/(ConstantType)6./L1/L1/L1
			+ L2*((ConstantType)(-12.)+(ConstantType)36.*L2-(ConstantType)22.*L2*L2+(ConstantType)3.*L2*L2*L2)/(ConstantType)12./L1/L1/L1/L1;


		return ans;
	}

	template<typename ArgumentType, typename ConstantType, typename IndexType>
	ArgumentType _lambw_asy(const ArgumentType &z, IndexType k = 0)
	{
		ArgumentType L1 = _log_k(z,k);
		ArgumentType L2 = std::log(L1);

		return _lambw_asy_body<ArgumentType,ConstantType>(L1,L2);
	}

	template<typename ArgumentType, typename ConstantType>
	ArgumentType _lambw_N1(const ArgumentType &z)
	{
		ArgumentType L1 = std::log(-z);
		ArgumentType L2 = std::log(-L1);

		return _lambw_asy_body<ArgumentType,ConstantType>(L1,L2);
	}

	template<typename ConstantType>
	boost::array<ConstantType,_lambw_series_L> _lw_midc();

	const boost::array<float,_lambw_series_L> _lw_midc_float = {(1.f/2.f), (1.f/8.f),-(1.f/12.f),(1.f/16.f),-(13.f/20.f),-(47.f/312.f),(73.f/1316.f),-(2447.f/2336.f)};

	template<>
	boost::array<float,_lambw_series_L> _lw_midc<float>()
	{
		return _lw_midc_float;
	}

	const boost::array<double,_lambw_series_L> _lw_midc_double = {(1./2.), (1./8.),-(1./12.),(1./16.),-(13./20.),-(47./312.),(73./1316.),-(2447./2336.)};

	template<>
	boost::array<double,_lambw_series_L> _lw_midc<double>()
	{
		return _lw_midc_double;
	}

	const boost::array<long double,_lambw_series_L> _lw_midc_long_double = {(1.l/2.l), (1.l/8.l),-(1.l/12.l),(1.l/16.l),-(13.l/20.l),-(47.l/312.l),(73.l/1316.l),-(2447.l/2336.l)};

	template<>
	boost::array<long double,_lambw_series_L> _lw_midc<long double>()
	{
		return _lw_midc_long_double;
	}

	template<typename ArgumentType, typename ConstantType>
	ArgumentType _lambw_mid(const ArgumentType &z)
	{
		ArgumentType x = std::log(z)-(ConstantType)1.;
		return _lambw_series_sum(x,_lw_midc<ConstantType>());
	}

	template<typename ConstantType>
	/*BOOST_CONSTEXPR*/ ConstantType _lambw_sing_coeffs_mu(std::size_t k);

	template<typename ConstantType>
	/*BOOST_CONSTEXPR*/ ConstantType _lambw_sing_coeffs_alpha(std::size_t k)
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
			ConstantType sum = 0;
			for(std::size_t j = 2; j <= (k-1); ++j)
			{
				sum+=_lambw_sing_coeffs_mu<ConstantType>(j)*_lambw_sing_coeffs_mu<ConstantType>(k+1-j);
			}
			return sum;
		}
	}

	template<typename ConstantType>
	/*BOOST_CONSTEXPR*/ ConstantType _lambw_sing_coeffs_mu(std::size_t k)
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
			return (_lambw_sing_coeffs_mu<ConstantType>(k-2)/2
					+_lambw_sing_coeffs_alpha<ConstantType>(k-2)/4)*(k-1)/(k+1)
					-_lambw_sing_coeffs_alpha<ConstantType>(k)/2
					-_lambw_sing_coeffs_mu<ConstantType>(k-1)/(k+1);
		}
	}

	template<typename ConstantType>
	/*BOOST_CONSTEXPR*/ ConstantType _lambw_sing_coeffs(std::size_t k)
	{
		ConstantType ans = _lambw_sing_coeffs_mu<ConstantType>(k+1)/_lambw_sing_coeffs_mu<ConstantType>(k);
		return ans;
	}

	template<typename ConstantType>
	boost::array<ConstantType,_lambw_series_L> _lw_singc();

	const boost::array<float,_lambw_series_L> _lw_singc_float = _lambw_coeff_array<float,_lambw_series_L>(_lambw_sing_coeffs<float>);

	template<>
	boost::array<float,_lambw_series_L> _lw_singc<float>()
	{
		return _lw_singc_float;
	}

	const boost::array<double,_lambw_series_L> _lw_singc_double = _lambw_coeff_array<double,_lambw_series_L>(_lambw_sing_coeffs<double>);

	template<>
	boost::array<double,_lambw_series_L> _lw_singc<double>()
	{
		return _lw_singc_double;
	}

	const boost::array<long double,_lambw_series_L> _lw_singc_long_double = _lambw_coeff_array<long double,_lambw_series_L>(_lambw_sing_coeffs<long double>);

	template<>
	boost::array<long double,_lambw_series_L> _lw_singc<long double>()
	{
		return _lw_singc_long_double;
	}

	template<typename ArgumentType, typename ConstantType, typename IndexType>
	ArgumentType _lambw_sing(const ArgumentType &z, IndexType k)
	{
		ArgumentType p;
		if(std::real(z) >= 0)
		{
			p = ((ConstantType)std::pow(-1.,k%2))*std::sqrt((ConstantType)2.*(std::exp((ConstantType)1.+std::log(z))+(ConstantType)1.));
		}
		else
		{
			p = ((ConstantType)std::pow(-1.,k%2))*std::sqrt((ConstantType)2.*(-std::exp((ConstantType)1.+std::log(-z))+(ConstantType)1.));
		}
		return -_lambw_series_sum(p,_lw_singc<ConstantType>());
	}

	template<typename ArgumentType, typename ConstantType, typename IndexType>
	ArgumentType _lambw_iter_frit(const ArgumentType &x, IndexType k, const ArgumentType &Wn)
	{
		if((k==0||k==1||k==-1)&&(x==(ConstantType)0.||x==-(ConstantType)std::exp(-1.)))
		{
			return Wn;
		}
		else
		{
			ArgumentType z_n;

			if(k==0)
			{
				if(std::imag(x)==0 &&(std::real(x)<=0)&&(std::real(x)>=-(ConstantType)std::exp(-1.l)) )
				{
					z_n = std::log(std::abs(x))-std::log(std::abs(Wn))-Wn;
				}
				else
				{
					z_n = std::log(x)-std::log(Wn)-Wn;
				}
			}
			else
			{
				if((k==-1)&&(std::imag(x)==0)&&(std::real(x)<=0)&&(std::real(x)>=-(ConstantType)std::exp(-1.l)))
				{
					z_n = std::log(std::abs(x))-std::log(std::abs(Wn))-Wn;
				}
				else
				{
					z_n = _log_k(x,k)+std::log((ConstantType)1./Wn)-Wn;
				}
			}

			ArgumentType q_n=(ConstantType)2.*((ConstantType)1.+Wn)*((ConstantType)1.+Wn+(ConstantType)2./(ConstantType)3.*z_n);

			ArgumentType e_n=z_n/((ConstantType)1.+Wn)*(q_n-z_n)/(q_n-(ConstantType)2.*z_n);

			return Wn*((ConstantType)1.+e_n);
		}
	}

	template<typename ArgumentType, typename ConstantType, typename IndexType>
	ArgumentType _lambw_iter_newt_nolog(const ArgumentType &x, const ArgumentType &fn)
	{
		return fn*((ConstantType)1.-(ConstantType)1./((ConstantType)1.+fn))+x/((ConstantType)1.+fn)*std::exp(-fn);
	}

	template<typename ArgumentType, typename ConstantType, typename IndexType>
	ArgumentType _lambw_iter_newt_wlog(const ArgumentType &x, const ArgumentType &fn)
	{
		if(std::real(x)<0)
		{
			if(std::real(fn)<0)
			{
				return fn-(ConstantType)1./((ConstantType)1.+(ConstantType)1./fn)+std::exp(std::log(-x)-std::log(-(ConstantType)1.-fn)-(fn));
			}
			else
			{
				return fn-(ConstantType)1./((ConstantType)1.+(ConstantType)1./fn)-std::exp(std::log(-x)-std::log((ConstantType)1.+fn)-(fn));
			}
		}
		else
		{
			if(std::real(fn)<0)
			{
				return fn-(ConstantType)1./((ConstantType)1.+(ConstantType)1./fn)-std::exp(std::log(x)-std::log(-(ConstantType)1.-fn)-(fn));
			}
			else
			{
				return fn-(ConstantType)1./((ConstantType)1.+(ConstantType)1./fn)+std::exp(std::log(x)-std::log((ConstantType)1.+fn)-(fn));
			}
		}
	}

	template<typename ArgumentType, typename ConstantType, typename IndexType>
	ArgumentType _lambw_iter_newt(const ArgumentType &x, IndexType k, const ArgumentType &fn)
	{
		if(((((k==0)||(k==-1))&&(std::imag(x)>=0))||((k==1)&&(std::imag(x)<0)))
			&&(std::abs(x+(ConstantType)std::exp(-1.l))<7e-4))
		{
			return fn;
		}
		else
		{
			if(((k==0)&&(std::abs(x)<(ConstantType)std::exp(-1.l)))
				||(k==-1&&std::imag(x)>=0&&std::abs(x+(ConstantType)std::exp(-1.l))<0.3)
				||(k==1&&std::imag(x)<0&&std::abs(x+(ConstantType)std::exp(-1.l))<0.3))
			{
				return _lambw_iter_newt_nolog<ArgumentType,ConstantType,IndexType>(x,fn);
			}
			else
			{
				return _lambw_iter_newt_wlog<ArgumentType,ConstantType,IndexType>(x,fn);
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

	template<typename ArgumentType, typename ConstantType, typename IndexType, int Precision>
	ArgumentType _lambw_raw(const ArgumentType &z, IndexType k = 0)
	{
		ArgumentType w;

		ConstantType z_imag = std::imag(z);

		if(k==-1)
		{
			if(std::abs(z+(ConstantType)0.04)<=0.14 && z_imag>=0)
			{
				w = _lambw_N1<ArgumentType,ConstantType>(z);
			}
			else
			{
				if(abs(z+(ConstantType)std::exp(-1.l))<0.3 && z_imag>=0)
				{
					w = _lambw_sing<ArgumentType,ConstantType,IndexType>(z,k);
				}
				else
				{
					w = _lambw_asy<ArgumentType,ConstantType,IndexType>(z,-1);
				}
			}
		}
		else
		{
			if(k==1)
			{
				if(std::abs(z+(ConstantType)0.04)<=0.14 && z_imag<0)
				{
					w = _lambw_N1<ArgumentType,ConstantType>(z);
				}
				else
				{
					if(std::abs(z+(ConstantType)std::exp(-1.l))<0.3 && z_imag<0)
					{
						w = _lambw_sing<ArgumentType,ConstantType,IndexType>(z,k);
					}
					else
					{
						w = _lambw_asy<ArgumentType,ConstantType,IndexType>(z,1);
					}
				}
			}
			else
			{
				if(k!=0)
				{
					w = _lambw_asy<ArgumentType,ConstantType,IndexType>(z,k);
				}
				else
				{
					if(std::abs(z-(ConstantType)0.02)<0.28)
					{
						w = _lambw_laur<ArgumentType,ConstantType>(z);
					}
					else
					{
						if(std::abs(z+(ConstantType)std::exp(-1.))<0.3)
						{
							w = _lambw_sing<ArgumentType,ConstantType,IndexType>(z,k);
						}
						else
						{
							if(std::abs(z-(ConstantType)8.)<13)
							{
								w = _lambw_mid<ArgumentType,ConstantType>(z);
							}
							else
							{
								w = _lambw_asy<ArgumentType,ConstantType,IndexType>(z,0);
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
			boost::array<ArgumentType,4> approxes = {_lambw_laur<ArgumentType,ConstantType>(z),
												_lambw_sing<ArgumentType,ConstantType,IndexType>(z,0),
												_lambw_mid<ArgumentType,ConstantType>(z),
												_lambw_asy<ArgumentType,ConstantType,IndexType>(z,0)};
			w = _lambw_best_approx(z,approxes);
		}
		else if(k == -1)
		{
			if(std::imag(z) >= 0)
			{
				boost::array<ArgumentType,3> approxes = {_lambw_sing<ArgumentType,ConstantType,IndexType>(z,-1),
													_lambw_N1<ArgumentType,ConstantType>(z),
													_lambw_asy<ArgumentType,ConstantType,IndexType>(z,-1)};
				w = _lambw_best_approx(z,approxes);
			}
			else
			{
				w = _lambw_asy<ArgumentType,ConstantType,IndexType>(z,-1);
			}
		}
		else if(k == 1)
		{
			if(std::imag(z) < 0)
			{
				boost::array<ArgumentType,3> approxes = {_lambw_sing<ArgumentType,ConstantType,IndexType>(z,1),
													_lambw_N1<ArgumentType,ConstantType>(z),
													_lambw_asy<ArgumentType,ConstantType,IndexType>(z,1)};
				w = _lambw_best_approx(z,approxes);
			}
			else
			{
				w = _lambw_asy<ArgumentType,ConstantType,IndexType>(z,1);
			}
		}
		else
		{
			w = _lambw_asy<ArgumentType,ConstantType,IndexType>(z,k);
		}
		*/

		if( (boost::math::isnan(std::real(w))||boost::math::isnan(std::imag(w))) || (std::real(z)==0 && std::imag(z)==0) )
		{
			return w;
		}

		int prec = _lambw_test_score(z,w);

		if(prec >= Precision)
		{
			return w;
		}

		while(prec <= Precision/4)
		{
			w = _lambw_iter_newt<ArgumentType,ConstantType,IndexType>(z,k,w);
			prec = _lambw_test_score(z,w);
		}

		return _lambw_iter_frit<ArgumentType,ConstantType,IndexType>(z,k,w);
	}

	float lambert_w(float z, int k = 0)
	{
		if( (k == 0) || (k == -1) )
		{
			if( (k == -1) )
			{
				if(z == 0)
				{
					return -std::numeric_limits<float>::infinity();
				}
				else if( z < -std::exp(-1.f) || z > 0 )
				{
					return std::numeric_limits<float>::quiet_NaN();
				}
				else
				{
					return _lambw_raw<float,float,int,8>(z, k);
				}
			}
			else if( (k == 0) && (std::isinf(z)) && (z > 0) )
			{
				return std::numeric_limits<float>::infinity();
			}
			else
			{
				return _lambw_raw<float,float,int,8>(z, k);
			}
		}
		else
		{
			return std::numeric_limits<float>::quiet_NaN();
		}
	}

	double lambert_w(double z, int k = 0)
	{
		if( (k == 0) || (k == -1) )
		{
			if( (k == -1) )
			{
				if(z == 0)
				{
					return -std::numeric_limits<double>::infinity();
				}
				else if( z < -std::exp(-1.) || z > 0 )
				{
					return std::numeric_limits<double>::quiet_NaN();
				}
				else
				{
					return _lambw_raw<double,double,int,16>(z, k);
				}
			}
			else if( (k == 0) && (std::isinf(z)) && (z > 0) )
			{
				return std::numeric_limits<double>::infinity();
			}
			else
			{
				return _lambw_raw<double,double,int,16>(z, k);
			}
		}
		else
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
	}

	long double lambert_w(long double z, int k = 0)
	{
		if( (k == 0) || (k == -1) )
		{
			if( (k == -1) )
			{
				if(z == 0)
				{
					return -std::numeric_limits<long double>::infinity();
				}
				else if( z < -std::exp(-1.l) || z > 0 )
				{
					return std::numeric_limits<long double>::quiet_NaN();
				}
				else
				{
					return _lambw_raw<long double,long double,int,24>(z, k);
				}
			}
			else if( (k == 0) && (std::isinf(z)) && (z > 0) )
			{
				return std::numeric_limits<long double>::infinity();
			}
			else
			{
				return _lambw_raw<long double,long double,int,24>(z, k);
			}
		}
		else
		{
			return std::numeric_limits<long double>::quiet_NaN();
		}
	}

	std::complex<float> lambert_w(const std::complex<float> &z)
	{
		return _lambw_raw<std::complex<float>,float,int,8>(z, 0);
	}

	std::complex<double> lambert_w(const std::complex<double> &z)
	{
		return _lambw_raw<std::complex<double>,double,int,16>(z, 0);
	}

	std::complex<long double> lambert_w(const std::complex<long double> &z)
	{
		return _lambw_raw<std::complex<long double>,long double,int,24>(z, 0);
	}

	template<typename IndexType>
	std::complex<float> lambert_w(const std::complex<float> &z, IndexType k)
	{
		return _lambw_raw<std::complex<float>,float,IndexType,8>(z, k);
	}

	template<typename IndexType>
	std::complex<double> lambert_w(const std::complex<double> &z, IndexType k)
	{
		return _lambw_raw<std::complex<double>,double,IndexType,16>(z, k);
	}

	template<typename IndexType>
	std::complex<long double> lambert_w(const std::complex<long double> &z, IndexType k)
	{
		return _lambw_raw<std::complex<long double>,long double,IndexType,24>(z, k);
	}

	/*
	template<typename ArgumentType>
	ArgumentType wright_omega(const ArgumentType &z)
	{
		return lambert_w(std::exp(z),(long long)boost::math::round(std::imag(z)/boost::math::constants::two_pi<BOOST_TYPEOF(std::imag(z))>()));
	}
	*/

	template<typename ArgumentType, typename IndexType>
	ArgumentType alt_lambert_w(const ArgumentType &z, IndexType k)
	{
		if(std::imag(z)<0)
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
