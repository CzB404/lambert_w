#ifndef LAMBW_HPP_INCLUDED
#define LAMBW_HPP_INCLUDED

#include <array>
#include <cmath>
#include <complex>
//#include <functional>
#include <limits>
#include <type_traits>
#include <boost/math/constants/constants.hpp>

namespace czb
{
	const std::size_t _series_L = 10;
	const std::size_t _asy_series_L = 4;
	const std::size_t _mid_series_L = 8;
	const std::size_t _sing_series_L = 10;
	const std::size_t _expm1_series_L = 64;

	template<typename Tc>
	Tc _imaginary_unit = std::numeric_limits<Tc>::quiet_NaN();

	template<>
	std::complex<float> _imaginary_unit<std::complex<float>> = std::complex<float>((float)0.,(float)1.);

	template<>
	std::complex<double> _imaginary_unit<std::complex<double>> = std::complex<double>((double)0.,(double)1.);

	template<>
	std::complex<long double> _imaginary_unit<std::complex<long double>> = std::complex<long double>((long double)0.,(long double)1.);

	template<typename Tz, typename Tc, std::size_t L>
	Tz _series_sum(const Tz &z, const std::array<Tc,L> c)
	{
		Tz sum = 1;

		for(std::size_t k = 0; k < L; ++k)
		{
			sum = ((Tc)1.)+c[L-k-1]*z*sum;
			//std::cout << sum << std::endl;
		}

		return sum;
	}

	template<typename Tc, std::size_t L>
	std::array<Tc,L> _coeff_array(Tc (&coeffs)(std::size_t) )
	{
		std::array<Tc,L> ans;
		for(std::size_t k = 0; k < L; ++k)
		{
			ans[k] = coeffs(k);
		}
		return ans;
	}

	template<typename Tc>
    constexpr Tc _expm1_coeffs(std::size_t k)
    {
    	Tc n = ((Tc)1.)/(k+2);
    	return n;
    }

    template<typename Tc>
	std::array<Tc,_expm1_series_L> _lw_expm1c = czb::_coeff_array<Tc,_expm1_series_L>(czb::_expm1_coeffs<Tc>);

	template<typename Tz>
	constexpr Tz _linstrips(const Tz &z)
	{
		if(std::is_arithmetic<Tz>::value)
		{
			return z;
		}
		else
		{
			return z-boost::math::constants::two_pi<decltype(std::real(z))>()*_imaginary_unit<Tz>*std::round(std::imag(z)/boost::math::constants::two_pi<decltype(std::real(z))>());
		}
	}

	template<typename Tz>
	Tz expm1(const Tz &z)
	{
		if(std::abs(std::real(z)) > 1.)
		{
			return std::exp(z)-(decltype(std::real(z)))1.;
		}
		else
		{
			Tz z0 = _linstrips(z);
			return z0*czb::_series_sum(z0,_lw_expm1c<decltype(std::real(z))>);
		}
	}

	template<typename Tz>
	Tz _lambw_test_s(const Tz &z, const Tz &w)
	{
		if(std::is_arithmetic<Tz>::value)
		{
			return _linstrips(std::log(std::abs(w))+w-std::log(std::abs(z)));
		}
		else
		{
			return _linstrips(std::log(w)+w-std::log(z));
		}
	}

	template<typename Tz>
	Tz _lambw_test(const Tz &z, const Tz &w)
	{
		return expm1(_lambw_test_s(z,w));
	}

	template<typename Tz>
	int _lambw_test_score(const Tz &z, const Tz &w)
	{
		auto score = -std::log10(std::abs(_lambw_test_s(z,w)));

		if(score == std::numeric_limits<decltype(std::real(z))>::infinity())
		{
			return std::numeric_limits<int>::max();
		}
		else
		{
			return std::round(score);
		}
	}

	/*
	template<typename Tz>
	class LambW_Test
	{
	public:
		LambW_Test(const std::function<Tz(Tz)> &&func,
					const std::function<Tz(Tz,Tz)> &&test_func = _lambw_test_s<Tz>)
					: _func(func), _test_func(test_func)
		{

		}

		Tz operator()(const Tz &z)
		{
			Tz w = _func(z);
			return _test_func(z,w);
		}

		Tz f(const Tz &z)
		{
			return _func(z);
		}
	private:
		std::function<Tz(Tz)> _func;
		std::function<Tz(Tz,Tz)> _test_func;
	};
	*/

    template<typename Tc>
    constexpr Tc _lambw_tay_coeffs(std::size_t k)
    {
    	Tc n = k+1;
    	return (-std::pow(((n+1)/n),n-1));
    }

    template<typename Tc>
    constexpr Tc _lambw_laur_coeffs(std::size_t k)
    {
    	Tc n = k;
    	switch(k)
    	{
		case 0:
			return 1.;
			break;
		case 1:
			return -0.5;
			break;
		default:
			return std::pow(n/(n-1.),n)*(1.-n)/(n+1.);
			break;
    	}
    }

	template<typename Tc>
	std::array<Tc,_series_L> _lw_tc = czb::_coeff_array<Tc,_series_L>(czb::_lambw_tay_coeffs<Tc>);

	template<typename Tc>
	std::array<Tc,_series_L> _lw_lc = czb::_coeff_array<Tc,_series_L>(czb::_lambw_laur_coeffs<Tc>);

	template<typename Tz = double,typename Tc = double>
	Tz _lambw_tay(const Tz &z)
	{
		return z*czb::_series_sum(z,_lw_tc<Tc>);
	}

	template<typename Tz = double,typename Tc = double>
	Tz _lambw_laur(const Tz &z)
	{
		return z/czb::_series_sum(z,_lw_lc<Tc>);
	}

	template<typename IntegralType>
	constexpr IntegralType _stirling(IntegralType n, IntegralType k)
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

	template<typename Tz, typename Ti>
	Tz _log_k(const Tz &z, Ti k = 0)
	{
		if(std::is_arithmetic<Tz>::value)
		{
			return (k==0)?(std::log(z)):(std::numeric_limits<Tz>::quiet_NaN());
		}
		else
		{
			return std::log(z) + (boost::math::constants::two_pi<decltype(std::real(z))>()*k)*_imaginary_unit<Tz>;
		}
	}

	template<typename IntegralType>
	constexpr IntegralType fac(IntegralType n)
	{
		IntegralType ans = 1;
		for(IntegralType k = 2; k <= n; ++k)
		{
			ans*=k;
		}
		return ans;
	}

	template<typename Tz, typename Tc>
	Tz _lambw_asy_body(const Tz &L1, const Tz &L2)
	{
		Tz ans = L1-L2;
		/*
		for(std::size_t l = 0; l < _asy_series_L; ++l)
		{
			for(std::size_t m = 1; m <= l+1; ++m)
			{
				ans+=std::pow(-1.,l%2)*_stirling(l+m,l+1)/fac(m)*std::pow(L1,-(Tc)l-(Tc)m)*std::pow(L2,(Tc)m);
			}
		}
		*/
		ans += L2/L1 + L2*(((Tc)(-2.))+L2)/(Tc)2./L1/L1 + L2*((Tc)6.-(Tc)9.*L2+(Tc)2.*L2*L2)/(Tc)6./L1/L1/L1
			+ L2*((Tc)(-12.)+(Tc)36.*L2-(Tc)22.*L2*L2+(Tc)3.*L2*L2*L2)/(Tc)12./L1/L1/L1/L1;


		return ans;
	}

	template<typename Tz = double, typename Tc = double, typename Ti = int>
	Tz _lambw_asy(const Tz &z, Ti k = 0)
	{
		Tz L1 = _log_k(z,k);
		Tz L2 = std::log(L1);

		return _lambw_asy_body<Tz,Tc>(L1,L2);
	}

	template<typename Tz = double, typename Tc = double>
	Tz _lambw_N1(const Tz &z)
	{
		Tz L1 = std::log(-z);
		Tz L2 = std::log(-L1);

		return _lambw_asy_body<Tz,Tc>(L1,L2);
	}
/*
	template<typename Ti>
	std::complex<double> _log_k<std::complex<double>,Ti>(const std::complex<double> &z, Ti k)
	{
		return std::log(z) + ((double)(2*k*boost::math::constants::pi<double>))*(std::complex<double>(0,1));
	}
*/

	template<typename Tc>
	std::array<Tc,_mid_series_L> _lw_midc = {(1./2.), (1./8.),-(1./12.),(1./16.),-(13./20.),-(47./312.),(73./1316.),-(2447./2336.)};

	template<typename Tz, typename Tc = double>
	Tz _lambw_mid(const Tz &z)
	{
		Tz x = std::log(z)-(Tc)1.;
		return _series_sum(x,_lw_midc<Tc>);
	}

	template<typename Tc>
	constexpr Tc _lambw_sing_coeffs_mu(std::size_t k);

	template<typename Tc>
	constexpr Tc _lambw_sing_coeffs_alpha(std::size_t k)
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
			Tc sum = 0;
			for(std::size_t j = 2; j <= (k-1); ++j)
			{
				sum+=_lambw_sing_coeffs_mu<Tc>(j)*_lambw_sing_coeffs_mu<Tc>(k+1-j);
			}
			return sum;
		}
	}

	template<typename Tc>
	constexpr Tc _lambw_sing_coeffs_mu(std::size_t k)
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
			return (_lambw_sing_coeffs_mu<Tc>(k-2)/2+_lambw_sing_coeffs_alpha<Tc>(k-2)/4)*(k-1)/(k+1)-_lambw_sing_coeffs_alpha<Tc>(k)/2-_lambw_sing_coeffs_mu<Tc>(k-1)/(k+1);
		}
	}

	template<typename Tc>
	constexpr Tc _lambw_sing_coeffs(std::size_t k)
	{
		Tc ans = _lambw_sing_coeffs_mu<Tc>(k+1)/_lambw_sing_coeffs_mu<Tc>(k);
		return ans;
	}

	template<typename Tc>
	std::array<Tc,_sing_series_L> _lw_singc = czb::_coeff_array<Tc,_sing_series_L>(czb::_lambw_sing_coeffs<Tc>);

	template<typename Tz = double, typename Tc = double, typename Ti = int>
	Tz _lambw_sing(const Tz &z, Ti k)
	{
		Tz p;
		if(std::real(z) >= 0)
		{
			p = ((Tc)std::pow(-1.,k%2))*std::sqrt((Tc)2.*(std::exp((Tc)1.+std::log(z))+(Tc)1.));
		}
		else
		{
			p = ((Tc)std::pow(-1.,k%2))*std::sqrt((Tc)2.*(-std::exp((Tc)1.+std::log(-z))+(Tc)1.));
		}
		//std:: cout << p << std::endl;
		return -_series_sum(p,_lw_singc<Tc>);
	}

	template<typename Tz, typename Tc, typename Ti>
	Tz _iter_frit(const Tz &x, Ti k, const Tz &Wn)
	{
		if((k==0||k==1||k==-1)&&(x==(Tc)0.||x==-(Tc)std::exp(-1.)))
		{
			return Wn;
		}
		else
		{
			Tz z_n;

			if(k==0)
			{
				if(std::imag(x)==0 &&(std::real(x)<=0)&&(std::real(x)>=-(Tc)std::exp(-1.l)) )
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
				if((k==-1)&&(std::imag(x)==0)&&(std::real(x)<=0)&&(std::real(x)>=-(Tc)std::exp(-1.l)))
				{
					z_n = std::log(std::abs(x))-std::log(std::abs(Wn))-Wn;
				}
				else
				{
					z_n = _log_k(x,k)+std::log((Tc)1./Wn)-Wn;
				}
			}

			Tz q_n=(Tc)2.*((Tc)1.+Wn)*((Tc)1.+Wn+(Tc)2./(Tc)3.*z_n);

			Tz e_n=z_n/((Tc)1.+Wn)*(q_n-z_n)/(q_n-(Tc)2.*z_n);

			return Wn*((Tc)1.+e_n);
		}
	}

	template<typename Tz, typename Tc, typename Ti>
	Tz _iter_newt_nolog(const Tz &x, const Tz &fn)
	{
		return fn*((Tc)1.-(Tc)1./((Tc)1.+fn))+x/((Tc)1.+fn)*std::exp(-fn);
	}

	template<typename Tz, typename Tc, typename Ti>
	Tz _iter_newt_wlog(const Tz &x, const Tz &fn)
	{
		if(std::real(x)<0)
		{
			if(std::real(fn)<0)
			{
				return fn-(Tc)1./((Tc)1.+(Tc)1./fn)+std::exp(std::log(-x)-std::log(-(Tc)1.-fn)-(fn));
			}
			else
			{
				return fn-(Tc)1./((Tc)1.+(Tc)1./fn)-std::exp(std::log(-x)-std::log((Tc)1.+fn)-(fn));
			}
		}
		else
		{
			if(std::real(fn)<0)
			{
				return fn-(Tc)1./((Tc)1.+(Tc)1./fn)-std::exp(std::log(x)-std::log(-(Tc)1.-fn)-(fn));
			}
			else
			{
				return fn-(Tc)1./((Tc)1.+(Tc)1./fn)+std::exp(std::log(x)-std::log((Tc)1.+fn)-(fn));
			}
		}
	}

	template<typename Tz, typename Tc, typename Ti>
	Tz _iter_newt(const Tz &x, Ti k, const Tz &fn)
	{
		if(((((k==0)||(k==-1))&&(std::imag(x)>=0))||((k==1)&&(std::imag(x)<0)))&&(std::abs(x+(Tc)std::exp(-1.l))<7e-4))
		{
			return fn;
		}
		else
		{
			if(((k==0)&&(std::abs(x)<(Tc)std::exp(-1.l)))||(k==-1&&std::imag(x)>=0&&std::abs(x+(Tc)std::exp(-1.l))<0.3)||(k==1&&std::imag(x)<0&&std::abs(x+(Tc)std::exp(-1.l))<0.3))
			{
				return _iter_newt_nolog<Tz,Tc,Ti>(x,fn);
			}
			else
			{
				return _iter_newt_wlog<Tz,Tc,Ti>(x,fn);
			}
		}
	}

	template<typename Tz = double, typename Tc = double, typename Ti = int, int Precision>
	Tz _lambw_raw(const Tz &z, Ti k = 0)
	{
		Tz w;

		Tc z_imag = std::imag(z);

		if(k==-1)
		{
			if(std::abs(z+(Tc)0.04)<=0.14 && z_imag>=0)
			{
				w = _lambw_N1<Tz,Tc>(z);
			}
			else
			{
				if(abs(z+(Tc)std::exp(-1.l))<0.3 && z_imag>=0)
				{
					w = _lambw_sing<Tz,Tc,Ti>(z,k);
				}
				else
				{
					w = _lambw_asy<Tz,Tc,Ti>(z,-1);
				}
			}
		}
		else
		{
			if(k==1)
			{
				if(std::abs(z+(Tc)0.04)<=0.14 && z_imag<0)
				{
					w = _lambw_N1<Tz,Tc>(z);
				}
				else
				{
					if(std::abs(z+(Tc)std::exp(-1.l))<0.3 && z_imag<0)
					{
						w = _lambw_sing<Tz,Tc,Ti>(z,k);
					}
					else
					{
						w = _lambw_asy<Tz,Tc,Ti>(z,1);
					}
				}
			}
			else
			{
				if(k!=0)
				{
					w = _lambw_asy<Tz,Tc,Ti>(z,k);
				}
				else
				{
					if(std::abs(z-(Tc)0.02)<0.28)
					{
						w = _lambw_laur<Tz,Tc>(z);
					}
					else
					{
						if(std::abs(z+(Tc)std::exp(-1.))<0.3)
						{
							w = _lambw_sing<Tz,Tc,Ti>(z,k);
						}
						else
						{
							if(std::abs(z-(Tc)8.)<13)
							{
								w = _lambw_mid<Tz,Tc>(z);
							}
							else
							{
								w = _lambw_asy<Tz,Tc,Ti>(z,0);
							}
						}
					}
				}
			}
		}

		if( (std::isnan(std::real(w))||std::isnan(std::imag(w))) || (std::real(z)==0 && std::imag(z)==0) )
		{
			return w;
		}

		int prec = _lambw_test_score(z,w);

		//std::cout << prec << std::endl;

		if(prec >= Precision)
		{
			return w;
		}

		while(prec <= Precision/4)
		{
			w = _iter_newt<Tz,Tc,Ti>(z,k,w);
			prec = _lambw_test_score(z,w);
		}

		return _iter_frit<Tz,Tc,Ti>(z,k,w);
	}

	float lambert_w(float z, int k = 0)
	{
		if( (k == 0) || (k == -1) )
		{
			if( (k == -1) && (z == 0) )
			{
				return -std::numeric_limits<float>::infinity();
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
			if( (k == -1) && (z == 0) )
			{
				return -std::numeric_limits<double>::infinity();
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
			if( (k == -1) && (z == 0) )
			{
				return -std::numeric_limits<long double>::infinity();
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

	template<typename Ti>
	std::complex<float> lambert_w(const std::complex<float> &z, Ti k)
	{
		return _lambw_raw<std::complex<float>,float,Ti,8>(z, k);
	}

	template<typename Ti>
	std::complex<double> lambert_w(const std::complex<double> &z, Ti k)
	{
		return _lambw_raw<std::complex<double>,double,Ti,16>(z, k);
	}

	template<typename Ti>
	std::complex<long double> lambert_w(const std::complex<long double> &z, Ti k)
	{
		return _lambw_raw<std::complex<long double>,long double,Ti,24>(z, k);
	}

	template<typename Tz>
	Tz wright_omega(const Tz &z)
	{
		return lambert_w(std::exp(z),(long long)std::round(std::imag(z)/boost::math::constants::two_pi<decltype(std::imag(z))>()));
	}

	template<typename Tz, typename Ti>
	Tz alt_lambert_w(const Tz &z, Ti k)
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

} // namespace czb

#endif // LAMBW_HPP_INCLUDED
