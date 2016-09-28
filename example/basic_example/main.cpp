
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <complex>
//#include <chrono>
#include <ctime>
//#include <random>
#include <cstdlib>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

//using RealType = long double;
//using IntegralType = int;

typedef double RealType;
typedef int IntegralType;
typedef boost::multiprecision::cpp_bin_float_100 MultiPrec;

RealType dist()
{
	RealType r = (RealType)rand()/(RealType)RAND_MAX;
	return std::tan(boost::math::constants::pi<RealType>()*(r-(RealType)0.5));
}

int main()
{
	//using namespace std::literals::complex_literals;

	//std::chrono::steady_clock clock;

	std::cout << "W(-1/e) = " << boost::math::lambert_w(-std::exp(-1)) << std::endl;
	std::cout << "W(-1/e,-1) = " << boost::math::lambert_w(-std::exp(-1),-1) << std::endl;
	std::cout << "W(0) = " << boost::math::lambert_w(0.) << std::endl;
	std::cout << "W(0,-1) = " << boost::math::lambert_w(0.,-1) << std::endl;
	std::cout << "W(-pi/2+0i) = " << boost::math::lambert_w(std::complex<double>(-boost::math::constants::half_pi<double>(),0)) << std::endl;
	std::cout << "W(1) = " << boost::math::lambert_w(1.) << std::endl;
	std::cout << "W(inf) = " << boost::math::lambert_w(std::numeric_limits<double>::infinity()) << std::endl;
	std::cout << "W(-1) = " << boost::math::lambert_w(-1.) << std::endl;
	std::cout << "W(-1,-1) = " << boost::math::lambert_w(-1.,-1) << std::endl;
	std::cout << "W(1,-1) = " << boost::math::lambert_w(1.,-1) << std::endl;
	std::cout << "W(1+0i,-1) = " << boost::math::lambert_w(std::complex<double>(1,0),-1) << std::endl;
	std::cout << "W(42) = " << boost::math::lambert_w(42.) << std::endl;
	std::cout << "W(42)*exp(W(42)) = " << boost::math::lambert_w(42.)*std::exp(boost::math::lambert_w(42.)) << std::endl;

	MultiPrec x = 42;
	MultiPrec w = boost::math::lambert_w(x);

	std::cout << "W(42) (Multiprecision) = " << w << std::endl;


	//std::random_device rd;
	//std::default_random_engine rnd;//(rd());
	//std::cauchy_distribution<double> dist;

	const std::size_t L = 1000000;

	//std::array<double,L> w;

	//auto t0 = std::chrono::steady_clock::now();
	std::time_t t0 = std::time(NULL);

	for(std::size_t k = 0; k < L; ++k)
	{
		boost::math::lambert_w(dist(/*rnd*/));
	}

	//auto t1 = std::chrono::steady_clock::now();
	std::time_t t1 = std::time(NULL);

	//std::chrono::duration<double> diff = t1-t0;
	std::time_t diff = t1-t0;

	std::cout << "Calculating W for " << L << " random values took " << diff/*.count()*/ << " seconds." << std::endl;

    return 0;
}
