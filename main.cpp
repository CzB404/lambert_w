
//          Copyright Balazs Cziraki 2016
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <chrono>
#include <random>
#include "lambw.hpp"
#include "datafile.hpp"

using RealType = long double;
using IntegralType = int;

int main()
{
	using namespace std::literals::complex_literals;

	std::chrono::steady_clock clock;

	std::cout << "W(-1/e) = " << czb::lambert_w(-std::exp(-1)) << std::endl;
	std::cout << "W(-1/e,-1) = " << czb::lambert_w(-std::exp(-1),-1) << std::endl;
	std::cout << "W(0) = " << czb::lambert_w(0.) << std::endl;
	std::cout << "W(0,-1) = " << czb::lambert_w(0.,-1) << std::endl;
	std::cout << "W(-pi/2) = " << czb::lambert_w(-boost::math::constants::half_pi<double>()+0i) << std::endl;
	std::cout << "W(1) = " << czb::lambert_w(1.) << std::endl;
	std::cout << "W(42) = " << czb::lambert_w(42.) << std::endl;
	std::cout << "W(42)*exp(W(42)) = " << czb::lambert_w(42.)*std::exp(czb::lambert_w(42.)) << std::endl;

	std::random_device rd;
	std::default_random_engine rnd(rd());
	std::cauchy_distribution<double> dist;

	const std::size_t L = 1000000;

	std::array<double,L> w;

	auto t0 = clock.now();

	for(std::size_t k = 0; k < L; ++k)
	{
		w[k] = czb::lambert_w(dist(rnd));
	}

	auto t1 = clock.now();

	std::chrono::duration<double> diff = t1-t0;

	std::cout << "Calculating W for " << L << " random values took " << diff.count() << " seconds." << std::endl;

	RealType d = 3;

	IntegralType samp = 129;

	czb::Plot<std::complex<RealType>> p
	(
		[&](std::complex<RealType> z)
		{
			return czb::_lambw_test_s(z,czb::alt_lambert_w(z,250));
		}
	,samp,d,"lin");

	std::ofstream file("test.txt");

	file.precision(12);

	file << p;

	file.close();

	file.open("settings.plt");

	file << p.settings();

	file.close();

	std::system("wgnuplot complex_datafile.plt");

    return 0;
}
