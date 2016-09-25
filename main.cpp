#include <iostream>
#include <complex>
#include <chrono>
#include <random>
#include "lambert_w.hpp"

using RealType = long double;
using IntegralType = int;

int main()
{
	//using namespace std::literals::complex_literals;

	std::chrono::steady_clock clock;

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

	std::random_device rd;
	std::default_random_engine rnd(rd());
	std::cauchy_distribution<double> dist;

	const std::size_t L = 10000;

	std::array<double,L> w;

	auto t0 = clock.now();

	for(std::size_t k = 0; k < L; ++k)
	{
		w[k] = boost::math::lambert_w(dist(rnd));
	}

	auto t1 = clock.now();

	std::chrono::duration<double> diff = t1-t0;

	std::cout << "Calculating W for " << L << " random values took " << diff.count() << " seconds." << std::endl;

    return 0;
}
