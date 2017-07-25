
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
#define BOOST_MATH_POLE_ERROR_POLICY ignore_error
#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error
#define BOOST_MATH_UNDERFLOW_ERROR_POLICY ignore_error
#define BOOST_MATH_DENORMALISATION_ERROR_POLICY ignore_error
#define BOOST_MATH_ROUNDING_ERROR_POLICY ignore_error
#define BOOST_MATH_EVALUATION_ERROR_POLICY ignore_error
#define BOOST_MATH_INDETERMINATE_RESULT_ERROR_POLICY ignore_error
/*
#define BOOST_MATH_DOMAIN_ERROR_POLICY throw_on_error
#define BOOST_MATH_POLE_ERROR_POLICY throw_on_error
#define BOOST_MATH_OVERFLOW_ERROR_POLICY throw_on_error
#define BOOST_MATH_UNDERFLOW_ERROR_POLICY throw_on_error
#define BOOST_MATH_DENORMALISATION_ERROR_POLICY throw_on_error
#define BOOST_MATH_ROUNDING_ERROR_POLICY throw_on_error
#define BOOST_MATH_EVALUATION_ERROR_POLICY throw_on_error
#define BOOST_MATH_INDETERMINATE_RESULT_ERROR_POLICY throw_on_error
*/

#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <ctime>
#include <sstream>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/tools/test_data.hpp>
#include <gsl/gsl_sf_lambert.h>
#include "testing_tools.hpp"

typedef boost::multiprecision::float128 MultiPrec;
typedef double RealType;
typedef int IntegralType;

const bool precision_tests = true;
const std::size_t test_count = 1000;

const bool benchmark = false;
const std::size_t benchmark_count = 1000000;

template<class CoeffType>
CoeffType p(CoeffType z, int k)
{
     return ((CoeffType)pow(-1.,k%2))*sqrt((CoeffType)2.*(boost::math::constants::e<CoeffType>()*z+(CoeffType)1.));
}

template<class ArgumentType>
ArgumentType test_lambw(const ArgumentType& z)
{
     using std::exp;
     using std::abs;
     ArgumentType w = boost::math::lambert_w(z);
     ArgumentType t = w*exp(w);

     return (z==t)?(ArgumentType)0.:std::max(abs((z-t)/z),abs((z-t)/t))/std::numeric_limits<ArgumentType>::epsilon();
}

int main()
{
     std::srand(std::time(NULL));

     if(precision_tests)
     {
          RealType x;
          RealType w_act;

          using std::sqrt;
          std::cout << comp(boost::math::lambw::_sing<RealType,RealType,IntegralType>(-boost::math::lambw::rec_e<RealType>(),(IntegralType)0),(RealType)(-1.)) << std::endl;
          std::cout << std::endl;

          //W(-pi/2) = i*pi/2
          test_exact_value(std::complex<RealType>(-boost::math::constants::half_pi<RealType>(),(RealType)0.),
                              std::complex<RealType>((RealType)0.,boost::math::constants::half_pi<RealType>()),0);
          std::cout << std::endl;

          //W(-1/e) = -1
          x = -boost::math::lambw::rec_e<RealType>();
          w_act = (RealType)(-1.);
          test_exact_value(x,w_act,0);
          std::cout << std::endl;

          //W(-1/e,-1) = -1
          x = -boost::math::lambw::rec_e<RealType>();
          w_act = (RealType)(-1.);
          test_exact_value(x,w_act,-1);
          std::cout << std::endl;

          //W(0) = 0
          x = (RealType)0.;
          w_act = (RealType)0.;
          test_exact_value(x,w_act,0);
          std::cout << std::endl;

          test_exact_value(std::complex<RealType>((RealType)0.,(RealType)0.),std::complex<RealType>((RealType)0.,(RealType)0.),0);
          std::cout << std::endl;

          //W(0,-1) = -inf
          x = (RealType)0.;
          w_act = -std::numeric_limits<RealType>::infinity();
          test_exact_value(x,w_act,-1);
          std::cout << std::endl;

          test_exact_value(std::complex<RealType>((RealType)0.,(RealType)0.),std::complex<RealType>(std::numeric_limits<RealType>::quiet_NaN(),std::numeric_limits<RealType>::quiet_NaN()),-1);
          std::cout << std::endl;

          //W(1) = Omega constant
          x = (RealType)1.;
          RealType omega;
          std::stringstream omega_text("0.567143290409783872999968662210355549753815787186512508135131079223045793086684566693219446961752294557638");
          omega_text >> omega;
          test_exact_value(x,omega,0);
          std::cout << std::endl;

          //Index bound test
          w_act = std::numeric_limits<RealType>::quiet_NaN();
          test_exact_value(x,w_act,(IntegralType)1);
          std::cout << std::endl;

          //W(e) = 1
          x = boost::math::constants::e<RealType>();
          w_act = (RealType)1.;
          test_exact_value(x,w_act,0);
          std::cout << std::endl;

          //W(x->inf) -> inf
          x = std::numeric_limits<RealType>::infinity();
          w_act = std::numeric_limits<RealType>::infinity();
          test_exact_value(x,w_act,0);
          std::cout << std::endl;

          //Range bound test
          x = (RealType)1.;
          w_act = std::numeric_limits<RealType>::quiet_NaN();
          test_exact_value(x,w_act,(IntegralType)(-1));
          std::cout << std::endl;

          std::cout << std::endl;

          ln_identity_test<RealType,test_count>();

          std::cout << std::endl;

          inverse_identity_test<RealType,test_count>(-1,10);

          std::cout << std::endl;

          gsl_inverse_identity_test<double,test_count>(-1,10);
     }

     if(benchmark)
     {
          std::clock_t t0 = std::clock();

          for(std::size_t k = 0; k < benchmark_count; ++k)
          {
               boost::math::lambert_w(random_unif_number<RealType>(0,10));
               //gsl_sf_lambert_W0(random_unif_number<RealType>(0,10));
          }

          std::clock_t t1 = std::clock();

          std::clock_t diff = t1-t0;

          std::cout << "Calculating W for " << benchmark_count << " random values took " << (double)diff/(double)CLOCKS_PER_SEC << " seconds. (CPU time)" << std::endl;
          std::cout << "Effectively " << (double)diff/(double)CLOCKS_PER_SEC/benchmark_count << " seconds/call." << std::endl;
     }

     return 0;
}
