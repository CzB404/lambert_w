
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
#include <complex>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
//#include <boost/math/tools/test.hpp>

typedef boost::multiprecision::cpp_bin_float_50 MultiPrec;
typedef double RealType;
typedef int IntegralType;

const bool real_type_tests = true;
const bool multiprec_tests = false;
const bool benchmark = false;

template<typename ArgumentType>
ArgumentType random_cauchy_number()
{
     using std::tan;

     ArgumentType r = (ArgumentType)rand()/(ArgumentType)RAND_MAX;
     return tan(boost::math::constants::pi<ArgumentType>()*(r-(ArgumentType)0.5));
}

template<typename ArgumentType>
ArgumentType random_unif_number(ArgumentType a, ArgumentType b)
{
     ArgumentType r = (ArgumentType)rand()/(ArgumentType)RAND_MAX;
     return a+r*(b-a);
}

template<typename ArgumentType>
ArgumentType comp(ArgumentType x1, ArgumentType x2)
{
     using std::abs;

     if(x1 == x2)
     {
          return (ArgumentType)0.;
     }
     else
     {
          return std::max( abs((x1-x2)/x1) , abs((x1-x2)/x1) )/std::numeric_limits<ArgumentType>::epsilon();
     }
}

template<typename ArgumentType>
ArgumentType comp(std::complex<ArgumentType> x1, std::complex<ArgumentType> x2)
{
     using std::abs;
     using std::arg;

     if(x1 == x2)
     {
          return (ArgumentType)0.;
     }
     else
     {
          return std::max( comp(abs(x1),abs(x2)) , comp(arg(x1),arg(x2)) );
     }
}

template<typename ArgumentType>
void test_exact_value(ArgumentType x, ArgumentType w_act, IntegralType k = 0)
{
     using std::abs;

     ArgumentType w;
     try
     {
          w = boost::math::lambert_w(x,k);
          BOOST_TYPEOF(abs(x)) err = comp(w,w_act);

          std::cout << "W(" << x << ',' << k << ") = " << w << std::endl;
          std::cout << "(Actual value: " << w_act << ", Error: " << err << " eps)" << std::endl;
     }
     catch(const std::exception& e)
     {
          std::cout << "Exception encountered while evaluating W(" << x << ',' << k << "):" << std::endl;
          std::cout << e.what() << std::endl;
     }

}

template<typename ArgumentType, std::size_t N>
void ln_identity_test()
{
     using std::log;
     using std::exp;

     ArgumentType sum_err = (ArgumentType)0.;
     ArgumentType max_err = (ArgumentType)0.;

     ArgumentType max_err_x = std::numeric_limits<ArgumentType>::quiet_NaN();

     for(std::size_t k = 0; k < N; ++k)
     {
          ArgumentType a = random_unif_number(
               boost::math::lambw::rec_e<ArgumentType>(),
               boost::math::constants::e<ArgumentType>());
          ArgumentType w = boost::math::lambert_w(-log(a)/a);
          ArgumentType w_act = -log(a);
          ArgumentType err = comp(w,w_act);
          sum_err += err;
          if(max_err < err)
          {
               max_err = err;
               max_err_x = -log(a)/a;
          }
     }

     std::cout << "W(-ln(a)/a) = -ln(a), 1/e <= a <= e" << std::endl;
     std::cout << "Average error: " << sum_err/N << " eps" << std::endl;
     std::cout << "Max error: " << max_err << " eps" << std::endl;
     std::cout << "(at x = " << max_err_x << " )" << std::endl;
     std::cout << "(which is -1/e + " << ( max_err_x+boost::math::lambw::rec_e<ArgumentType>() )/std::numeric_limits<ArgumentType>::epsilon() << " eps )" << std::endl;

     ArgumentType max_err_w = boost::math::lambert_w(max_err_x);

     std::cout << "(The error of this tested with W(x)*exp(W(x)) = x is "
          << comp(max_err_x,max_err_w*exp(max_err_w) ) << " eps )" << std::endl;
}

template<typename ArgumentType, std::size_t N>
std::pair<ArgumentType,ArgumentType> inverse_identity_test(ArgumentType a, ArgumentType b, bool output = true)
{
     using std::log;
     using std::exp;

     ArgumentType avg_err = (ArgumentType)0.;
     ArgumentType max_err = (ArgumentType)0.;

     for(std::size_t k = 0; k < N; ++k)
     {
          ArgumentType x_test = random_unif_number(a,b);
          ArgumentType w = x_test*exp(x_test);
          ArgumentType x = boost::math::lambert_w(w);
          ArgumentType err = comp(x,x_test);
          avg_err += err;
          if(max_err < err)
          {
               max_err = err;
          }
     }

     avg_err /= N;

     if(output)
     {
          std::cout << "W(x)*exp(W(x)) = x, " << a << " <= x <= " << b << std::endl;
          std::cout << "Average error: " << avg_err << " eps" << std::endl;
          std::cout << "Max error: " << max_err << " eps" << std::endl;
     }

     return std::make_pair(avg_err,max_err);
}

template<typename ArgumentType>
void test_entire_range()
{

}

int main()
{
     std::srand(std::time(NULL));

     if(real_type_tests)
     {
          //W(-pi/2) = i*pi/2
          test_exact_value(std::complex<RealType>(-boost::math::constants::half_pi<RealType>(),(RealType)0.),
                              std::complex<RealType>((RealType)0.,boost::math::constants::half_pi<RealType>()));
          std::cout << std::endl;
          //W(-1/e) = -1
          test_exact_value(-boost::math::lambw::rec_e<RealType>(),(RealType)(-1.));
          std::cout << std::endl;
          //W(-1/e,-1) = -1
          test_exact_value(-boost::math::lambw::rec_e<RealType>(),(RealType)(-1.),-1);
          std::cout << std::endl;
          //W(0) = 0
          test_exact_value((RealType)0.,(RealType)0.);
          std::cout << std::endl;
          //W(0,-1) = -inf
          test_exact_value((RealType)0.,-std::numeric_limits<RealType>::infinity(),-1);
          std::cout << std::endl;
          test_exact_value(std::complex<RealType>((RealType)0.,(RealType)0.),std::complex<RealType>(std::numeric_limits<RealType>::quiet_NaN(),std::numeric_limits<RealType>::quiet_NaN()),-1);
          std::cout << std::endl;
          //W(1) = Omega constant
          RealType omega;
          std::stringstream omega_text("0.567143290409783872999968662210355549753815787186512508135131079223045793086684566693219446961752294557638");
          omega_text >> omega;
          test_exact_value((RealType)1.,omega);
          std::cout << std::endl;
          test_exact_value((RealType)1.,std::numeric_limits<RealType>::quiet_NaN(),(IntegralType)1);
          std::cout << std::endl;
          //W(e) = 1
          test_exact_value(boost::math::constants::e<RealType>(),(RealType)1.);
          std::cout << std::endl;
          //W(x->inf) -> inf
          test_exact_value(std::numeric_limits<RealType>::infinity(),std::numeric_limits<RealType>::infinity());
          std::cout << std::endl;
          test_exact_value((RealType)1.,std::numeric_limits<RealType>::quiet_NaN(),(IntegralType)(-1));
          std::cout << std::endl;

          std::cout << std::endl;

          ln_identity_test<RealType,10000>();

          std::cout << std::endl;

          inverse_identity_test<RealType,10000>(-boost::math::lambw::rec_e<RealType>(),100);

          /*
          std::cout << "W(-1) = " << boost::math::lambert_w(-1.) << std::endl;
          std::cout << "W(-1,-1) = " << boost::math::lambert_w(-1.,-1) << std::endl;
          std::cout << "W(1,-1) = " << boost::math::lambert_w(1.,-1) << std::endl;
          std::cout << "W(1+0i,-1) = " << boost::math::lambert_w(std::complex<RealType>(1,0),-1) << std::endl;
          std::cout << "W(42) = " << boost::math::lambert_w(42.) << std::endl;
          std::cout << "W(42)*exp(W(42)) = " << boost::math::lambert_w(42.)*std::exp(boost::math::lambert_w(42.)) << std::endl;
          */
     }

     if(multiprec_tests)
     {
          MultiPrec x = -exp(-(MultiPrec)1.);
          MultiPrec w = boost::math::lambert_w(x,-1);
          MultiPrec err = boost::math::lambw::_test_s(x,w);

          std::cout << "W(" << x << ") (Multiprecision) = " << w << " (With an error of " << err << " )" << std::endl;
     }

     if(benchmark)
     {
          const std::size_t L = 1000000;

          std::time_t t0 = std::time(NULL);

          for(std::size_t k = 0; k < L; ++k)
          {
               boost::math::lambert_w(random_unif_number<RealType>(0,10));
          }

          std::time_t t1 = std::time(NULL);

          std::time_t diff = t1-t0;

          std::cout << "Calculating W for " << L << " random values took " << diff << " seconds." << std::endl;
     }

    return 0;
}
