
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef TESTING_TOOLS_HPP_INCLUDED
#define TESTING_TOOLS_HPP_INCLUDED

#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>

template<class ArgumentType>
ArgumentType random_cauchy_number()
{
     using std::tan;

     ArgumentType r = (ArgumentType)rand()/(ArgumentType)RAND_MAX;
     return tan(boost::math::constants::pi<ArgumentType>()*(r-(ArgumentType)0.5));
}

template<class ArgumentType>
ArgumentType random_unif_number(ArgumentType a, ArgumentType b)
{
     ArgumentType r = (ArgumentType)rand()/(ArgumentType)RAND_MAX;
     return a+r*(b-a);
}

template<class ArgumentType>
ArgumentType comp(ArgumentType x1, ArgumentType x2)
{
     using std::abs;

     if(x1 == x2)
     {
          return (ArgumentType)0.;
     }
     else
     {
          return std::max( x1!=0.?abs((x1-x2)/x1):0. , x2!=0.?abs((x1-x2)/x1):0. )/std::numeric_limits<ArgumentType>::epsilon();
     }
}

template<class ArgumentType>
ArgumentType comp(std::complex<ArgumentType> x1, std::complex<ArgumentType> x2)
{
     using std::abs;

     if(x1 == x2)
     {
          return (ArgumentType)0.;
     }
     else
     {
          return std::max( x1!=std::complex<ArgumentType>(0.,0.)?abs((x1-x2)/x1):0. , x2!=std::complex<ArgumentType>(0.,0.)?abs((x1-x2)/x1):0. )/std::numeric_limits<ArgumentType>::epsilon();
     }
}

template<class ArgumentType, class IndexType>
void test_exact_value(const ArgumentType& x, const ArgumentType& w_act, const IndexType& k)
{
     ArgumentType w = (ArgumentType)0.;
     try
     {
          w = boost::math::lambert_w(x,k);
          ArgumentType err = comp(w,w_act);

          std::cout << "W(" << x << ',' << k << ") = " << w << std::endl;
          std::cout << "(Actual value: " << w_act << ", Error: " << err << " eps)" << std::endl;
     }
     catch(const std::exception& e)
     {
          std::cout << "Exception encountered while evaluating W(" << x << ',' << k << "):" << std::endl;
          std::cout << e.what() << std::endl;
     }
}

template<class ArgumentType, class IndexType>
void test_exact_value(const std::complex<ArgumentType>& x, const std::complex<ArgumentType>& w_act, const IndexType& k)
{
     std::complex<ArgumentType> w;
     try
     {
          w = boost::math::lambert_w(x,k);
          ArgumentType err = comp(w,w_act);

          std::cout << "W(" << x << ',' << k << ") = " << w << std::endl;
          std::cout << "(Actual value: " << w_act << ", Error: " << err << " eps)" << std::endl;
     }
     catch(const std::exception& e)
     {
          std::cout << "Exception encountered while evaluating W(" << x << ',' << k << "):" << std::endl;
          std::cout << e.what() << std::endl;
     }
}

template<class ArgumentType, std::size_t N>
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
          ArgumentType w = boost::math::lambert_w((ArgumentType)(-log(a)/a));
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

     ArgumentType inv_test_err = comp(max_err_x, (ArgumentType)(max_err_w*exp(max_err_w)) );

     std::cout << "(The error of this tested with W(x)*exp(W(x)) = x is "
          << inv_test_err << " eps )" << std::endl;
}

template<class ArgumentType, std::size_t N>
std::pair<ArgumentType,ArgumentType> inverse_identity_test(ArgumentType a, ArgumentType b, bool output = true)
{
     using std::log;
     using std::exp;

     ArgumentType avg_err = (ArgumentType)0.;
     ArgumentType max_err = (ArgumentType)0.;
     ArgumentType max_err_x = (ArgumentType)0.;

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
               max_err_x = w;
          }
     }

     avg_err /= N;

     if(output)
     {
          std::cout << "Inverse identity test results:" << std::endl;
          std::cout << "W(x)*exp(W(x)) = x, " << a << " <= x <= " << b << std::endl;
          std::cout << "Average error: " << avg_err << " eps" << std::endl;
          std::cout << "Max error: " << max_err << " eps" << std::endl;
          std::cout << "(which occurs at -1/e + " << ( max_err_x+boost::math::lambw::rec_e<ArgumentType>() )/std::numeric_limits<ArgumentType>::epsilon() << " eps )" << std::endl;
     }

     return std::make_pair(avg_err,max_err);
}

template<class ArgumentType, std::size_t N>
std::pair<ArgumentType,ArgumentType> gsl_inverse_identity_test(ArgumentType a, ArgumentType b, bool output = true)
{
     using std::log;
     using std::exp;

     ArgumentType avg_err = (ArgumentType)0.;
     ArgumentType max_err = (ArgumentType)0.;
     ArgumentType max_err_x = (ArgumentType)0.;

     for(std::size_t k = 0; k < N; ++k)
     {
          ArgumentType x_test = random_unif_number(a,b);
          ArgumentType w = x_test*exp(x_test);
          ArgumentType x = gsl_sf_lambert_W0(w);
          ArgumentType err = comp(x,x_test);
          avg_err += err;
          if(max_err < err)
          {
               max_err = err;
               max_err_x = w;
          }
     }

     avg_err /= N;

     if(output)
     {
          std::cout << "Comparison with the GSL Lambert W function:" << std::endl;
          std::cout << "W(x)*exp(W(x)) = x, " << a << " <= x <= " << b << std::endl;
          std::cout << "Average error: " << avg_err << " eps" << std::endl;
          std::cout << "Max error: " << max_err << " eps" << std::endl;
          std::cout << "(which occurs at -1/e + " << ( max_err_x+boost::math::lambw::rec_e<ArgumentType>() )/std::numeric_limits<ArgumentType>::epsilon() << " eps )" << std::endl;
     }

     return std::make_pair(avg_err,max_err);
}

void make_complex_gnuplottable(const std::string& infilename, const std::string& outfilename)
{
     std::ifstream infile(infilename.c_str());
     std::ofstream outfile(outfilename.c_str());

     std::string x;
     std::string y;
     std::string u;
     std::string v;

     std::string x_prev;

     infile >> x >> y >> u >> v;

     while(infile)
     {
          x_prev = x;
          outfile << x << '\t' << y << '\t' << u << '\t' << v << std::endl;
          infile >> x >> y >> u >> v;
          if(x_prev != x)
          {
               outfile << std::endl;
          }
     }

     infile.close();
     outfile.close();
}

#endif // TESTING_TOOLS_HPP_INCLUDED
