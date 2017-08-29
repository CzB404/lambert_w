
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
This example creates two files, each containing the complex numeric results of
the Lambert W function on certain branches.

"lw_0.txt" will contain W_{0}, and "lw_pm1.txt" will contain W_{-1} on the
positive imaginary half plane and W_{1} on the negative imaginary half plane,
showing that they join smoothly on the line [-e^{-1}:0].

The supplied "plot_complex_0pm1.plt" gnuplot script creates a png image from the
two files, showing how the three branches smoothly join together and how the
surface of the complex Lambert W function passes through itself.
*/

#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
#define BOOST_MATH_POLE_ERROR_POLICY ignore_error
#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error
#define BOOST_MATH_UNDERFLOW_ERROR_POLICY ignore_error
#define BOOST_MATH_DENORMALISATION_ERROR_POLICY ignore_error
#define BOOST_MATH_ROUNDING_ERROR_POLICY ignore_error
#define BOOST_MATH_EVALUATION_ERROR_POLICY ignore_error
#define BOOST_MATH_INDETERMINATE_RESULT_ERROR_POLICY ignore_error


#include <iostream>
#include <fstream>
#include <boost/math/tools/test_data.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/lambert_w.hpp>

typedef boost::multiprecision::float128 MultiPrec;//Select a multiprecision type
typedef double RealType; //Select a real type to use (can be MultiPrec)
typedef int IntegralType; //Select an integral type to use

//Complex Lambert W tester function
template<class ArgumentType>
std::pair<ArgumentType,ArgumentType>
test_complex_lambw(const ArgumentType& x, const ArgumentType& y)
{
     using std::exp;
     using std::abs;
     std::complex<ArgumentType> z(x,y);
     std::complex<ArgumentType> w = boost::math::lambert_w(z);
     std::complex<ArgumentType> t = w*exp(w);

     ArgumentType ans = (z==t)?(ArgumentType)0.:
          std::max(abs((z-t)/z),abs((z-t)/t))
          /std::numeric_limits<ArgumentType>::epsilon();
     return std::make_pair(ans,(ArgumentType)0.);
}

//Wraps the std::complex boost::math::lambert_w
//with index 0 to use in boost::math::tools::test_data
template<class ArgumentType>
std::pair<ArgumentType,ArgumentType>
plot_complex_lambw(const ArgumentType& x, const ArgumentType& y)
{
     std::complex<ArgumentType> w =
          boost::math::lambert_w(std::complex<ArgumentType>(x,y));
     return std::make_pair(w.real(),w.imag());
}

//Similar to plot_complex_lambw, but instead it cuts W_{-1} and W_{1} together
//so that they join smoothly in the [-e^{-1}:0] line on the complex plane.
template<class ArgumentType>
std::pair<ArgumentType,ArgumentType>
plot_complex_alt_lambw(const ArgumentType& x, const ArgumentType& y)
{
     static const IntegralType k = -1;
     std::complex<ArgumentType> w =
          (y<0)?
               boost::math::lambert_w(std::complex<ArgumentType>(x,y),-k)
               :boost::math::lambert_w(std::complex<ArgumentType>(x,y),k);
     return std::make_pair(w.real(),w.imag());
}

//Makes the content of the file with name infilename 3D plottable by gnuplot
//and writes the result in a file with name outfilename.
void make_complex_gnuplottable(const std::string& infilename,
                               const std::string& outfilename)
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

int main()
{
     //Plot parameters
     //d1: Positive axis end
     //d0: Negative axis end
     //samp: Sampling count
     using std::exp;
     RealType d1 = 1.0;
     RealType d0 = -d1;
     const int samp = 1024;

     boost::math::tools::test_data<RealType> data;
     std::pair<RealType,RealType> (*pf)(const RealType&, const RealType&) =
          plot_complex_lambw;


     data.insert(pf, boost::math::tools::make_periodic_param(d0,d1,samp),
                 boost::math::tools::make_periodic_param(d0,d1,samp));

     std::ofstream file("test_raw.txt");
     file.precision(std::numeric_limits<RealType>::digits10+2);
     boost::math::tools::write_csv(file, data, "\t");
     file.close();

     make_complex_gnuplottable("test_raw.txt","lw_0.txt");

     boost::math::tools::test_data<RealType> data2;
     std::pair<RealType,RealType> (*pg)(const RealType&, const RealType&) =
          plot_complex_alt_lambw;
     data2.insert(pg, boost::math::tools::make_periodic_param(d0,d1,samp),
                  boost::math::tools::make_periodic_param(d0,d1,samp));

     file.open("test_raw.txt");
     file.precision(std::numeric_limits<RealType>::digits10+2);
     boost::math::tools::write_csv(file, data2, "\t");
     file.close();

     make_complex_gnuplottable("test_raw.txt","lw_pm1.txt");

     return 0;
}
