#include <iostream>
#include <fstream>
#include <limits>
#include <boost/math/special_functions/lambert_w.hpp>

double plog(double x)
{
  if (x == 0) {
    return 0;
  }

  double w0, w1;
  if (x > 0) {
    w0 = log(1.2 * x / log(2.4 * x / log1p(2.4 * x)));
  } else {
    double v = 1.4142135623730950488 * sqrt(1 + 2.7182818284590452354 * x);
    double N2 = 10.242640687119285146 + 1.9797586132081854940 * v;
    double N1 = 0.29289321881345247560 * (1.4142135623730950488 + N2);
    w0 = -1 + v * (N2 + v) / (N2 + v + N1 * v);
  }

  int cntr = 1000;
  while (cntr > 0) {
    double e = exp(w0);
    double f = w0 * e - x;
    w1 = w0 - f / ((e * (w0 + 1) - (w0 + 2) * f / (w0+w0 + 2)));
    if (fabs(w0 / w1 - 1) < 1.4901161193847656e-8) {
      break;
    }
    w0 = w1;
    --cntr;
  }
  return w1;
}

int main()
{
    double x_min = -std::exp(-1.);
    double x_max = 0;
    int samp = 1025;
    double dx = (x_max-x_min)/(samp-1);

    std::ofstream file("data.txt");

    for(int k = 0; k < samp; ++k)
	{
		double x = x_min + k*dx;
		file << x << '\t' << std::abs(boost::math::_lambw_test_s(x,boost::math::lambert_w(x)))
					<< '\t' << std::abs(boost::math::_lambw_test_s(x,plog(x))) << std::endl;
	}

    file.close();

    return 0;
}
