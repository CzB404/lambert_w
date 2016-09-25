
//          Copyright Balazs Cziraki 2016
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef DATAFILE_HPP_INCLUDED
#define DATAFILE_HPP_INCLUDED

#include <ostream>
#include <complex>
#include <string>
#include <sstream>
#include <functional>

namespace czb
{
	template<typename Tz = std::complex<double>>
	class Plot
	{
	public:
		Plot(const std::function<Tz(Tz)> &&func, int x_samp, int y_samp,
			double xrange_min, double xrange_max,
			double yrange_min, double yrange_max,
			//double zrange_min, double zrange_max,
			const std::string &&scale = "lin") :
				_func(func), _x_samp(x_samp), _y_samp(y_samp),
				_xrange_min(xrange_min), _xrange_max(xrange_max),
				_yrange_min(yrange_min), _yrange_max(yrange_max),
				//_zrange_min(zrange_min), _zrange_max(zrange_max),
				_scale(scale)
		{

		}

		Plot(const std::function<Tz(Tz)> &&func, int samp = 129, double d = 10, /*double dz = 20,*/ const std::string &&scale = "lin") :
				_func(func), _x_samp(samp), _y_samp(samp),
				_xrange_min(-d), _xrange_max(d),
				_yrange_min(-d), _yrange_max(d),
				//_zrange_min(0), _zrange_max(dz),
				_scale(scale)
		{

		}

		void print(std::ostream &out)
		{
			double x_step = (_xrange_max-_xrange_min)/(_x_samp-1);
			double y_step = (_yrange_max-_yrange_min)/(_y_samp-1);

			for(int k = 0; k < _x_samp; ++k)
			{
				double x = _xrange_min+x_step*k;
				for(int l = 0; l < _y_samp; ++l)
				{
					double y = _yrange_min+y_step*l;

					Tz z(x,y);

					Tz w = _func(z);

					out << x << '\t' << y << '\t' << w.real() << '\t' << w.imag() << "\r\n";
				}
				out << "\r\n";
			}
		}

		std::string settings()
		{
			std::stringstream text;
			text << "set xrange[" << _xrange_min << ':' << _xrange_max << "]\r\n";
			text << "set yrange[" << _yrange_min << ':' << _yrange_max << "]\r\n";

			if(_scale == "lin")
			{
				text << "unset logscale z\r\n";
				//text << "set zrange[" << _zrange_min << ':' << _zrange_max << "]\r\n";
				text << "set autoscale z\r\n";
			}
			else if(_scale == "log")
			{
				text << "set logscale z\r\n";
				text << "set autoscale z\r\n";
			}

			return text.str();
		}

	private:
		std::function<Tz(Tz)> _func;

		int _x_samp;
		int _y_samp;

		double _xrange_min;
		double _xrange_max;
		double _yrange_min;
		double _yrange_max;

		double _zrange_min;
		double _zrange_max;

		std::string _scale;
	};

	template<typename Tz>
	std::ostream& operator<<(std::ostream& out, Plot<Tz> p)
	{
		p.print(out);
		return out;
	}
}

#endif // DATAFILE_HPP_INCLUDED
