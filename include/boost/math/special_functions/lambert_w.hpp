
//  Copyright Balazs Cziraki 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//Sources:
//
//Lambert W function - Wikipedia, the free encyclopedia
//
//Lagrange inversion theorem - Wikipedia, the free encyclopedia
//
//Numerical Evaluation of the Lambert W Function
//and Application to Generation of Generalized
//Gaussian Noise With Exponent 1/2
// - François Chapeau-Blondeau, Member, IEEE, and Abdelilah Monir
//(p2160 IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 50, NO. 9, SEPTEMBER 2002)
//
//Having Fun with Lambert W(x) Function
// - Darko Veberic a,b,c
//   a University of Nova Gorica, Slovenia
//   b IK, Forschungszentrum Karlsruhe, Germany
//   c J. Stefan Institute, Ljubljana, Slovenia
//(arXiv:1003.1628v1 [cs.MS] 8 Mar 2010)

#ifndef BOOST_MATH_LAMBERT_W_HPP_INCLUDED
#define BOOST_MATH_LAMBERT_W_HPP_INCLUDED

#include <complex>
#include <boost/math/tools/promotion.hpp>
#include "detail/lambert_w_raw.hpp"

namespace boost
{
namespace math
{
    //Returns W(z) with the default policy for real argument z.
    template<class ArgumentType>
    inline typename tools::promote_args<ArgumentType,int>::type
    lambert_w(const ArgumentType& z)
    {
        return lambw::_with_real_range_checks<ArgumentType,int>(z,0,
                                                 policies::policy<>());
    }

    //Returns W_k(z) with the default policy for real argument z.
    template<class ArgumentType, class IndexType>
    inline typename tools::promote_args<ArgumentType,IndexType>::type
    lambert_w(const ArgumentType& z, const IndexType& k)
    {
        return lambw::_with_real_range_checks<ArgumentType,IndexType>(z,k,
                                                 policies::policy<>());
    }

    //Returns W_k(z) with the policy pol for real argument z.
    template<class ArgumentType, class IndexType, class Policy>
    inline typename tools::promote_args<ArgumentType,IndexType>::type
    lambert_w(const ArgumentType& z, const IndexType& k, const Policy& pol)
    {
        BOOST_FPU_EXCEPTION_GUARD

        typedef typename tools::promote_args<ArgumentType,IndexType>::type
            result_type;

        typedef typename policies::evaluation<result_type, Policy>::type
            value_type;

        typedef typename policies::normalise<
            Policy,
            policies::promote_float<false>,
            policies::promote_double<false>,
            policies::discrete_quantile<>,
            policies::assert_undefined<> >::type forwarding_policy;

        return
            policies::checked_narrowing_cast<result_type, forwarding_policy>(
                lambw::_with_real_range_checks<result_type,result_type,Policy>(
                    static_cast<value_type>(z),
                    static_cast<value_type>(k),
                forwarding_policy()),
            "boost::math::lambert_w<%1%,%1%>(%1%,%1%)");
    }

    //Returns W(z) with the default policy for complex argument z.
    template<class ArgumentType>
    inline std::complex<ArgumentType>
    lambert_w(const std::complex<ArgumentType> &z)
    {
        return lambw::_with_complex_range_checks<ArgumentType,int>(z,0,
                                                 policies::policy<>());
    }

    //Return W_k(z) with the default policy for complex argument z.
    template<class ArgumentType, class IndexType>
    inline std::complex<ArgumentType>
    lambert_w(const std::complex<ArgumentType> &z, const IndexType& k)
    {
        return lambw::_with_complex_range_checks<ArgumentType,IndexType>(z,k,
                                                 policies::policy<>());
    }

    //Return W_k(z) with the policy pol for complex argument z.
    template<class ArgumentType, class IndexType, class Policy>
    inline std::complex<ArgumentType>
    lambert_w(const std::complex<ArgumentType>& z, const IndexType& k,
              const Policy& pol)
    {
        BOOST_FPU_EXCEPTION_GUARD

        typedef typename tools::promote_args<ArgumentType, IndexType>::type
            result_type;

        typedef typename policies::evaluation<result_type, Policy>::type
            value_type;

        typedef typename policies::normalise<
            Policy,
            policies::promote_float<false>,
            policies::promote_double<false>,
            policies::discrete_quantile<>,
            policies::assert_undefined<> >::type forwarding_policy;

        return
            policies::checked_narrowing_cast<result_type, forwarding_policy>(
                lambw::_with_complex_range_checks<std::complex<result_type>,
                result_type,Policy>(
                    static_cast<std::complex<value_type> >(z),
                    static_cast<value_type>(k),
                forwarding_policy()),
        "boost::math::lambert_w<std::complex<%1%>,%1%>(std::complex<%1%>,%1%)"
            );
    }
} // namespace math
} // namespace boost

#endif // BOOST_MATH_LAMBERT_W_HPP_INCLUDED
