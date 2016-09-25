# lambert_w
A C++ library that can compute the Lambert W function for floating point and complex number types.
The Lambert W function is a trancendental function that can solve the equation W*exp(W)=z for W.

More information can be found here:
https://en.wikipedia.org/wiki/Lambert_W_function

For now the library is slowly being prepared to be submitted to the Boost review board.

The library uses the C++11 standard and requires the compiler to support that standard.
It also uses the Boost Math Constants library.

The function itself can be called with `boost::math::lambert_w(ArgumentType z, IndexType k)`,
where `ArgumentType` can be `float`, `double` or `long double`, and their `std::complex` variants.
`IndexType` is expected to be an integer type.
