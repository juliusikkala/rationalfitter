#ifndef RATIONALFITTER_RATIONAL_HH
#define RATIONALFITTER_RATIONAL_HH
#include "polynomial.hh"

struct rational
{
    polynomial numerator;
    polynomial denominator;
};

std::optional<rational> differentiate(const rational& r, variable id);
rational simplify(const rational& r);

rational assign(const rational& r, variable id, double value);
std::set<variable> live_variables(const rational& r);

polynomial get_zero_polynomial(const rational& r, double right_side);

#endif
