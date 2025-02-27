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

rational assign(const rational& r, variable id, number value);
std::set<variable> live_variables(const rational& r);

polynomial get_zero_polynomial(const rational& r, number right_side);
std::optional<number> evaluate(const rational& r, const std::vector<number>& variable_values);

#endif
