#ifndef RATIONALFITTER_OPTIMIZATION_HH
#define RATIONALFITTER_OPTIMIZATION_HH
#include "matrix.hh"
#include "rational.hh"

matrix least_squares(const matrix& A, const matrix& b);

std::optional<double> l2_loss(
    const rational& r,
    variable right_side_variable,
    const std::map<variable, const double*>& data,
    size_t data_size
);

// Fits remaining variables to given data. The returned rational function should
// no longer depend on any other variables than the ones in 'data'.
std::optional<rational> fit(
    const rational& r,
    variable right_side_variable,
    const std::map<variable, const double*>& data,
    size_t data_size,
    const std::map<variable, double>& initial_guesses = {}
);

#endif
