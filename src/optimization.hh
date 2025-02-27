#ifndef RATIONALFITTER_OPTIMIZATION_HH
#define RATIONALFITTER_OPTIMIZATION_HH
#include "matrix.hh"
#include "rational.hh"
#include <functional>

matrix least_squares(const matrix& A, const matrix& b);

// Fits remaining variables to given data. The returned rational function should
// no longer depend on any other variables than the ones in 'data'.
struct fit_params
{
    variable right_side_variable;
    std::map<variable, const number*> data;
    size_t data_size;
    std::map<variable, number> initial_guesses;
    double nlls_step_size;
    unsigned nlls_max_iterations;
    double nlls_convergence_limit;
};

// Takes coefficients from the given vector and 'data'.
std::optional<number> l2_loss(
    const rational& r,
    const fit_params& params,
    const std::vector<number>& coefficients = {}
);

std::optional<rational> fit(const rational& func, const fit_params& params);
std::optional<rational> fit_eliminate_variable(const rational& func, const fit_params& params, std::vector<variable>& eliminated_variables);

#endif
