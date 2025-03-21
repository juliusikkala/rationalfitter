#ifndef RATIONALFITTER_POLYNOMIAL_HH
#define RATIONALFITTER_POLYNOMIAL_HH
#include <vector>
#include <optional>
#include <variant>
#include <map>
#include <set>
#include <tuple>
#include "number.hh"

struct polynomial;

using variable = int;

struct term;
struct polynomial
{
    std::vector<term> terms;

    static polynomial zero();
    static polynomial create(number value);
    static polynomial create(term single_term);
    static polynomial roots(variable var, polynomial expr);
    static polynomial create(
        const variable* indeterminates,
        size_t dimension,
        unsigned degree,
        variable& variable_counter,
        bool normalized = false// If true, one less variable is used and the constant is 1.
    );
};
bool operator<(const polynomial& a, const polynomial& b);
bool operator==(const polynomial& a, const polynomial& b);

struct var_power;
struct term
{
    number coefficient = 0;
    std::vector<var_power> mul;
};
bool operator<(const term& a, const term& b);
bool operator==(const term& a, const term& b);

// The value of a 'roots' expr are the roots of 'expr' in terms of 'var'.
struct roots_expr
{
    variable var;
    polynomial expr;
};
bool operator<(const roots_expr& a, const roots_expr& b);
bool operator==(const roots_expr& a, const roots_expr& b);

// This normally represents a variable to some power, e.g. x^4. However, if
// roots.has_value(), the expression is instead roots(expr)^degree. This form
// allows elimination of variables stemming from complicated constraints.
struct var_power
{
    variable id = -1;
    unsigned degree = 0;
    std::optional<roots_expr> roots;
};
bool operator<(const var_power& a, const var_power& b);
bool operator==(const var_power& a, const var_power& b);

std::optional<term> differentiate(const term& t, variable id);
term sort(const term& t); // Sorts multipliers by variable index.

// Many 'term' functions may end up resolving 'roots' for whatever reason, and
// when that happens, the term may get split up into multiple terms. Which is
// why many of these return a polynomial instead of a single term.
polynomial multiply(const term& a, const term& b);

// Merges multiplicands with common variables and resolves roots if possible.
polynomial simplify(const term& t);
polynomial assign(const term& t, variable id, number value);
std::optional<term> try_sum(const term& a, const term& b); // Succeeds only if multipliers are common.
polynomial assign(const term& t, variable id, const polynomial& equivalent);
bool compatible(const var_power& a, const var_power& b);

std::optional<polynomial> differentiate(const polynomial& p, variable id);
polynomial multiply(const polynomial& a, const polynomial& b);
polynomial multiply(const polynomial& a, number v);
polynomial sum(const polynomial& a, const polynomial& b);
polynomial sort(const polynomial& p); // Sorts terms by variable indices and powers.
polynomial simplify(const polynomial& p); // Merges terms with same variables and removes zero terms
polynomial assign(const polynomial& p, variable id, const polynomial& equivalent);
polynomial assign(const polynomial& p, variable id, number value);


// There's two ways this can fail: either 't' contains 'roots' expressions that
// cannot be solved or has multiple solutions, or 'variable_values' doesn't
// actually contain a value for each live variable.
std::optional<number> evaluate(const term& t, const std::vector<number>& variable_values);
std::optional<number> evaluate(const polynomial& p, const std::vector<number>& variable_values);

// Try to factor a polynomial to the form
// (var-root) * (result_polynomial)
// This is necessary to simplify rationals and solve nth derivative constraints
// on them.
std::optional<polynomial> try_factor(const polynomial& p, variable id, const polynomial& root);

std::optional<number> try_get_constant_value(const polynomial& p);
bool depends_on_var(const polynomial& p, variable id);
std::set<variable> live_variables(const polynomial& p);
// Tries to solve the roots raised to 'exponent', and returns a value if it
// found a single root. The exponent is meaningful to do here, because e.g. if
// the root is +-2 but exponent is 2, there's only one possible value anyway
// (4). Bigger number == less severe issue.
enum class solve_failure_reason
{
    NO_ROOTS = 0,
    INFINITE_ROOTS = 1,
    CANT_SOLVE = 2,
    CANT_REPRESENT = 3,
    MULTIPLE_ROOTS = 4
};
std::variant<polynomial, solve_failure_reason> try_solve_single_root(const polynomial& p, variable id, unsigned exponent);
// This function may return _any_ of the roots of polynomial 'p' if it can.
std::optional<polynomial> try_find_any_root(const polynomial& p, variable id);
// List of roots + residual
struct roots_result
{
    bool found_all = false;
    std::vector<polynomial> roots;
    polynomial residual;
};
roots_result try_find_all_roots(const polynomial& p, variable id);
std::optional<polynomial> try_get_nth_root(const polynomial& constant, unsigned N);
void find_common_variables(const polynomial& p, std::map<variable, int>& common, bool continue_from_previous = false);
polynomial factor_common_variables(const polynomial& p, const std::map<variable, int>& common);

struct indeterminate_group
{
    std::vector<var_power> indeterminates;
};
bool operator<(const indeterminate_group& a, const indeterminate_group& b);

std::map<indeterminate_group, polynomial> group_by_indeterminates(
    const polynomial& p,
    const variable* indeterminates,
    size_t indeterminate_count
);

polynomial get_zero_polynomial(const polynomial& a, number right_side);

// Reduces the number of variables with the knowledge that 'zero' must always be
// zero regardless of the values of 'indeterminates'. The resulting equivalences
// between variables are then assigned to every polynomial in 'target'.
//
// Additionally, this function ensures that 'nonzero' is non-zero. This is
// necessary to filter out solutions that break a division. If you don't need
// this feature, you can set it to polynomial::create(1.0).
bool pin(
    const polynomial& zero,
    const polynomial& nonzero,
    const variable* indeterminates,
    size_t indeterminate_count,
    std::vector<polynomial>& target
);

#endif
