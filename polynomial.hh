#ifndef POLYNOMIALPINNER_POLYNOMIAL_HH
#define POLYNOMIALPINNER_POLYNOMIAL_HH
#include <vector>
#include <optional>
#include <variant>
#include <map>

struct polynomial;

using variable = int;

struct term;
struct polynomial
{
    std::vector<term> terms;

    static polynomial zero();
    static polynomial create(double value);
    static polynomial create(term single_term);
    static polynomial roots(variable var, polynomial expr);
    static polynomial create(
        const variable* indeterminates,
        size_t dimension,
        unsigned degree,
        variable& variable_counter
    );
};
bool operator<(const polynomial& a, const polynomial& b);
bool operator==(const polynomial& a, const polynomial& b);

struct var_power;
struct term
{
    double coefficient = 0;
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
polynomial assign(const term& t, variable id, double value);
std::optional<term> try_sum(const term& a, const term& b); // Succeeds only if multipliers are common.
polynomial assign(const term& t, variable id, const polynomial& equivalent);
bool compatible(const var_power& a, const var_power& b);

std::optional<polynomial> differentiate(const polynomial& p, variable id);
polynomial multiply(const polynomial& a, const polynomial& b);
polynomial sum(const polynomial& a, const polynomial& b);
polynomial sort(const polynomial& p); // Sorts terms by variable indices and powers.
polynomial simplify(const polynomial& p, double zero_epsilon = 1e-10); // Merges terms with same variables and removes zero terms
polynomial assign(const polynomial& p, variable id, const polynomial& equivalent);
polynomial assign(const polynomial& p, variable id, double value);

std::optional<double> try_get_constant_polynomial_value(const polynomial& p);
bool depends_on_var(const polynomial& p, variable id);
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

// Reduces the number of variables with the knowledge that 'zero' must always be
// zero regardless of the values of 'indeterminates'. The resulting equivalences
// between variables are then assigned to every polynomial in 'target'.
bool pin(
    const polynomial& zero,
    const variable* indeterminates,
    size_t indeterminate_count,
    std::vector<polynomial>& target
);

#endif
