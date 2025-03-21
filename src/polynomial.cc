#include "polynomial.hh"
#include "rational.hh"
#include <algorithm>
#include <cmath>
#include <climits>
#include <cstdint>
#include <cstdio>

bool operator<(const polynomial& a, const polynomial& b)
{
    if(a.terms.size() < b.terms.size())
        return true;
    if(b.terms.size() < a.terms.size())
        return false;

    for(unsigned i = 0; i < a.terms.size(); ++i)
    {
        if(a.terms[i] < b.terms[i])
            return true;
        if(b.terms[i] < a.terms[i])
            return false;
    }
    return false;
}

bool operator==(const polynomial& a, const polynomial& b)
{
    return a.terms == b.terms;
}

bool operator<(const term& a, const term& b)
{
    if(a.mul.size() < b.mul.size())
        return true;
    if(b.mul.size() < a.mul.size())
        return false;

    for(unsigned i = 0; i < a.mul.size(); ++i)
    {
        if(a.mul[i] < b.mul[i])
            return true;
        if(b.mul[i] < a.mul[i])
            return false;
    }

    return a.coefficient < b.coefficient;
}

bool operator==(const term& a, const term& b)
{
    return a.coefficient == b.coefficient && a.mul == b.mul;
}

bool operator<(const roots_expr& a, const roots_expr& b)
{
    if(a.var < b.var)
        return true;
    if(b.var < a.var)
        return false;
    return a.expr < b.expr;
}

bool operator==(const roots_expr& a, const roots_expr& b)
{
    return a.var == b.var && a.expr == b.expr;
}

bool operator<(const var_power& a, const var_power& b)
{
    if(a.id < b.id)
       return true;
    if(b.id < a.id)
        return false;
    if(a.roots.has_value() && !b.roots.has_value())
        return true;
    if(b.roots.has_value() && !a.roots.has_value())
        return false;
    if(a.roots.has_value() && b.roots.has_value())
    {
        if(*a.roots < *b.roots)
            return true;
        else
            return false;
    }
    return a.degree < b.degree;
}

bool operator==(const var_power& a, const var_power& b)
{
    if(a.id != b.id || a.degree != b.degree)
        return false;
    if(a.roots.has_value() != b.roots.has_value())
        return false;
    if(a.roots.has_value() && !(*a.roots == *b.roots))
        return false;
    return true;
}

polynomial polynomial::zero()
{
    return polynomial{};
}

polynomial polynomial::create(number value)
{
    polynomial p;
    p.terms.push_back(term{value, {}});
    return p;
}

polynomial polynomial::create(term single_term)
{
    polynomial p;
    p.terms.push_back(single_term);
    return p;
}

polynomial polynomial::roots(variable var, polynomial expr)
{
    polynomial p;
    p.terms.push_back(term{1.0, {var_power{-1, 1, roots_expr{var, expr}}}});
    return p;
}

polynomial polynomial::create(
    const variable* indeterminates,
    size_t dimension,
    unsigned degree,
    variable& variable_counter,
    bool normalized
){
    polynomial p;
    unsigned term_count = 1;
    for(unsigned i = 0; i < dimension; ++i)
        term_count *= (degree+1);
    p.terms.reserve(term_count);
    for(unsigned i = 0; i < term_count; ++i)
    {
        term t;
        t.coefficient = 1;
        if(i != 0 || !normalized)
            t.mul.push_back({variable_counter++, 1});

        unsigned it = i;
        for(unsigned d = 0; d < dimension; ++d)
        {
            unsigned exponent = it % (degree+1);

            if(exponent != 0)
            {
                t.mul.push_back({
                    indeterminates[d],
                    exponent
                });
            }
            it /= degree+1;
        }
        p.terms.push_back(t);
    }

    return sort(p);
}

std::optional<term> differentiate(const term& t, variable id)
{
    term res;
    res.coefficient = 0;

    for(const var_power& p: t.mul)
    {
        // If there's a roots() expression that contains the variable we
        // want to differentiate over, just give up. I ain't dealing with
        // that shit.
        if(p.roots && depends_on_var(p.roots->expr, id))
            return {};

        if(p.id == id)
        {
            // d(c*x^n)/dx = c*n*x^(n-1)
            res.coefficient = t.coefficient * p.degree;
            var_power dp = p;
            dp.degree--;
            if(dp.degree > 0)
                res.mul.push_back(dp);
        }
        else res.mul.push_back(p);
    }

    // 0 coefficient -> the rest of the multipliers don't matter.
    if(res.coefficient == 0)
        res.mul.clear();
    return res;
}

term sort(const term& t)
{
    term res = t;
    std::stable_sort(res.mul.begin(), res.mul.end());
    return res;
}

polynomial multiply(const term& a, const term& b)
{
    term res;
    res.coefficient = a.coefficient * b.coefficient;
    res.mul = a.mul;
    res.mul.insert(res.mul.end(), b.mul.begin(), b.mul.end());

    // If 'a' and 'b' have common variables, they're merged in 'simplify'.
    return simplify(res);
}

polynomial simplify(const term& t)
{
    // If multiplier is zero, it doesn't matter what the rest is.
    if(t.coefficient == 0)
        return polynomial::zero();

    term res = t;

    // Simplify roots-expressions first so that the sorting and merging below
    // works.
    bool has_roots = false;
    for(var_power& p: res.mul)
    {
        if(p.roots.has_value())
        {
            p.roots->expr = simplify(p.roots->expr);
            has_roots = true;
        }
    }

    // Sort so that common variables are in successive entries.
    res = sort(res);
    for(size_t i = 0; i < res.mul.size();)
    {
        var_power& p = res.mul[i];
        if(p.degree == 0)
        { // 0 exponent -> always 1 -> variable does not matter.
            res.mul.erase(res.mul.begin()+i);
            continue;
        }
        if(i+1 < res.mul.size())
        {
            var_power& next = res.mul[i+1];
            if(compatible(p, next))
            { // Next has the same variable, merge.
                next.degree += p.degree;
                res.mul.erase(res.mul.begin()+i);
                continue;
            }
        }
        // If we're here, the variable didn't need to be cut.
        ++i;
    }

    if(!has_roots)
    { // Early exit, there are no roots to be resolved.
        return polynomial::create(res);
    }

    // 'res' should now have as simple as possible multiplicands. Let's resolve
    // roots now if possible.
    polynomial pres = polynomial::create(1.0);

    bool solved_root = false;
    for(size_t i = 0; i < res.mul.size();)
    {
        const var_power& p = res.mul[i];
        if(p.roots.has_value())
        {
            std::variant<polynomial, solve_failure_reason> root = try_solve_single_root(p.roots->expr, p.roots->var, p.degree);
            if(polynomial* p = std::get_if<polynomial>(&root))
            { // We're able to solve this root now!
                pres = multiply(pres, *p);
                res.mul.erase(res.mul.begin()+i);
                solved_root = true;
                continue;
            }
        }
        ++i;
    }
    if(solved_root)
        return multiply(pres, polynomial::create(res));
    return polynomial::create(res);
}

polynomial assign(const term& t, variable id, number value)
{
    term res = t;
    for(auto it = res.mul.begin(); it != res.mul.end();)
    {
        if(it->roots.has_value())
        {
            roots_expr& re = it->roots.value();
            // The variable whose roots we're looking for cannot be assigned
            // anymore, but otherwise it's fair game.
            if(re.var != id)
                re.expr = assign(re.expr, id, value);
        }
        else if(it->id == id)
        {
            res.coefficient *= pow(value, number(it->degree));
            it = res.mul.erase(it);
            continue;
        }
        ++it;
    }
    return simplify(res);
}

bool compatible(const var_power& a, const var_power& b)
{
    if(a.roots.has_value() && b.roots.has_value())
        return *a.roots == *b.roots;
    return a.id == b.id;
}

std::optional<term> try_sum(const term& a, const term& b)
{
    if(a.mul.size() != b.mul.size())
        return {};
    for(size_t i = 0; i < a.mul.size(); ++i)
    {
        if(!compatible(a.mul[i], b.mul[i]) || a.mul[i].degree != b.mul[i].degree)
            return {};
    }
    term res = a;
    res.coefficient += b.coefficient;
    return res;
}

polynomial assign(const term& t, variable id, const polynomial& equivalent)
{
    polynomial p;
    p.terms.push_back(t);
    term& res = p.terms.back();
    polynomial assignas = polynomial::create(1.0);

    for(auto it = res.mul.begin(); it != res.mul.end();)
    {
        if(it->id == id)
        {
            for(unsigned i = 0; i < it->degree; ++i)
                assignas = multiply(assignas, equivalent);
            it = res.mul.erase(it);
        }
        else ++it;
    }
    return multiply(p, assignas);
}

std::optional<polynomial> differentiate(const polynomial& p, variable id)
{
    polynomial res;
    for(const term& t: p.terms)
    {
        std::optional<term> dt = differentiate(t, id);
        if(dt.has_value())
            res.terms.push_back(*dt);
        else return {};
    }
    return simplify(res);
}

polynomial multiply(const polynomial& a, const polynomial& b)
{
    polynomial res;
    for(const term& ta: a.terms)
    for(const term& tb: b.terms)
    {
        polynomial p = multiply(ta, tb);
        res.terms.insert(res.terms.end(), p.terms.begin(), p.terms.end());
    }
    return simplify(res);
}

polynomial multiply(const polynomial& a, number v)
{
    polynomial res = a;
    for(term& t: res.terms)
        t.coefficient *= v;
    return simplify(res);
}

polynomial sum(const polynomial& a, const polynomial& b)
{
    polynomial res = a;
    res.terms.insert(res.terms.end(), b.terms.begin(), b.terms.end());
    return simplify(res);
}

polynomial sort(const polynomial& p)
{
    polynomial res = p;
    // Sort the variables in each term first.
    for(term& t: res.terms)
        t = sort(t);

    // Sort the terms by their variables.
    std::stable_sort(res.terms.begin(), res.terms.end());
    return res;
}

polynomial simplify(const polynomial& p)
{
    polynomial res;

    // Simplify all terms individually first.
    for(const term& t: p.terms)
    {
        polynomial p = simplify(t);
        res.terms.insert(res.terms.end(), p.terms.begin(), p.terms.end());
    }
    res = sort(res);

    // Try to merge successive terms and remove zero terms. Compatible terms
    // should be successively thanks to the sort above.
    for(unsigned i = 0; i < res.terms.size();)
    {
        term& cur = res.terms[i];
        if(cur.coefficient == 0.0)
        {
            res.terms.erase(res.terms.begin()+i);
            continue;
        }
        if(i+1 < res.terms.size())
        {
            term& next = res.terms[i+1];
            std::optional<term> sum_term = try_sum(cur, next);
            if(sum_term.has_value())
            {
                next = *sum_term;
                res.terms.erase(res.terms.begin()+i);
                continue;
            }
        }
        ++i;
    }
    // The changed sums may have changed the order here.
    return res;
}

polynomial assign(const polynomial& p, variable id, const polynomial& equivalent)
{
    polynomial res;
    for(const term& t: p.terms)
    {
        polynomial ta = assign(t, id, equivalent);
        res.terms.insert(res.terms.end(), ta.terms.begin(), ta.terms.end());
    }
    return simplify(res);
}

polynomial assign(const polynomial& p, variable id, number value)
{
    polynomial res;
    for(const term& t: p.terms)
    {
        polynomial ta = assign(t, id, value);
        res.terms.insert(res.terms.end(), ta.terms.begin(), ta.terms.end());
    }
    return simplify(res);
}

std::optional<number> evaluate(const term& t, const std::vector<number>& variable_values)
{
    number product = t.coefficient;
    for(const var_power& vp: t.mul)
    {
        number value = 0.0;
        if(vp.roots.has_value())
        {
            // Uh oh, a roots expression... This needs symbolic evaluation and
            // can be very costly.
            polynomial assigned = vp.roots->expr;
            for(unsigned id = 0; id < variable_values.size(); ++id)
            {
                // Don't replace the variable we're solving for.
                if(id == vp.roots->var)
                    continue;
                assigned = assign(assigned, id, variable_values[id]);
            }

            assigned = simplify(assigned);
            std::variant<polynomial, solve_failure_reason> result = try_solve_single_root(assigned, vp.roots->var, vp.degree);
            if(polynomial* p = std::get_if<polynomial>(&result))
            {
                std::optional<number> opt_value = try_get_constant_value(*p);
                if(!opt_value.has_value())
                    return {};
                value = *opt_value;
            }
            else return {};
        }
        else if(vp.id < variable_values.size())
        {
            value = pow(variable_values[vp.id], vp.degree);
        }
        else return {};

        product *= value;
    }
    return product;
}

std::optional<number> evaluate(const polynomial& p, const std::vector<number>& variable_values)
{
    number sum = 0.0;
    for(const term& t: p.terms)
    {
        std::optional<number> value = evaluate(t, variable_values);
        if(!value.has_value())
            return {};
        sum += *value;
    }
    return sum;
}

std::optional<polynomial> try_factor(const polynomial& p, variable id, const polynomial& root)
{
    polynomial factored;

    if(try_get_constant_value(root) == 0)
    {
        // Special case: factor as x*(factored)
        // => Try to remove a degree of x from all terms, fail if this is not
        // possible.
        factored = p;
        for(term& t: factored.terms)
        {
            bool found_var = false;
            for(var_power& vp: t.mul)
            {
                if(vp.id == id)
                {
                    vp.degree--;
                    found_var = true;
                    break;
                }
            }
            if(!found_var)
                return {};
        }
        return simplify(factored);
    }

    // Factor as (x-r)*(factored)
    // 'Coefficients' contains coefficients by degree of the original
    // polynomial.
    std::vector<polynomial> coefficients;
    for(const term& t: p.terms)
    {
        unsigned degree = 0;
        term without_t;
        without_t.coefficient = t.coefficient;
        for(const var_power& vp: t.mul)
        {
            if(vp.id == id)
                degree = std::max(degree, vp.degree);
            else without_t.mul.push_back(vp);
        }

        if(coefficients.size() <= degree)
            coefficients.resize(degree+1);
        coefficients[degree].terms.push_back(without_t);
    }

    if(coefficients.size() <= 1)
    { // Doesn't even depend on 'id', so can't factor it out.
        return {};
    }

    // Simplify the coefficients for good measure
    for(polynomial& c: coefficients)
        c = simplify(c);

    polynomial prev;

    for(unsigned i = 0; i+1 < coefficients.size(); ++i)
    {
        rational rat;
        rat.numerator = sum(prev, multiply(coefficients[i], -1.0));
        rat.denominator = root;
        rat = simplify(rat);
        if(try_get_constant_value(rat.denominator) != 1)
            return {}; // Fractional coefficient, not supported yet.

        term base_term = term{1.0, {}};
        if(i > 0) base_term.mul.push_back(var_power{id, i});
        prev = rat.numerator;

        polynomial new_term = multiply(rat.numerator, polynomial::create(base_term));
        factored.terms.insert(factored.terms.end(), new_term.terms.begin(), new_term.terms.end());
    }

    if(prev == coefficients.back())
        return simplify(factored);
    else return {}; // Not a root for this polynomial.
}

std::optional<number> try_get_constant_value(const polynomial& p)
{
    if(p.terms.size() == 0)
        return 0;
    if(p.terms.size() == 1 && p.terms[0].mul.size() == 0)
        return p.terms[0].coefficient;
    return {};
}

bool depends_on_var(const polynomial& p, variable id)
{
    for(const term& t: p.terms)
    {
        for(const var_power& vp: t.mul)
        {
            if(vp.id == id)
                return true;
            if(vp.roots && depends_on_var(vp.roots->expr, id))
                return true;
        }
    }
    return false;
}

std::set<variable> live_variables(const polynomial& p)
{
    std::set<variable> variables;
    for(const term& t: p.terms)
    for(const var_power& vp: t.mul)
    {
        if(vp.roots.has_value())
            continue;
        variables.insert(vp.id);
    }
    return variables;
}

std::variant<polynomial, solve_failure_reason> try_solve_single_root(const polynomial& p, variable id, unsigned exponent)
{
    roots_result roots = try_find_all_roots(p, id);
    if(!roots.found_all)
        return solve_failure_reason::CANT_SOLVE;
    if(roots.roots.size() == 0)
        return solve_failure_reason::NO_ROOTS;

    bool even_exponent = (exponent&1) == 0;
    if(!even_exponent && roots.roots.size() > 1)
        return solve_failure_reason::MULTIPLE_ROOTS;

    // Raise all roots to nth power
    for(polynomial& root: roots.roots)
    {
        polynomial res = polynomial::create(1.0);
        for(int i = 0; i < exponent; ++i)
            res = multiply(res, root);
        root = res;
    }

    // Ensure roots are all same.
    for(unsigned i = 1; i < roots.roots.size(); ++i)
    {
        if(!(roots.roots[i-1] == roots.roots[i]))
            return solve_failure_reason::MULTIPLE_ROOTS;
    }
    return roots.roots[0];
}

std::optional<polynomial> try_find_any_root(const polynomial& p, variable id)
{
    // Find coefficients for each degree of 'id'.
    std::vector<polynomial> coefficients;
    for(const term& t: p.terms)
    {
        unsigned degree = 0;
        term without_t;
        without_t.coefficient = t.coefficient;
        for(const var_power& vp: t.mul)
        {
            if(vp.id == id)
                degree = std::max(degree, vp.degree);
            else without_t.mul.push_back(vp);
        }

        if(coefficients.size() <= degree)
            coefficients.resize(degree+1);
        coefficients[degree].terms.push_back(without_t);
    }
    if(coefficients.size() <= 1)
    { // Never found the variable, so there isn't really a root for it.
        return {};
    }

    if(try_get_constant_value(coefficients[0]) == 0)
    { // First coefficient is zero => zero is a root.
        return polynomial::zero();
    }

    unsigned N = coefficients.size()-1;
    bool only_last_is_nonzero = true;
    for(unsigned i = 1; i < coefficients.size()-1; ++i)
    {
        if(try_get_constant_value(coefficients[i]) != 0)
            only_last_is_nonzero = false;
    }

    { // Remove common factors from coefficients, if present.
        std::map<variable, int> common;
        for(unsigned i = 0; i < coefficients.size(); ++i)
            find_common_variables(coefficients[i], common, i != 0);
        for(unsigned i = 0; i < coefficients.size(); ++i)
            coefficients[i] = factor_common_variables(coefficients[i], common);
    }

    if(only_last_is_nonzero)
    {
        // Nth degree equation, but of the form ax^N+b=0 => trivial.
        // This also includes 1st degree equations.

        // x = Nthroot(-b/a);
        rational rat;
        rat.numerator = multiply(coefficients[0], -1.0);
        rat.denominator = coefficients.back();
        rat = simplify(rat);

        std::optional<number> a = try_get_constant_value(rat.denominator);
        if(!a.has_value())
        {
            // TODO: Return rational functions instead of polynomials. This
            // entails large structural changes.
            return {};
        }

        // 'simplify' above should've already normalized a = 1, so no need to
        // divide by it.
        return try_get_nth_root(rat.numerator, N);
    }
    else if(coefficients.size() == 3)
    { // Second-degree equation, solvable.
        // ax^2 + bx + c = 0
        //One root is:
        // x=(-b+sqrt(b^2-4ac))/(2a)
        polynomial discriminant = sum(
            multiply(coefficients[1], coefficients[1]),
            multiply(multiply(coefficients[2], coefficients[0]), -4.0)
        );
        std::optional<polynomial> sqrt_discriminant = try_get_nth_root(discriminant, 2);
        if(!sqrt_discriminant.has_value())
            return {};

        rational rat;
        rat.numerator = sum(multiply(coefficients[1], -1.0), sqrt_discriminant.value());
        rat.denominator = multiply(coefficients[2], 2.0);
        rat = simplify(rat);

        std::optional<number> a = try_get_constant_value(rat.denominator);
        if(!a.has_value())
        {
            // TODO: Return rational functions instead of polynomials. This
            // entails large structural changes.
            return {};
        }
        return rat.numerator;
    }
    else
    {
        // Find extrema by finding roots of derivative, then search for
        // roots between zero-crossing extrema via binary search or something.
        // This only works if deriv_roots are constant values :(
        std::optional<polynomial> deriv = differentiate(p, id);
        if(!deriv.has_value())
            return {};
        roots_result deriv_roots = try_find_all_roots(*deriv, id);
        std::vector<number> extrema;
        for(const polynomial& p: deriv_roots.roots)
        {
            std::optional<number> val = try_get_constant_value(p);
            if(val.has_value())
                extrema.push_back(*val);
        }

        auto eval = [&](const polynomial& p, number value)
        {
            return try_get_constant_value(assign(p, id, value));
        };

        auto binary_search = [&](number start, number end) -> std::optional<number>
        {
            std::optional<number> start_y = eval(p, start);
            std::optional<number> end_y = eval(p, end);
            if(!start_y.has_value() || !end_y.has_value())
                return {};

            if(*start_y < *end_y)
            {
                std::swap(start_y, end_y);
                std::swap(start, end);
            }
            if(*start_y == 0) return start;
            if(*end_y == 0) return end;

            if(*start_y > 0 || *end_y < 0)
                return {};

            // Great, there should be a root here.
            number mid = 0;
            std::optional<number> mid_y;
            for(int iter = 0; iter < 64; ++iter)
            {
                number mid = (start + end) * 0.5;
                mid_y = eval(p, start);
                if(!mid_y.has_value())
                    break;
                if(*mid_y > 0) end = mid;
                if(*mid_y < 0) start = mid;
                if(*mid_y == 0) return mid;
            }
            if(mid_y.has_value() && fabs(*mid_y) < 1e-16)
                return mid;
            return {};
        };

        // Search between extrema first
        for(unsigned i = 1; i < extrema.size(); ++i)
        {
            number start = extrema[i-1];
            number end = extrema[i];
            std::optional<number> root = binary_search(start, end);
            if(root.has_value())
                return polynomial::create(*root);
        }

        // Search outside of extrema next
        auto exterior_search = [&](number start_from) -> std::optional<number>
        {
            std::optional<number> start_y = eval(p, start_from);
            std::optional<number> slope = eval(*deriv, start_from);
            if(!slope.has_value() || !start_y.has_value())
                return {};

            if(*start_y == 0)
                return start_from;

            number search_direction = sign(*slope);
            if(*start_y > 0) search_direction = -search_direction;

            // Find search interval
            number step_size = search_direction * 1.0;
            number start_x = *start_y;
            number end_x = *start_y;
            for(int iter = 0; iter < 64; ++iter, step_size *= 4)
            {
                number x = start_x + step_size;
                std::optional<number> y = eval(p, x);
                if(!y.has_value())
                    break;
                if(sign(*y) != sign(*start_y))
                { // Crossed over zero!
                    end_x = x;
                    break;
                }
                else
                {
                    start_x = x;
                    start_y = y;
                }
            }

            return binary_search(start_x, end_x);
        };

        std::optional<number> pos = exterior_search(extrema.size() == 0 ? 0 : extrema[0]-1.0);
        if(pos.has_value()) return polynomial::create(*pos);

        pos = exterior_search(extrema.size() == 0 ? 0 : extrema.back()+1.0);
        if(pos.has_value()) return polynomial::create(*pos);
        return {};
    }
}

std::optional<polynomial> try_get_nth_root(const polynomial& constant, unsigned N)
{
    if(N == 1) return constant;

    bool even_degree = (N&1) == 0;
    std::optional<number> b = try_get_constant_value(constant);
    if(b == 0) return polynomial::zero();
    else if(b.has_value())
    { // Constant value.
        if(even_degree && *b < 0)
            return {}; // Complex roots
        return polynomial::create(sign(*b) * pow(fabs(*b), 1.0/N));
    }

    // Time to get desperate. If we can factor 'b' by any of its variables
    // such that all roots are equal and the number of roots is a multiple
    // of 'N', it's still solvable.
    // TODO: This is good at exploding the amount of work and stalling the
    // program.
    /*
    std::set<variable> live = live_variables(constant);
    for(variable v: live)
    {
        roots_result roots = try_find_all_roots(constant, v);
        if(!roots.found_all)
            continue;

        std::optional<polynomial> residual = try_get_nth_root(roots.residual, N);
        if(!residual.has_value())
            continue;

        bool equal = true;
        for(unsigned i = 1; i < roots.roots.size(); ++i)
            if(!(roots.roots[i-1] == roots.roots[i]))
                equal = false;
        if(!equal)
            continue;

        if(roots.roots.size() > 0 && (roots.roots.size()%N) == 0)
        {
            // Compatible!
            polynomial p = *residual;
            polynomial root_expr = sum(
                polynomial::create(term{1.0, {var_power{v, 1}}}),
                multiply(roots.roots[0], -1.0)
            );
            int exp = roots.roots.size() / N;
            for(int i = 0; i < exp; ++i)
                p = multiply(p, root_expr);
            return p;
        }
    }
    */

    return {};
}

roots_result try_find_all_roots(const polynomial& p, variable id)
{
    roots_result result;
    result.found_all = false;
    result.roots = {};
    result.residual = p;
    while(depends_on_var(result.residual, id))
    {
        std::optional<polynomial> root = try_find_any_root(result.residual, id);
        if(!root.has_value())
            return result;
        std::optional<polynomial> next_residual = try_factor(result.residual, id, *root);
        if(!next_residual.has_value())
            return result;
        result.roots.push_back(*root);
        result.residual = *next_residual;
    }
    std::sort(result.roots.begin(), result.roots.end());
    result.found_all = true;
    return result;
}

void find_common_variables(const polynomial& p, std::map<variable, int>& common, bool continue_from_previous)
{
    if(!continue_from_previous)
        common.clear();

    bool first = !continue_from_previous;
    for(const term& t: p.terms)
    {
        if(first)
        {
            for(const var_power& vp: t.mul)
            {
                if(vp.id >= 0)
                    common[vp.id] = vp.degree;
            }
            first = false;
        }
        else
        {
            for(auto it = common.begin(); it != common.end();)
            {
                bool found = false;
                for(const var_power& vp: t.mul)
                {
                    if(vp.id == it->first)
                    {
                        int& degree = common[vp.id];
                        degree = degree < vp.degree ? degree : vp.degree;
                        found = true;
                    }
                }
                if(!found) it = common.erase(it);
                else ++it;
            }
        }
    }
}

polynomial factor_common_variables(const polynomial& p, const std::map<variable, int>& common)
{
    polynomial res = p;
    for(term& t: res.terms)
    {
        for(var_power& vp: t.mul)
        {
            if(common.count(vp.id))
                vp.degree -= common.at(vp.id);
        }
    }
    return simplify(res);
}

bool operator<(const indeterminate_group& a, const indeterminate_group& b)
{
    if(a.indeterminates.size() < b.indeterminates.size())
        return true;
    if(b.indeterminates.size() < a.indeterminates.size())
        return false;
    for(int i = 0; i < a.indeterminates.size(); ++i)
    {
        const var_power& ma = a.indeterminates[i];
        const var_power& mb = b.indeterminates[i];
        if(ma < mb)
            return true;
        else if(mb < ma)
            return false;
    }
    return false;
}

std::map<indeterminate_group, polynomial> group_by_indeterminates(
    const polynomial& p,
    const variable* indeterminates,
    size_t indeterminate_count
){
    std::map<indeterminate_group, polynomial> groups;
    for(const term& t: p.terms)
    {
        indeterminate_group group;
        term r = t;

        for(size_t i = 0; i < indeterminate_count; ++i)
        for(auto it = r.mul.begin(); it != r.mul.end();)
        {
            if(it->id == indeterminates[i])
            {
                group.indeterminates.push_back(*it);
                it = r.mul.erase(it);
            }
            else ++it;
        }

        groups[group].terms.push_back(r);
    }

    for(auto it = groups.begin(); it != groups.end();)
    {
        it->second = simplify(it->second);
        if(it->second.terms.empty())
            it = groups.erase(it);
        else ++it;
    }
    return groups;
}

polynomial get_zero_polynomial(const polynomial& a, number right_side)
{
    polynomial res = a;
    res.terms.push_back(term{-right_side, {}});
    return res;
}

bool pin(
    const polynomial& zero,
    const polynomial& nonzero,
    const variable* indeterminates,
    size_t indeterminate_count,
    std::vector<polynomial>& target
){
    polynomial target_nonzero = nonzero;
    std::map<indeterminate_group, polynomial> zero_groups =
        group_by_indeterminates(zero, indeterminates, indeterminate_count);

    for(auto it = zero_groups.begin(); it != zero_groups.end();)
    {
        polynomial& zero_polynomial = it->second;
        std::optional<number> value = try_get_constant_value(zero_polynomial);
        if(value.has_value())
        {
            // If it's already zero, there's nothing to be done here.
            if(*value == 0)
            {
                ++it;
                continue;
            }
            else
            {
                // Can't set a constant factor to zero in any way.
                printf("Constraint is unachievable due to a conflicting constant.\n");
                return false;
            }
        }

        std::set<variable> candidates = live_variables(zero_polynomial);

        // Find simplest equivalent expression.
        variable best_id = -1;
        solve_failure_reason best_failure = solve_failure_reason::NO_ROOTS;
        std::optional<polynomial> best_equivalent;
        for(variable id: candidates)
        {
            roots_result roots = try_find_all_roots(zero_polynomial, id);

            if(!roots.found_all)
            { // Didn't find all roots, so this _could_ be available later.
                solve_failure_reason reason = solve_failure_reason::CANT_SOLVE;
                if(!best_equivalent.has_value() && reason > best_failure)
                {
                    best_id = id;
                    best_failure = reason;
                    // fallthrough intentional - we still want to try the roots
                    // that were found, if any.
                }
            }

            // Remove roots that are incompatible with 'nonzero'.
            for(auto it = roots.roots.begin(); it != roots.roots.end();)
            {
                polynomial condition = assign(target_nonzero, id, *it);
                if(try_get_constant_value(condition) == 0)
                    it = roots.roots.erase(it);
                else ++it;
            }

            if(roots.roots.size() == 0)
            { // No suitable roots available.
                continue;
            }

            // Ensure roots are all same.
            bool all_equal = true;
            for(unsigned i = 1; i < roots.roots.size(); ++i)
            {
                if(!(roots.roots[i-1] == roots.roots[i]))
                    all_equal = false;
            }

            if(!all_equal)
            { // Multiple roots, so we won't decide which one to take here.
                solve_failure_reason reason = solve_failure_reason::MULTIPLE_ROOTS;
                if(!best_equivalent.has_value() && reason > best_failure)
                {
                    best_id = id;
                    best_failure = reason;
                    continue;
                }
            }

            polynomial& root = roots.roots[0];

            if(!best_equivalent.has_value() || root < *best_equivalent)
            {
                best_equivalent = std::move(root);
                best_id = id;
            }
        }

        if(!best_equivalent.has_value())
        { // No good equivalent polynomial found, so resort to roots().
            // Literally no roots found or zero candidates, this sucks.
            if(best_failure <= solve_failure_reason::INFINITE_ROOTS || best_id < 0)
            {
                if(best_failure == solve_failure_reason::NO_ROOTS)
                    printf("Constraint is unachievable; requires roots from equation which has none.\n");
                if(best_failure == solve_failure_reason::INFINITE_ROOTS)
                    printf("Constraint is unachievable; has infinitely many options.\n");
                return false;
            }
            best_equivalent = polynomial::roots(best_id, zero_polynomial);
        }

        it = zero_groups.erase(it);

        // Replace found equivalencies
        for(auto& pair: zero_groups)
            pair.second = assign(pair.second, best_id, *best_equivalent);
        for(polynomial& p: target)
            p = assign(p, best_id, *best_equivalent);
        target_nonzero = assign(target_nonzero, best_id, *best_equivalent);
    }

    return true;
}
