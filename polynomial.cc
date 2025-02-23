#include "polynomial.hh"
#include "math.hh"
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
            return true;
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
            return true;
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

polynomial polynomial::create(double value)
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
    variable& variable_counter
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

    return p;
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
    std::sort(res.mul.begin(), res.mul.end());
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
                continue;
            }
        }
        ++i;
    }
    return multiply(pres, polynomial::create(res));
}

polynomial assign(const term& t, variable id, double value)
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
            res.coefficient *= ipow(value, it->degree);
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
    std::sort(res.terms.begin(), res.terms.end());
    return res;
}

polynomial simplify(const polynomial& p, double zero_epsilon)
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
        if(fabs(cur.coefficient) <= zero_epsilon)
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

polynomial assign(const polynomial& p, variable id, double value)
{
    polynomial res;
    for(const term& t: p.terms)
    {
        polynomial ta = assign(t, id, value);
        res.terms.insert(res.terms.end(), ta.terms.begin(), ta.terms.end());
    }
    return simplify(res);
}

std::optional<double> try_get_constant_polynomial_value(const polynomial& p)
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
        variables.insert(vp.id);
    return variables;
}

std::variant<polynomial, solve_failure_reason> try_solve_single_root(const polynomial& p, variable id, unsigned exponent)
{
    unsigned max_degree = 0;
    unsigned min_degree = UINT32_MAX;
    for(const term& t: p.terms)
    {
        for(const var_power& vp: t.mul)
        {
            if(vp.id == id)
            {
                max_degree = max_degree > vp.degree ? max_degree : vp.degree;
                min_degree = min_degree < vp.degree ? min_degree : vp.degree;
            }
        }
    }

    // Never found the variable. Any value is possible, so it's kinda useless.
    if(min_degree > max_degree) return solve_failure_reason::INFINITE_ROOTS;
    if(max_degree == min_degree)
    { // Simple to solve.
        polynomial num;
        polynomial denom;
        for(const term& t: p.terms)
        {
            term filtered_t = t;
            bool has_var = false;
            for(auto it = filtered_t.mul.begin(); it != filtered_t.mul.end();)
            {
                if(it->id == id)
                {
                    has_var = true;
                    it = filtered_t.mul.erase(it);
                }
                else ++it;
            }
            if(has_var)
            {
                denom.terms.push_back(filtered_t);
            }
            else
            {
                filtered_t.coefficient = -filtered_t.coefficient;
                num.terms.push_back(filtered_t);
            }
        }

        num = simplify(num);
        denom = simplify(denom);

        std::optional<double> denom_value = try_get_constant_polynomial_value(denom);
        if(!denom_value.has_value())
        {
            // TODO: Rational functions support! Necessary! To make this work!
            return solve_failure_reason::CANT_REPRESENT;
        }

        double inv_denom = 1.0 / *denom_value;
        for(term& t: num.terms)
            t.coefficient *= inv_denom;

        std::optional<double> num_value = try_get_constant_polynomial_value(num);

        // If the root is zero, it's always the same, no questions asked.
        if(num_value.has_value() && *num_value == 0)
            return polynomial::zero();

        // If the degree is even, there's two roots unless the exponent is also
        // even.
        bool even_degree = (max_degree & 1) == 0;
        bool even_exponent = (exponent & 1) == 0;
        if(even_degree && !even_exponent)
            return solve_failure_reason::MULTIPLE_ROOTS;

        if(num_value.has_value())
        { // Constant value.
            if(even_degree && *num_value < 0)
                return solve_failure_reason::NO_ROOTS;
            if(max_degree == 1)
            { // Fast track
                return polynomial::create(ipow(*num_value, exponent));
            }
            double s = even_exponent ? 1 : sign(*num_value);
            double value = s * pow(fabs(*num_value), double(exponent)/max_degree);
            return polynomial::create(value);
        }

        if(max_degree == 1)
        {
            if(exponent == 1)
            { // Fast track
                return num;
            }
            polynomial res = polynomial::create(1.0);
            for(int i = 0; i < exponent; ++i)
                res = multiply(res, num);
            return res;
        }

        // Can't represent nth roots containing variables yet.
        return solve_failure_reason::CANT_REPRESENT;
    }
    // Has multiple roots or can't be represented, so that sucks.
    else return solve_failure_reason::CANT_SOLVE;
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

bool pin(
    const polynomial& zero,
    const variable* indeterminates,
    size_t indeterminate_count,
    std::vector<polynomial>& target
){
    std::map<indeterminate_group, polynomial> zero_groups = 
        group_by_indeterminates(zero, indeterminates, indeterminate_count);

    for(auto it = zero_groups.begin(); it != zero_groups.end();)
    {
        polynomial& zero_polynomial = it->second;
        std::optional<double> value = try_get_constant_polynomial_value(zero_polynomial);
        if(value.has_value())
        {
            // If it's already zero, there's nothing to be done here.
            if(*value == 0)
                continue;
            else
            {
                // Can't set a constant factor to zero in any way.
                return false;
            }
        }

        std::set<variable> candidates;
        for(term& t: zero_polynomial.terms)
        for(var_power& vp: t.mul)
        {
            // roots-expressions are treated as constants.
            if(vp.roots.has_value())
                continue;
            candidates.insert(vp.id);
        }

        // Find simplest equivalent expression.
        variable best_id = -1;
        solve_failure_reason best_failure = solve_failure_reason::NO_ROOTS;
        std::optional<polynomial> best_equivalent;
        for(variable id: candidates)
        {
            std::variant<polynomial, solve_failure_reason> root = try_solve_single_root(zero_polynomial, id, 1);
            if(polynomial* p = std::get_if<polynomial>(&root))
            {
                if(!best_equivalent.has_value() || *p < *best_equivalent)
                {
                    best_equivalent = std::move(*p);
                    best_id = id;
                }
            }
            else if(solve_failure_reason* r = std::get_if<solve_failure_reason>(&root))
            {
                if(!best_equivalent.has_value() && *r > best_failure)
                    best_id = id;
            }
        }

        if(!best_equivalent.has_value())
        { // No good equivalent polynomial found, so resort to roots().
            // Literally no roots found or zero candidates, this sucks.
            if(best_failure <= solve_failure_reason::INFINITE_ROOTS || best_id < 0)
                return false;
            best_equivalent = polynomial::roots(best_id, zero_polynomial);
        }

        it = zero_groups.erase(it);

        // Replace found equivalencies
        for(auto& pair: zero_groups)
            pair.second = assign(pair.second, best_id, *best_equivalent);
        for(polynomial& p: target)
            p = assign(p, best_id, *best_equivalent);
    }

    return true;
}
