#include "rational.hh"
#include <cmath>
#include <algorithm>

std::optional<rational> differentiate(const rational& r, variable id)
{
    std::optional<polynomial> dnum = differentiate(r.numerator, id);
    std::optional<number> denom_value = try_get_constant_value(r.denominator);
    if(denom_value.has_value() && *denom_value == 1)
    { // Fast track for pure polynomials.
        if(dnum.has_value())
            return rational{*dnum, r.denominator};
        else return {};
    }

    std::optional<polynomial> ddenom = differentiate(r.denominator, id);
    if(!dnum.has_value() || !ddenom.has_value())
        return {};

    rational res;
    res.numerator = r.numerator;
    for(term& t: res.numerator.terms)
        t.coefficient = -t.coefficient;

    res.numerator = sum(
        multiply(r.denominator, *dnum),
        multiply(res.numerator, *ddenom)
    );
    res.denominator = multiply(r.denominator, r.denominator);

    return simplify(res);
}

rational simplify(const rational& r)
{
    rational res;

    res.numerator = simplify(r.numerator);
    res.denominator = simplify(r.denominator);

    number denom_constant = 0;
    number denom_scale = 0;
    number num_scale = 0;
    for(term& t: res.denominator.terms)
    {
        if(t.mul.size() == 0)
            denom_constant += t.coefficient;
        denom_scale = fabs(t.coefficient) > fabs(denom_scale) ?
            t.coefficient : denom_scale;
    }
    for(term& t: res.numerator.terms)
    {
        num_scale = fabs(t.coefficient) > fabs(num_scale) ?
            t.coefficient : num_scale;
    }

    if(denom_scale != 0)
    {
        for(term& t: res.denominator.terms)
            t.coefficient /= denom_scale;
    }
    if(num_scale != 0)
    {
        for(term& t: res.numerator.terms)
            t.coefficient /= num_scale;
    }

    // Identical numerator and denominator, trivial to simplify.
    if(res.numerator == res.denominator)
        return rational{polynomial::create(num_scale/denom_scale), polynomial::create(1.0)};

    // Normalize by denominator constant.
    if(denom_constant == 0) denom_constant = 1;
    for(term& t: res.numerator.terms)
        t.coefficient *= num_scale / denom_constant;
    for(term& t: res.denominator.terms)
        t.coefficient *= denom_scale / denom_constant;

    // Otherwise, try to factorize and simplify that way.
    /*
    std::set<variable> num_live = live_variables(res.numerator);
    std::set<variable> denom_live = live_variables(res.denominator);
    std::vector<variable> live;
    std::set_intersection(
        denom_live.begin(), denom_live.end(), num_live.begin(), num_live.end(),
        std::inserter(live, live.begin())
    );

    bool factored = false;
    for(variable id: live)
    {
        roots_result num_roots = try_find_all_roots(res.numerator, id);
        for(polynomial& num_root: num_roots.roots)
        {
            std::optional<polynomial> denom_factored = try_factor(res.denominator, id, num_root);
            if(!denom_factored.has_value()) continue;
            std::optional<polynomial> num_factored = try_factor(res.numerator, id, num_root);
            if(!num_factored.has_value()) continue;
            res.denominator = *denom_factored;
            res.numerator = *num_factored;
            factored = true;
        }
        roots_result denom_roots = try_find_all_roots(res.denominator, id);
        for(polynomial& denom_root: denom_roots.roots)
        {
            std::optional<polynomial> denom_factored = try_factor(res.denominator, id, denom_root);
            if(!denom_factored.has_value()) continue;
            std::optional<polynomial> num_factored = try_factor(res.numerator, id, denom_root);
            if(!num_factored.has_value()) continue;
            res.denominator = *denom_factored;
            res.numerator = *num_factored;
            factored = true;
        }
    }

    if(factored) return simplify(res);
    */
    return res;
}

rational assign(const rational& r, variable id, number value)
{
    rational res;
    res.numerator = assign(r.numerator, id, value);
    res.denominator = assign(r.denominator, id, value);
    return simplify(res);
}

std::set<variable> live_variables(const rational& r)
{
    std::set<variable> live = live_variables(r.numerator);
    live.merge(live_variables(r.denominator));
    return live;
}

polynomial get_zero_polynomial(const rational& r, number right_side)
{
    // If the right side is NaN, this means that our denominator is zero.
    if(std::isnan((double)right_side))
    {
        return r.denominator;
    }
    else
    {
        polynomial res = r.numerator;

        for(const term& t: r.denominator.terms)
        {
            term nt = t;
            nt.coefficient = -nt.coefficient * right_side;
            res.terms.push_back(nt);
        }

        return simplify(res);
    }
}

std::optional<number> evaluate(const rational& r, const std::vector<number>& variable_values)
{
    std::optional<number> num_value = evaluate(r.numerator, variable_values);
    std::optional<number> denom_value = evaluate(r.denominator, variable_values);
    if(!num_value.has_value() || !denom_value.has_value())
        return {};
    if(*denom_value == 0)
        return {};
    return *num_value / *denom_value;
}
