#include "rational.hh"
#include <cmath>

std::optional<rational> differentiate(const rational& r, variable id)
{
    std::optional<polynomial> dnum = differentiate(r.numerator, id);
    std::optional<double> denom_value = try_get_constant_value(r.denominator);
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

    // Normalize to a Padé approximant if possible.
    double factor = 0;
    for(term& t: res.denominator.terms)
    {
        if(t.mul.size() == 0)
            factor += t.coefficient;
    }

    if(fabs(factor) < 1e-10)
    {
        // Yikes; looks like we don't have a constant in the denominator.
        // This does luckily present an opportunity: it's possible that there
        // is a common variable in all terms that we can factor out.
        std::map<variable, int /*min_degree*/> common;

        bool first_check = true;
        auto check_term = [&](term& t){
            if(first_check)
            {
                for(var_power& vp: t.mul)
                {
                    if(vp.id >= 0)
                        common[vp.id] = vp.degree;
                }
                first_check = false;
            }
            else
            {
                for(auto it = common.begin(); it != common.end();)
                {
                    bool found = false;
                    for(var_power& vp: t.mul)
                    {
                        if(vp.id == it->first)
                        {
                            int degree = common[vp.id];
                            degree = degree < vp.degree ? degree : vp.degree;
                            found = true;
                        }
                    }
                    if(!found) it = common.erase(it);
                    else ++it;
                }
            }
        };
        auto factor_term = [&](term& t){
            for(var_power& vp: t.mul)
            {
                if(common.count(vp.id))
                    vp.degree -= common[vp.id];
            }
        };
        for(term& t: res.numerator.terms) check_term(t);
        for(term& t: res.denominator.terms) check_term(t);
        for(term& t: res.numerator.terms) factor_term(t);
        for(term& t: res.denominator.terms) factor_term(t);

        // Factoring out the common variables may have revealed a new constant,
        // so let's just re-simplify the whole thing.
        return common.size() > 0 ? simplify(res) : res;
    }
    else
    {
        for(term& t: res.numerator.terms)
            t.coefficient /= factor;
        for(term& t: res.denominator.terms)
            t.coefficient /= factor;
        return res;
    }
}

rational assign(const rational& r, variable id, double value)
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

polynomial get_zero_polynomial(const rational& r, double right_side)
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
