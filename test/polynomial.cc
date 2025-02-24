#include "polynomial.hh"
#include "test.hh"

int main()
{
    { // Test everything that can be done on a single var_power.
        var_power a = {0, 5};
        var_power b = {0, 4};
        var_power c = {2, 0};
        CASE(a<c && !(c<a));
        CASE(b<a && !(a<b));
        CASE(b<c && !(c<b));
        CASE(a==a);
        CASE(b==b);
        CASE(c==c);

        CASE(compatible(a, b));
        CASE(!compatible(a, c));
    }

    { // Test terms
        term a = term{2.0, {var_power{0, 2}, var_power{1,2}}};
        term b = term{-3.0, {var_power{1, 4}}};
        // Comparison
        CASE(b<a && !(a<b));
        CASE(a==a);
        CASE(b==b);
        CASE(!(a==b));

        // Diffentiation
        std::optional<term> da = differentiate(a, 0);
        CASE(da.has_value() && *da == term{4.0f, {var_power{0,1}, var_power{1,2}}});
        std::optional<term> db = differentiate(b, 0);
        CASE(db.has_value() && *db == term{0.0f, {}});

        // Sorting
        CASE(sort(term{2.0, {var_power{1, 2}, var_power{0,2}}}) == a);

        // Multiplication
        CASE(multiply(a, b) == polynomial::create(term{-6.0, {var_power{0,2}, var_power{1,6}}}));
        CASE(multiply(a, term{1, {}}) == polynomial::create(a));
        CASE(multiply(a, term{0, {}}) == polynomial::zero());
    }
    FINISH;
}
