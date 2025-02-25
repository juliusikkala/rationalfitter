#include "optimization.hh"
#include <cstdio>

matrix least_squares(const matrix& A, const matrix& b)
{
    printf("Running QR decomposition\n");
    auto[R, v] = lstsq_r_decompose(A, b);

    //printf("A = %s\n", to_string(A, "    ").c_str());
    //printf("b = %s\n", to_string(b, "    ").c_str());
    //printf("R = %s\n", to_string(R, "    ").c_str());
    //printf("v = %s\n", to_string(v, "    ").c_str());

    matrix x = matrix::zeroes(1, A.w);

    printf("Back-solving\n");
    for(int i = R.w-1; i >= 0; --i)
    {
        double right_side = v(0, i);
        for(int j = i+1; j < R.w; ++j)
            right_side -= R(j,i) * x(0, j);
        x(0, i) = right_side / R(i,i);
    }
    return x;
}

std::optional<double> l2_loss(
    const rational& r,
    variable right_side_variable,
    const std::map<variable, const double*>& data,
    size_t data_size
){
    double sum_loss = 0;
    size_t evaluated_points = 0;

    std::vector<double> variable_values;
    for(size_t j = 0; j < data_size; ++j)
    {
        // Fetch values for defined variables
        for(const auto& pair: data)
        {
            if(pair.first >= variable_values.size())
                variable_values.resize(pair.first+1);
            variable_values[pair.first] = pair.second[j];
        }

        // Then, evaluate groups.
        std::optional<double> value = evaluate(r, variable_values);
        if(!value.has_value())
            continue;

        evaluated_points++;

        double right_side = data.at(right_side_variable)[j];
        double delta = *value - right_side;
        sum_loss += delta * delta;
    }

    if(evaluated_points == 0)
        return {};
    return sum_loss/evaluated_points;
}

std::optional<rational> fit(
    const rational& r,
    variable right_side_variable,
    const std::map<variable, const double*>& data,
    size_t data_size,
    const std::map<variable, double>& initial_guesses
){
    std::set<variable> live_vars = live_variables(r);
    std::vector<variable> fit_parameters;
    for(variable var: live_vars)
    {
        if(var != right_side_variable && data.count(var) == 0)
            fit_parameters.push_back(var);
    }

    matrix coefficients = matrix::zeroes(1, fit_parameters.size());
    for(size_t i = 0; i < fit_parameters.size(); ++i)
    {
        variable var = fit_parameters[i];
        if(initial_guesses.count(var) != 0)
            coefficients(0, i) = initial_guesses.at(var);
        else coefficients(0, i) = 1;
    }

    auto assign_coefficients = [&](){
        rational rat = r;
        for(size_t i = 0; i < fit_parameters.size(); ++i)
            rat = assign(rat, fit_parameters[i], coefficients(0,i));
        std::optional<double> loss = l2_loss(rat, right_side_variable, data, data_size);
        if(loss.has_value())
            printf("L2 loss: %f\n", *loss);
        else
            printf("Failed to compute loss for some reason\n");
        return rat;
    };

    // Try to start with linear least squares. This may fail, and that's OK.
    // Turn the rational into something fittable with least squares:
    // P(x)/Q(x) = R
    // => P(X)-R*Q(x) = 0
    polynomial linear = r.numerator;

    for(const term& t: r.denominator.terms)
    {
        term nt = t;
        nt.coefficient = -nt.coefficient;
        nt.mul.push_back(var_power{right_side_variable, 1});
        linear.terms.push_back(nt);
    }

    std::map<indeterminate_group, polynomial> groups = group_by_indeterminates(
        linear,
        fit_parameters.data(),
        fit_parameters.size()
    );
    // Make sure all groups' indeterminates are individual. That means that
    // this is a simple linear combination, which we need for linear least
    // squares.
    bool linear_combination = true;
    for(const auto& pair: groups)
    {
        if(pair.first.indeterminates.size() > 1)
            linear_combination = false;
    }
    // Ensure constant group exists.
    groups[indeterminate_group()];

    if(linear_combination && groups.size() == fit_parameters.size()+1)
    { // Great, linear least squares should work.
        printf("Trying least squares\n");
        matrix A;
        matrix b;

        printf("Collecting terms at data points\n");
        std::vector<double> variable_values;
        std::vector<double> group_values;
        for(size_t j = 0; j < data_size; ++j)
        {
            // Fetch values for defined variables
            for(const auto& pair: data)
            {
                if(pair.first >= variable_values.size())
                    variable_values.resize(pair.first+1);
                variable_values[pair.first] = pair.second[j];
            }

            // Then, evaluate groups.
            bool cant_eval = false;
            group_values.clear();
            for(const auto& pair: groups)
            {
                std::optional<double> value = evaluate(pair.second, variable_values);
                if(!value.has_value())
                {
                    cant_eval = true;
                    break;
                }
                group_values.push_back(*value);
            }
            // Skip this datapoint if we can't eval something, it may just be
            // an extremity.
            if(cant_eval) continue;
            b.values.push_back(-group_values[0]);
            A.values.insert(A.values.end(), group_values.begin()+1, group_values.end());
        }

        A.w = groups.size()-1;
        A.h = b.values.size();
        b.w = 1;
        b.h = b.values.size();

        printf("Running linear least squares\n");
        coefficients = least_squares(A, b);

        std::optional<rational> res = assign_coefficients();
        if(try_get_constant_value(r.denominator) == 1.0)
        { // If this was just a polynomial, we're done.
            return res;
        }

        // Otherwise, we continue on to non-linear optimization with the above
        // coefficients as a guess.
    }

    // TODO: Non-linear optimization if we fall through here.
    return {};
}
