#include "optimization.hh"
#include <cstdio>

matrix least_squares(const matrix& A, const matrix& b)
{
    auto[R, v] = lstsq_r_decompose(A, b);

    //printf("A = %s\n", to_string(A, "    ").c_str());
    //printf("b = %s\n", to_string(b, "    ").c_str());
    //printf("R = %s\n", to_string(R, "    ").c_str());
    //printf("v = %s\n", to_string(v, "    ").c_str());

    matrix x = matrix::zeroes(1, A.w);

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
    size_t data_size,
    const std::vector<double>& coefficients
){
    double sum_loss = 0;
    size_t evaluated_points = 0;

    std::vector<double> variable_values;
    for(size_t j = 0; j < data_size; ++j)
    {
        variable_values = coefficients;
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
    const rational& func,
    variable right_side_variable,
    const std::map<variable, const double*>& data,
    size_t data_size,
    const std::map<variable, double>& initial_guesses,
    double nlls_step_size,
    unsigned nlls_max_iterations,
    double nlls_convergence_limit
){
    std::set<variable> live_vars = live_variables(func);
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

    auto assign_coefficients = [&](const rational& r, const matrix& coefficients){
        rational rat = r;
        for(size_t i = 0; i < fit_parameters.size(); ++i)
            rat = assign(rat, fit_parameters[i], coefficients(0,i));
        return rat;
    };

    auto get_loss = [&](const rational& r, const matrix& coefficients)
    {
        std::vector<double> variable_values;
        for(size_t i = 0; i < fit_parameters.size(); ++i)
        {
            variable var = fit_parameters[i];
            if(var >= variable_values.size())
                variable_values.resize(var+1);
            variable_values[var] = coefficients(0,i);
        }
        std::optional<double> loss = l2_loss(r, right_side_variable, data, data_size, variable_values);
        return loss.has_value() ? *loss : -1.0f;
    };

    // Try to start with linear least squares. This may fail, and that's OK.
    // Turn the rational into something fittable with least squares:
    // P(x)/Q(x) = R
    // => P(X)-R*Q(x) = 0
    polynomial linear = func.numerator;

    for(const term& t: func.denominator.terms)
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
        printf("Starting least squares\n");
        matrix A;
        matrix b;

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

        coefficients = least_squares(A, b);
        printf("Least squares loss: %f\n", get_loss(func, coefficients));

        std::optional<rational> res = assign_coefficients(func, coefficients);
        if(try_get_constant_value(func.denominator) == 1.0)
        { // If this was just a polynomial, we're done.
            return res;
        }

        // Otherwise, we continue on to non-linear optimization with the above
        // coefficients as a guess.
    }

    // This implementation uses the Gauss-Newton algorithm, since we're usually
    // able to get derivatives from the rational.
    // https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm

    printf("Starting Gauss-Newton NLLS\n");

    // Compute derivatives along parameters.
    std::vector<rational> derivatives(fit_parameters.size());
    for(size_t i = 0; i < fit_parameters.size(); ++i)
    {
        std::optional<rational> deriv = differentiate(func, fit_parameters[i]);
        if(!deriv.has_value())
        { // Can't differentiate, so this is a lost cause.
            return {};
        }
        derivatives[i] = *deriv;
    }

    double best_loss = get_loss(func, coefficients);
    matrix best_coefficients = coefficients;
    bool converged = false;
    for(unsigned iter = 0; iter < nlls_max_iterations && !converged; ++iter)
    {
        matrix J;
        matrix r;

        std::vector<double> variable_values;
        std::vector<double> derivative_values(fit_parameters.size());
        for(size_t i = 0; i < fit_parameters.size(); ++i)
        {
            variable var = fit_parameters[i];
            if(var >= variable_values.size())
                variable_values.resize(var+1);
            variable_values[var] = coefficients(0,i);
        }
        for(size_t j = 0; j < data_size; ++j)
        {
            // Fetch values for defined variables
            for(const auto& pair: data)
            {
                if(pair.first >= variable_values.size())
                    variable_values.resize(pair.first+1);
                variable_values[pair.first] = pair.second[j];
            }

            std::optional<double> func_value = evaluate(func, variable_values);
            if(!func_value.has_value())
                continue;

            bool cant_eval = false;
            for(size_t i = 0; i < fit_parameters.size(); ++i)
            {
                std::optional<double> deriv_value = evaluate(derivatives[i], variable_values);
                if(!deriv_value.has_value())
                {
                    cant_eval = true;
                    break;
                }
                derivative_values[i] = *deriv_value;
            }
            // Skip this datapoint if we can't eval something, it may just be
            // an extremity.
            if(cant_eval) continue;

            J.values.insert(J.values.end(), derivative_values.begin(), derivative_values.end());
            r.values.push_back(data.at(right_side_variable)[j] - *func_value);
        }

        J.w = fit_parameters.size();
        J.h = r.values.size();
        r.w = 1;
        r.h = r.values.size();

        matrix delta_coefficients = least_squares(J, r);

        coefficients = add(coefficients, mul(delta_coefficients, nlls_step_size)).value();

        // Just to print the loss
        double loss = get_loss(func, coefficients);
        printf("Iteration %u loss: %f\n", iter, loss);
        if(loss < best_loss)
        {
            best_coefficients = coefficients;
            // Convergence condition
            if(loss == 0.0 || (best_loss - loss) / best_loss < nlls_convergence_limit * nlls_step_size)
                converged = true;
            best_loss = loss;
        }
        else converged = true;
    }
    if(!converged)
        return {};
    return assign_coefficients(func, best_coefficients);
}
