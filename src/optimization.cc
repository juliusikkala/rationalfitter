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

std::optional<number> l2_loss(
    const rational& r,
    const fit_params& params,
    const std::vector<number>& coefficients
){
    number sum_loss = 0;
    size_t evaluated_points = 0;

    std::vector<number> variable_values;
    for(size_t j = 0; j < params.data_size; ++j)
    {
        variable_values = coefficients;
        // Fetch values for defined variables
        for(const auto& pair: params.data)
        {
            if(pair.first >= variable_values.size())
                variable_values.resize(pair.first+1);
            variable_values[pair.first] = pair.second[j];
        }

        // Then, evaluate groups.
        std::optional<number> value = evaluate(r, variable_values);
        if(!value.has_value())
            continue;

        evaluated_points++;

        number right_side = params.data.at(params.right_side_variable)[j];
        number delta = *value - right_side;
        sum_loss += delta * delta;
    }

    if(evaluated_points == 0)
        return {};
    return sum_loss/evaluated_points;
}

std::optional<rational> fit(const rational& func, const fit_params& params)
{
    std::set<variable> live_vars = live_variables(func);
    std::vector<variable> fit_parameters;
    for(variable var: live_vars)
    {
        if(var != params.right_side_variable && params.data.count(var) == 0)
            fit_parameters.push_back(var);
    }

    matrix coefficients = matrix::zeroes(1, fit_parameters.size());
    for(size_t i = 0; i < fit_parameters.size(); ++i)
    {
        variable var = fit_parameters[i];
        if(params.initial_guesses.count(var) != 0)
            coefficients(0, i) = (double)params.initial_guesses.at(var);
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
        std::vector<number> variable_values;
        for(size_t i = 0; i < fit_parameters.size(); ++i)
        {
            variable var = fit_parameters[i];
            if(var >= variable_values.size())
                variable_values.resize(var+1);
            variable_values[var] = coefficients(0,i);
        }
        std::optional<number> loss = l2_loss(r, params, variable_values);
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
        nt.mul.push_back(var_power{params.right_side_variable, 1});
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

        std::vector<number> variable_values;
        std::vector<number> group_values;
        for(size_t j = 0; j < params.data_size; ++j)
        {
            // Fetch values for defined variables
            for(const auto& pair: params.data)
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
                std::optional<number> value = evaluate(pair.second, variable_values);
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
            b.values.push_back((double)-group_values[0]);
            for(size_t i = 1; i < group_values.size(); ++i)
                A.values.push_back((double)group_values[i]);
        }

        A.w = groups.size()-1;
        A.h = b.values.size();
        b.w = 1;
        b.h = b.values.size();

        coefficients = least_squares(A, b);
        printf("Least squares loss: %f\n", (double)get_loss(func, coefficients));

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

    number best_loss = get_loss(func, coefficients);
    matrix best_coefficients = coefficients;
    bool converged = false;
    for(unsigned iter = 0; iter < params.nlls_max_iterations && !converged; ++iter)
    {
        matrix J;
        matrix r;

        std::vector<number> variable_values;
        std::vector<number> derivative_values(fit_parameters.size());
        for(size_t i = 0; i < fit_parameters.size(); ++i)
        {
            variable var = fit_parameters[i];
            if(var >= variable_values.size())
                variable_values.resize(var+1);
            variable_values[var] = coefficients(0,i);
        }
        for(size_t j = 0; j < params.data_size; ++j)
        {
            // Fetch values for defined variables
            for(const auto& pair: params.data)
            {
                if(pair.first >= variable_values.size())
                    variable_values.resize(pair.first+1);
                variable_values[pair.first] = pair.second[j];
            }

            std::optional<number> func_value = evaluate(func, variable_values);
            if(!func_value.has_value())
                continue;

            bool cant_eval = false;
            for(size_t i = 0; i < fit_parameters.size(); ++i)
            {
                std::optional<number> deriv_value = evaluate(derivatives[i], variable_values);
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

            for(size_t i = 1; i < derivative_values.size(); ++i)
                J.values.push_back((double)derivative_values[i]);
            r.values.push_back((double)(params.data.at(params.right_side_variable)[j] - *func_value));
        }

        J.w = fit_parameters.size();
        J.h = r.values.size();
        r.w = 1;
        r.h = r.values.size();

        matrix delta_coefficients = least_squares(J, r);

        coefficients = add(coefficients, mul(delta_coefficients, params.nlls_step_size)).value();

        // Just to print the loss
        number loss = get_loss(func, coefficients);
        printf("Iteration %u loss: %f\n", iter, (double)loss);
        if(loss < best_loss)
        {
            best_coefficients = coefficients;
            // Convergence condition
            if(loss == 0.0 || (best_loss - loss) / best_loss < params.nlls_convergence_limit * params.nlls_step_size)
                converged = true;
            best_loss = loss;
        }
        else converged = true;
    }
    if(!converged)
        return {};
    return assign_coefficients(func, best_coefficients);
}

std::optional<rational> fit_eliminate_variable(const rational& func, const fit_params& params, std::vector<variable>& eliminated_variables){
    rational fit_src = func;
    for(variable var: eliminated_variables)
        fit_src = assign(fit_src, var, 0.0);

    std::set<variable> live_vars = live_variables(fit_src);
    if(live_vars.size() == 0)
        return fit_src;
    std::vector<variable> live_vars_vec;
    for(variable var: live_vars)
        live_vars_vec.push_back(var);

    variable best_loss_var = 0;
    number best_loss = -1.0;
    std::optional<rational> best_result;
#pragma omp parallel for
    for(size_t i = 0; i < live_vars_vec.size(); ++i)
    {
        variable var = live_vars_vec[i];
        rational without_var = assign(fit_src, var, 0.0);
        std::optional<rational> result = fit(without_var, params);
        if(!result.has_value())
            continue;

        std::optional<number> loss = l2_loss(*result, params);
        if(!loss.has_value())
            continue;

#pragma omp critical
        {
            if(*loss < best_loss || best_loss < 0)
            {
                best_loss_var = var;
                best_loss = *loss;
                best_result = result;
            }
        }
    }
    eliminated_variables.push_back(best_loss_var);
    return best_result;
}
