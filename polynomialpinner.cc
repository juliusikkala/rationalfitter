#include "polynomial.hh"
#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <functional>
#include <variant>
#include <map>
#include <set>

const char* command_list_format = R"(Command list format:
'#' starts a comment row and is not parsed.

'polynomial <degree> <in-axes> <out-axis>' initializes a polynomial.
'pin x=<a> y=<b>' makes the polynomial equal to 'b' at 'a'.
'print [multiline] [lc]' prints the current state of the polynomial.
    'multiline' makes it split each term on a new line, and 'lc' reformats the
    polynomial as a linear combination.
'differentiate <axis>' differentiates the expression along given axis.
'let <axis>=<value>' replaces <axis> from expression using given <value>.
)";

struct context
{
    std::map<variable, std::string> var_names;
    std::map<std::string, variable> name_vars;
    std::vector<variable> axes;
    std::vector<variable> coefficients;
    polynomial p;
};

bool read_text_file(const char* path, std::string& data)
{
    FILE* f = fopen(path, "rb");

    if(!f) return false;

    fseek(f, 0, SEEK_END);
    size_t bytes = ftell(f);
    fseek(f, 0, SEEK_SET);

    data.resize(bytes);
    if(fread(data.data(), 1, bytes, f) != bytes)
    {
        data.clear();
        return false;
    }
    fclose(f);

    return true;
}

bool is_integer(double d)
{
    return d == double(int64_t(d));
}

/*
struct polynomial
{
    polynomial reassign_variable_names() const
    {
        polynomial copy = *this;
        std::map<char, char> old_to_new;
        for(coefficient& c: copy.coefficients)
        for(var_weight& w: c.sum)
            if(w.name)
                old_to_new[w.name];

        char name_counter = 'a';
        for(auto& pair: old_to_new)
        {
            if(name_counter == 'e') ++name_counter;
            pair.second = name_counter++;
        }

        for(coefficient& c: copy.coefficients)
        for(var_weight& w: c.sum)
            if(w.name)
                w.name = old_to_new[w.name];
        return copy;
    }

    void print_c()
    {
        std::set<char> vars;
        for(const coefficient& c: coefficients)
        for(const var_weight& w: c.sum)
            vars.insert(w.name);

        printf("float func(");
        std::vector<unsigned> max_degrees(output_axis, 0);
        for(const coefficient& c: coefficients)
        {
            for(unsigned i = 0; i < max_degrees.size(); ++i)
                max_degrees[i] = max_degrees[i] > c.degrees[i] ? max_degrees[i] : c.degrees[i];
        }

        bool first_param = true;
        for(unsigned i = 0; i < max_degrees.size(); ++i)
        {
            if(max_degrees[i] > 0)
            {
                if(!first_param) printf(", ");
                printf("float %c", 'x'+i);
                first_param = false;
            }
        }

        for(char var: vars)
        {
            if(var != 0)
            {
                if(!first_param) printf(", ");
                printf("float %c", var);
                first_param = false;
            }
        }
        printf(")\n{\n");
        for(unsigned i = 0; i < max_degrees.size(); ++i)
        {
            int degree = max_degrees[i];
            char name = 'x'+i;
            std::string prevname(1, name);
            for(unsigned j = 2; j <= degree; ++j)
            {
                std::string dname(1, name);
                dname += std::to_string(j);
                printf("    float %s = %s * %c;\n", dname.c_str(), prevname.c_str(), name);
                prevname = dname;
            }
        }

        printf("    float res = 0;\n");
        for(unsigned i = 0; i < coefficients.size(); ++i)
        {
            const std::vector<unsigned>& degrees = coefficients[i].degrees;

            std::string name = coefficient_name(i);
            std::string term;
            bool first_d = true;
            for(unsigned d = 0; d < degrees.size(); ++d)
            {
                char letter = 'x' + d;
                if(degrees[d] == 0)
                    continue;
                if(!first_d)
                    term += " * ";
                first_d = false;
                term += letter;
                if(degrees[d] > 1)
                    term += std::to_string(degrees[d]);
            }
            if(name[0] == '-')
            {
                printf("    res -= ");
                name.erase(name.begin());
            }
            else
            {
                printf("    res += ");
            }
            bool skip_name = name == "1" && !term.empty();
            if(!skip_name)
                printf("%s", name.c_str());
            if(term.size() > 0)
            {
                if(!skip_name)
                    printf(" * ");
                printf("%s", term.c_str());
            }
            printf(";\n");
        }
        printf("    return res;\n}\n");
    }

    void print_numpy()
    {
        std::set<char> vars;
        for(const coefficient& c: coefficients)
        for(const var_weight& w: c.sum)
            vars.insert(w.name);

        // We always add the constant offset even if it's not being used.
        vars.insert(0);

        std::string constant;
        std::vector<std::string> multipliers;
        for(char var: vars)
        {
            std::string sum;
            int sum_entries = 0;
            bool sum_has_axes = false;
            for(const coefficient& c: coefficients)
            {
                std::string term;
                bool first_d = true;
                for(unsigned d = 0; d < c.degrees.size(); ++d)
                {
                    char letter = 'x' + d;
                    if(c.degrees[d] == 0)
                        continue;
                    if(!first_d)
                        term += " * ";
                    first_d = false;
                    term += letter;
                    if(c.degrees[d] > 1)
                        term += "**"+std::to_string(c.degrees[d]);
                }
                for(const var_weight& w: c.sum)
                {
                    if(w.name == var)
                    {
                        float weight = w.weight;

                        if(sum_entries != 0)
                        {
                            sum += " ";
                            if(weight >= 0)
                                sum += "+ ";
                            else
                            {
                                sum += "- ";
                                weight = fabs(weight);
                            }
                        }
                        else if(weight < 0)
                        {
                            sum += "-";
                            weight = fabs(weight);
                        }

                        if(weight != 1 || term.empty())
                        {
                            sum += is_integer(weight) ?
                                std::to_string(int64_t(weight)) :
                                std::to_string(weight);
                            if(!term.empty())
                                sum += " * ";
                        }

                        if(!term.empty())
                            sum_has_axes = true;

                        sum += term;
                        sum_entries++;
                    }
                }
            }

            if(var == 0)
            {
                if(sum_has_axes)
                    constant = "offset = "+sum;
                else if(sum_entries != 0)
                    constant = "offset = np.ones(x.shape) * ("+sum+")";
                else
                    constant = "offset = np.zeros(x.shape)";
            }
            else
            {
                multipliers.push_back(sum + " # " + var);
            }
        }

        std::string multipliers_str;
        for(unsigned i = 0; i < multipliers.size(); ++i)
        {
            if(i != 0)
                multipliers_str += "        , ";
            else
                multipliers_str += "        ";
            multipliers_str += multipliers[i] + "\n";
        }

        std::string params;
        std::string flats;
        for(unsigned i = 0; i <= output_axis; ++i)
        {
            if(i != 0)
                params += ", ";
            params += ('X' + i);
            flats += "    ";
            flats += ('x' + i);
            flats += " = ";
            flats += ('X' + i);
            flats += ".flatten()\n";
        }

        printf(R"(import numpy as np

def fit(%s):
%s
    %s
    b = %c - offset
    A = np.array([
%s    ]).T
    coeff, r, rank, s = np.linalg.lstsq(A, b, rcond=None)
    out = (np.dot(A, coeff) + offset).reshape(X.shape)
    return (coeff, out)
)",
            params.c_str(),
            flats.c_str(),
            constant.c_str(),
            'x'+output_axis,
            multipliers_str.c_str()
        );
    }

    unsigned output_axis;
    std::vector<coefficient> coefficients;
};
*/
std::string polynomial_to_string(
    context& ctx,
    const polynomial& p,
    bool group_by_axes,
    bool multiline,
    const char* pow_symbol = "^"
);

std::string term_to_string(
    context& ctx,
    const term& t,
    const char* pow_symbol
){
    if(t.coefficient == 0)
        return "0";

    float coefficient = fabs(t.coefficient);

    std::string ret = t.coefficient < 0 ? "- " : "";
    bool show_coefficient = false;
    if(coefficient != 1 || t.mul.size() == 0)
    {
        ret += is_integer(coefficient) ? 
            std::to_string(int64_t(coefficient)) :
            std::to_string(coefficient);
        show_coefficient = true;
    }

    for(size_t i = 0; i < t.mul.size(); ++i)
    {
        if(i != 0 || show_coefficient)
            ret += " * ";
        const var_power& vp = t.mul[i];
        if(vp.roots.has_value())
        {
            ret += "roots(" + polynomial_to_string(ctx, vp.roots->expr, true, false, pow_symbol) + ", " + ctx.var_names[vp.roots->var] + ")";
        }
        else
        {
            ret += ctx.var_names[vp.id];
        }
        if(vp.degree != 1)
        {
            ret += pow_symbol;
            ret += std::to_string(vp.degree);
        }
    }

    return ret;
}

std::string polynomial_to_string(
    context& ctx,
    const polynomial& p,
    bool group_by_axes,
    bool multiline,
    const char* pow_symbol
){
    if(p.terms.size() == 0)
        return "0";

    std::vector<variable>& indeterminates = group_by_axes ? ctx.axes : ctx.coefficients;
    std::map<indeterminate_group, polynomial> groups = 
        group_by_indeterminates(p, indeterminates.data(), indeterminates.size());

    bool first_line = true;
    std::string polynomial_string;
    for(auto& [group, poly]: groups)
    {
        std::string poly_string;
        bool first = true;
        for(const term& t: poly.terms)
        {
            std::string tstr = term_to_string(ctx, t, pow_symbol);
            if(tstr.empty())
                continue;

            if(!first)
            {
                if(tstr[0] == '-')
                    poly_string += " ";
                else
                    poly_string += " + ";
            }

            poly_string += tstr;
            first = false;
        }

        term indeterminate_term;
        indeterminate_term.coefficient = 1;
        indeterminate_term.mul = group.indeterminates;
        std::string indeterminate_string = term_to_string(ctx, indeterminate_term, pow_symbol);

        bool implicit_coefficient = false;
        if(poly.terms.size() > 1)
        {
            poly_string = "(" + poly_string + ")";
        }

        std::string joined_string;
        if(indeterminate_string == "1")
            joined_string = poly_string;
        else
        {
            if(poly_string == "- 1")
                joined_string = "- " + indeterminate_string;
            else if(poly_string == "1")
                joined_string = indeterminate_string;
            else joined_string = poly_string + " * " + indeterminate_string;
        }

        if(joined_string == "0" || joined_string.empty())
            continue;

        if(!first_line)
        {
            if(joined_string[0] == '-')
                polynomial_string += multiline ? "\n" : " ";
            else polynomial_string += multiline ? "\n+ " : " + ";
        }
        polynomial_string += joined_string;
        first_line = false;
    }
    return polynomial_string;
}

void skip_whitespace(const char*& str)
{
    while(*str != 0 && strchr(" \t\n", *str)) str++;
}

void skip_spaces(const char*& str)
{
    while(*str != 0 && strchr(" =\t", *str)) str++;
}

bool read_double(const char*& str, double& d)
{
    char* out = nullptr;
    d = strtod(str, &out);
    if(*out == 0 || strchr(" =\t\n", *out))
    {
        str = out;
        return true;
    }
    else return false;
}

std::string read_token(const char*& str)
{
    std::string token;
    while(*str && strchr(" =\t\n", *str) == nullptr)
    {
        token += *str;
        str++;
    }
    return token;
}

using parameter = std::variant<std::string, double>;

int get_params(const parameter* parameters, size_t param_count)
{
    if(param_count != 0)
        return 1; // Argument count mismatch
    return 0;
}

template<typename U, typename... T>
int get_params(const parameter* parameters, size_t param_count, U& first, T&... rest)
{
    if(param_count == 0)
        return 1; // Argument count mismatch.

    if constexpr(
        std::is_same_v<U, int> ||
        std::is_same_v<U, float> ||
        std::is_same_v<U, double>
    ){
        if(const double* val = std::get_if<double>(parameters))
            first = *val;
        else return 2; // Argument type mismatch.
    }
    else if constexpr(std::is_same_v<U, const char*>)
    {
        if(const std::string* val = std::get_if<std::string>(parameters))
            first = val->c_str();
        else return 2; // Argument type mismatch.
    }
    else if constexpr(std::is_same_v<U, std::string>)
    {
        if(const std::string* val = std::get_if<std::string>(parameters))
            first = *val;
        else return 2; // Argument type mismatch.
    }
    else if constexpr(std::is_same_v<U, std::vector<std::string>>)
    {
        for(size_t i = 0; i < param_count; ++i)
        {
            if(const std::string* val = std::get_if<std::string>(parameters+i))
                first.push_back(*val);
            else return 2; // Argument type mismatch.
        }
        return get_params(nullptr, 0, rest...);
    }
    else if constexpr(std::is_same_v<U, parameter>)
    {
        first = *parameters;
    }

    return get_params(parameters+1, param_count-1, rest...);
}

template<typename... T>
int va_sizeof(T&... rest) {return sizeof...(rest);}

#define PARAMS(...) \
    { \
        int ret = get_params(parameters.data(), parameters.size(), __VA_ARGS__);\
        if(ret == 1)\
        {\
            fprintf(stderr, "Argument count mismatch: expected %d, got %lu\n", va_sizeof(__VA_ARGS__), parameters.size());\
            return false; \
        }\
        else if(ret == 2)\
        {\
            fprintf(stderr, "Argument type mismatch\n"); \
            return false; \
        }\
    }\

#define NO_PARAMS() \
    { \
        int ret = get_params(parameters.data(), parameters.size());\
        if(ret == 1)\
        {\
            fprintf(stderr, "Argument count mismatch: expected 0, got %lu\n", parameters.size());\
            return false; \
        }\
    }\

using command_handler = std::function<bool(context&, const std::vector<parameter>& parameters)>;

const std::unordered_map<std::string, command_handler> command_handlers = {
    {"print", [](context& ctx, const std::vector<parameter>& parameters)->bool
    {
        bool linear_combination = false;
        bool multiline = false;
        bool numpy = false;
        bool c_code = false;
        for(const parameter& p: parameters)
        {
            if(const std::string* val = std::get_if<std::string>(&p))
            {
                if(*val == "lc") linear_combination = true;
                else if(*val == "multiline") multiline = true;
                else if(*val == "numpy") numpy = true;
                else if(*val == "c") c_code = true;
            }
            else
            {
                printf("Unrecognized argument\n");
                return false;
            }
        }
        /*
        if(c_code)
            p.print_c();
        else if(numpy)
            p.print_numpy();
            */

        std::string str = polynomial_to_string(
            ctx,
            ctx.p,
            !linear_combination,
            multiline,
            "^"
        );
        printf("%s\n", str.c_str());
        return true;
    }},
    {"polynomial", [](context& ctx, const std::vector<parameter>& parameters)->bool
    {
        int degree;
        std::vector<std::string> axis_names;
        PARAMS(degree, axis_names);
        if(degree < 0)
        {
            printf("Degree must be a positive integer!\n");
            return false;
        }

        int dimension = axis_names.size()-1;
        if(dimension < 0 || dimension > 4)
        {
            printf("Currently, only dimensions [0, 4] are supported.\n");
            return false;
        }
        ctx.axes.clear();
        ctx.coefficients.clear();
        ctx.var_names.clear();
        ctx.name_vars.clear();
        for(size_t i = 0; i < axis_names.size(); ++i)
        {
            variable id = (1<<20)+i;
            ctx.axes.push_back(id);
            ctx.var_names[id] = axis_names[i];
            if(ctx.name_vars.count(axis_names[i]))
            {
                printf("Axis names must be unique.\n");
                return false;
            }
            ctx.name_vars[axis_names[i]] = id;
        }

        variable var_counter = 0;
        ctx.p = polynomial::create(ctx.axes.data(), dimension, degree, var_counter);
        bool alphabet_names = var_counter < 'z'-'a'-axis_names.size();
        char alphabet_counter = 'a';
        for(variable i = 0; i < var_counter; ++i)
        {
            std::string name;
            if(alphabet_names)
            {
                // Skip 'e' for desmos compatibility :D
                while(ctx.name_vars.count(std::string(1, alphabet_counter)) || alphabet_counter == 'e')
                    alphabet_counter++;
                name = std::string(1, alphabet_counter++);
            }
            else
            {
                name = "c"+std::to_string(i);
            }

            ctx.name_vars[name] = i;
            ctx.var_names[i] = name;
            ctx.coefficients.push_back(i);
        }
        return true;
    }},
    {"let", [](context& ctx, const std::vector<parameter>& parameters)->bool
    {
        std::string name;
        double value;
        PARAMS(name, value);
        if(!ctx.name_vars.count(name))
        {
            printf("No such variable: %s\n", name.c_str());
            return false;
        }
        variable id = ctx.name_vars[name];
        ctx.p = assign(ctx.p, id, value);
        return true;
    }},
    {"differentiate", [](context& ctx, const std::vector<parameter>& parameters)->bool
    {
        std::string name;
        PARAMS(name);

        if(!ctx.name_vars.count(name))
        {
            printf("No such variable: %s\n", name.c_str());
            return false;
        }
        variable id = ctx.name_vars[name];
        std::optional<polynomial> result = differentiate(ctx.p, id);
        if(!result.has_value())
        {
            printf("Differentiation failed for polynomial %s\n", polynomial_to_string(ctx, ctx.p, true, false).c_str());
            return false;
        }
        ctx.p = *result;
        return true;
    }},
    {"pin", [](context& ctx, const std::vector<parameter>& parameters)->bool
    {
        struct input_constraint
        {
            double value;
            variable id;
        };
        struct output_constraint
        {
            double value;
            int derivative;
            variable id;
        };
        std::vector<input_constraint> input_constraints;
        std::vector<output_constraint> output_constraints;

        std::string output_axis_name = ctx.var_names[ctx.axes.back()];
        for(int i = 0; i < parameters.size(); i+=2)
        {
            std::string label;
            double value;
            int res = get_params(parameters.data()+i, std::min(parameters.size()-i, 2lu), label, value);
            if(res != 0)
            {
                printf("Pin parameters must be in the format [']<axis>=<value>.\n");
                return false;
            }

            const char* name = label.c_str();
            if(*name == '\'' || label == output_axis_name)
            {
                output_constraint oc;
                oc.value = value;
                oc.derivative = 0;
                while(*name == '\'')
                {
                    oc.derivative++;
                    name++;
                }
                if(!ctx.name_vars.count(name))
                {
                    printf("Can't differentiate over unknown variable %s!\n", name);
                    return false;
                }
                oc.id = ctx.name_vars[name];
                output_constraints.push_back(oc);
            }
            else
            {
                input_constraint ic;
                ic.value = value;
                if(!ctx.name_vars.count(name))
                {
                    printf("Can't constrain over unknown variable %s!\n", name);
                    return false;
                }
                ic.id = ctx.name_vars[name];
                input_constraints.push_back(ic);
            }
        }

        // Output constraints have to be applied one-by-one.
        std::vector<polynomial> target = {ctx.p};
        for(const output_constraint& o: output_constraints)
        {
            polynomial zero = target[0];
            for(int i = 0; i < o.derivative; ++i)
            {
                std::optional<polynomial> res = differentiate(zero, o.id);
                if(!res.has_value())
                {
                    printf("Differentiation failed!\n");
                    return false;
                }
                zero = *res;
            }
            for(const input_constraint& ic: input_constraints)
                zero = assign(zero, ic.id, ic.value);
            zero.terms.push_back(term{-o.value, {}});

            if(!pin(zero, ctx.axes.data(), ctx.axes.size()-1, target))
            {
                printf("Pinning failed!\n");
                return false;
            }
        }
        ctx.p = target[0];
        return true;
    }},
    /*
    {"reassign-names", [](polynomial& p, const std::vector<parameter>& parameters)->bool
    {
        NO_PARAMS();

        p = p.reassign_variable_names();
        return true;
    }},
    */
};

int main(int argc, char** argv)
{
    if(argc != 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "Usage: %s <path-to-command-list>\n%s", argv[0], command_list_format);
        return 1;
    }

    context ctx;

    std::string src_str;
    bool success = read_text_file(argv[1], src_str);
    if(!success)
    {
        fprintf(stderr, "Failed to open %s\n", argv[1]);
        return 2;
    }
    const char* command_list = src_str.c_str();
    skip_whitespace(command_list);
    while(*command_list)
    {
        if(*command_list == '#')
        {
            // Read until newline.
            while(*command_list != '\n' && *command_list) command_list++;
            skip_whitespace(command_list);
            continue;
        }

        // Read command name.
        std::string command_name;
        while(!strchr(" \t\n", *command_list) && *command_list)
        {
            command_name += *command_list;
            command_list++;
        }

        if(command_handlers.count(command_name) == 0)
        {
            fprintf(stderr, "No such command: %s\n", command_name.c_str());
            return 3;
        }

        skip_spaces(command_list);
        // Read parameters.
        std::vector<parameter> parameters;
        while(*command_list && *command_list != '\n')
        {
            if((*command_list >= '0' && *command_list <= '9') || *command_list == '-')
            {
                double d;
                if(!read_double(command_list, d))
                {
                    fprintf(stderr, "Failed to parse number: %s\n", read_token(command_list).c_str());
                    return 3;
                }
                else parameters.push_back(d);
            }
            else
            {
                std::string token = read_token(command_list);
                parameters.push_back(token);
            }
            skip_spaces(command_list);
        }

        auto it = command_handlers.find(command_name);
        if(!it->second(ctx, parameters))
        {
            fprintf(stderr, "Command \"%s\" failed, exiting.\n", command_name.c_str());
            return 4;
        }
        skip_whitespace(command_list);
    }
    return 0;
}
