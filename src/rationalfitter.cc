#include "polynomial.hh"
#include "rational.hh"
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
    rational r;
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
        if(poly.terms.size() > 1 && (groups.size() > 1 || indeterminate_string != "1"))
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
    return polynomial_string.size() == 0 ? "0" : polynomial_string;
}

std::string rational_to_string(
    context& ctx,
    const rational& r,
    bool group_by_axes,
    bool multiline,
    const char* pow_symbol="^"
){
    std::string num = polynomial_to_string(ctx, r.numerator, group_by_axes, multiline, pow_symbol);
    std::string denom = polynomial_to_string(ctx, r.denominator, group_by_axes, multiline, pow_symbol);

    if(denom == "1") return num;

    if(r.numerator.terms.size() >= 2) num = "("+num+")";
    if(r.denominator.terms.size() >= 2) denom = "("+denom+")";

    return num + "/" + denom;
}

std::string rational_to_c(
    context& ctx, const rational& r
){
    std::map<std::string, int> vars_max_degrees;
    bool has_roots = false;
    auto check_vars = [&](const term& t){
        for(const var_power& vp: t.mul)
        {
            if(vp.roots.has_value())
            {
                has_roots = true;
                continue;
            }
            int& degree = vars_max_degrees[ctx.var_names[vp.id]];
            degree = degree > vp.degree ? degree : vp.degree;
        }
    };
    for(const term& t: r.numerator.terms) check_vars(t);
    for(const term& t: r.denominator.terms) check_vars(t);
    if(has_roots)
        return "(cannot represent roots() in C)";
    std::string code = "float func(";

    bool first_param = true;
    for(auto& pair: vars_max_degrees)
    {
        if(pair.second > 0)
        {
            if(!first_param) code += ", ";
            code += "float " + pair.first;
            first_param = false;
        }
    }
    code += ")\n{\n";

    for(auto& pair: vars_max_degrees)
    {
        int degree = pair.second;
        std::string prevname = pair.first;
        for(unsigned j = 2; j <= degree; ++j)
        {
            std::string dname = pair.first + std::to_string(j);
            code += "    float "+dname+" = "+prevname+" * "+pair.first+";\n";
            prevname = dname;
        }
    }
    code += "    float num = 0;\n";
    for(const term& t: r.numerator.terms)
    {
        std::string tstr = term_to_string(ctx, t, "");
        if(tstr[0] == '-')
        {
            code += "    num -=";
            tstr.erase(tstr.begin());
        }
        else code += "    num += ";
        code += tstr + ";\n";
    }

    if(try_get_constant_value(r.denominator) != 1.0)
    {
        code += "    float denom = 0;\n";
        for(const term& t: r.denominator.terms)
        {
            std::string tstr = term_to_string(ctx, t, "");
            if(tstr[0] == '-')
            {
                code += "    denom -=";
                tstr.erase(tstr.begin());
            }
            else code += "    denom += ";
            code += tstr + ";\n";
        }
        code += "    return num*(1.0f/denom);\n}";
    }
    else
    {
        code += "    return num;\n}";
    }
    return code;
}

std::string rational_to_numpy(
    context& ctx, const rational& r
){
    if(try_get_constant_value(r.denominator) != 1.0)
        return "(numpy fitting code output for rationals is not supported yet)";
    const polynomial& p = r.numerator;
    std::map<indeterminate_group, polynomial> groups = 
        group_by_indeterminates(p, ctx.coefficients.data(), ctx.coefficients.size());
    std::string code;

    code += "import numpy as np\n\n";
    code += "def fit(";

    for(unsigned i = 0; i < ctx.axes.size(); ++i)
    {
        if(i != 0)
            code += ", ";
        code += ctx.var_names[ctx.axes[i]]+"_data";
    }
    code += "):\n";

    for(unsigned i = 0; i < ctx.axes.size(); ++i)
    {
        const std::string& name = ctx.var_names[ctx.axes[i]];
        code += "    "+name+" = "+name+"_data.flatten()\n";
    }

    code += "    offset = ";
    indeterminate_group constant_group;
    if(groups.count(constant_group))
    {
        polynomial& group_polynomial = groups[constant_group];
        std::set<variable> live = live_variables(group_polynomial);
        std::string poly_str = polynomial_to_string(ctx, group_polynomial, false, false, "**");

        if(live.size() != 0)
            code += poly_str+"\n";
        else
            code += "np.ones("+ctx.var_names[ctx.axes[0]]+".shape) * ("+poly_str+")\n";
    }
    else
    {
        code += "np.zeros("+ctx.var_names[ctx.axes[0]]+".shape)\n";
    }
    code += "    b = "+ctx.var_names[ctx.axes.back()]+" - offset\n";
    code += "    A = np.array([\n";

    for(auto it = groups.begin(); it != groups.end();)
    {
        // Skip constants here.
        if(it->first.indeterminates.size() == 0)
        {
            ++it;
            continue;
        }
        std::string poly_str = polynomial_to_string(ctx, it->second, false, false, "**");
        term indeterminate_term = term{1, it->first.indeterminates};
        code += "        " + poly_str;
        ++it;
        bool last = it == groups.end();
        if(!last) code += ",";
        code += " # " + term_to_string(ctx, indeterminate_term, "**");
        code += "\n";
    }

    code += "    ).T\n";

    code += "    coeff, r, rank, s = np.linalg.lstsq(A, b, rcond=None)\n";
    code += "    out = (np.dot(A, coeff) + offset).reshape("+ctx.var_names[ctx.axes[0]]+"_data.shape)\n";
    code += "    return (coeff, out)\n";
    return code;
}

bool assign_axis_names(context& ctx, const std::vector<std::string>& axis_names)
{
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
    return true;
}

void assign_variable_names(context& ctx)
{
    std::set<variable> live = live_variables(ctx.r);

    // Remove non-axis coefficient names.
    for(auto it = ctx.var_names.begin(); it != ctx.var_names.end();)
    {
        bool is_axis = false;
        for(variable axis: ctx.axes)
        {
            if(it->first == axis)
                is_axis = true;
        }

        if(!is_axis)
        {
            ctx.name_vars.erase(it->second);
            it = ctx.var_names.erase(it);
        }
        else ++it;
    }

    bool alphabet_names = live.size() < 'z'-'a'-ctx.axes.size();
    char alphabet_counter = 'a';
    unsigned i = 0;
    for(variable var: live)
    {
        if(ctx.var_names.count(var))
            continue;

        std::string name;
        if(alphabet_names)
        {
            while(ctx.name_vars.count(std::string(1, alphabet_counter)) || alphabet_counter == 'e')
                alphabet_counter++;
            name = std::string(1, alphabet_counter++);
        }
        else
        {
            name = "c"+std::to_string(i++);
        }

        ctx.name_vars[name] = var;
        ctx.var_names[var] = name;
    }
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

        std::string str;
        if(c_code)
            str = rational_to_c(ctx, ctx.r);
        else if(numpy)
        {
            str = rational_to_numpy(ctx, ctx.r);
        }
        else str = rational_to_string(
            ctx,
            ctx.r,
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

        int dimension = int(axis_names.size())-1;
        if(dimension <= 0)
        {
            printf("Dimensions must be greater than zero!\n");
            return false;
        }
        bool success = assign_axis_names(ctx, axis_names);
        if(!success)
            return false;

        variable var_counter = 0;
        ctx.r.numerator = polynomial::create(ctx.axes.data(), dimension, degree, var_counter);
        ctx.r.denominator = polynomial::create(1);
        for(variable i = 0; i < var_counter; ++i)
            ctx.coefficients.push_back(i);
        assign_variable_names(ctx);
        return true;
    }},
    {"rational", [](context& ctx, const std::vector<parameter>& parameters)->bool
    {
        int num_degree;
        int denom_degree;
        std::vector<std::string> axis_names;
        PARAMS(num_degree, denom_degree, axis_names);
        if(num_degree < 0 || denom_degree < 0)
        {
            printf("Degree must be a positive integer!\n");
            return false;
        }

        int dimension = int(axis_names.size())-1;
        if(dimension <= 0)
        {
            printf("Dimensions must be greater than zero!\n");
            return false;
        }
        bool success = assign_axis_names(ctx, axis_names);
        if(!success)
            return false;

        variable var_counter = 0;
        ctx.r.numerator = polynomial::create(ctx.axes.data(), dimension, num_degree, var_counter);
        ctx.r.denominator = polynomial::create(ctx.axes.data(), dimension, denom_degree, var_counter, true);
        for(variable i = 0; i < var_counter; ++i)
            ctx.coefficients.push_back(i);
        assign_variable_names(ctx);
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
        ctx.r = assign(ctx.r, id, value);
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
        std::optional<rational> result = differentiate(ctx.r, id);
        if(!result.has_value())
        {
            printf("Differentiation failed for %s\n", rational_to_string(ctx, ctx.r, true, false).c_str());
            return false;
        }
        ctx.r = *result;
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
        rational target = ctx.r;
        for(const output_constraint& o: output_constraints)
        {
            rational zero = target;
            for(int i = 0; i < o.derivative; ++i)
            {
                std::optional<rational> res = differentiate(zero, o.id);
                if(!res.has_value())
                {
                    printf("Differentiation failed!\n");
                    return false;
                }
                zero = *res;
            }

            for(const input_constraint& ic: input_constraints)
                zero = assign(zero, ic.id, ic.value);
            polynomial zero_poly = get_zero_polynomial(zero, o.value);

            std::vector<polynomial> vec = {target.numerator, target.denominator};
            if(!pin(zero_poly, zero.denominator, ctx.axes.data(), ctx.axes.size()-1, vec))
            {
                printf("Pinning failed!\n");
                return false;
            }
            target.numerator = vec[0];
            target.denominator = vec[1];
            target = simplify(target);
        }
        ctx.r = target;
        return true;
    }},
    {"reassign-names", [](context& ctx, const std::vector<parameter>& parameters)->bool
    {
        NO_PARAMS();
        assign_variable_names(ctx);
        return true;
    }},
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
