#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <functional>
#include <algorithm>
#include <variant>
#include <optional>

const char* command_list_format = R"(Command list format:
'#' starts a comment row and is not parsed.

'domain <degree> <dimension>' initializes the polynomial.
'pin x=<a> y=<b>' makes the polynomial equal to 'b' at 'a'.
'printl' prints the current state of the polynomial, splitting each coefficient on a new line.
'print' prints the current state of the polynomial.
)";

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

int get_axis(const char* label)
{
    int axis;
    switch(*label)
    {
    case 'x':
    case 'X':
        axis = 0;
        break;
    case 'y':
    case 'Y':
        axis = 1;
        break;
    case 'z':
    case 'Z':
        axis = 2;
        break;
    case 'w':
    case 'W':
        axis = 3;
        break;
    default:
        return -1;
    }
    label++;
    return *label == 0 ? axis : -1;
}

bool is_integer(double d)
{
    return d == double(int64_t(d));
}

struct pin_constraint
{
    static bool create(
        const char* label,
        double value,
        pin_constraint& constraint
    ){
        constraint.derivative = false;
        constraint.value = value;
        if(*label == 'd')
        {
            constraint.derivative = true;
            label++;
        }
        if(!*label) return false;
        constraint.axis = get_axis(label);
        return constraint.axis >= 0;
    }

    int axis;
    bool derivative;
    double value;
};

struct coefficient;
struct var_weight
{
    char name; // 0 is a special variable marking constant values.
    double weight;
};

struct coefficient
{
    std::vector<unsigned> degrees;
    std::vector<var_weight> sum;

    void regroup()
    {
        for(unsigned i = 0; i < sum.size(); ++i)
        for(unsigned j = i+1; j < sum.size();)
        {
            var_weight& parent = sum[i];
            var_weight& child = sum[j];
            if(parent.name == child.name)
            {
                parent.weight += child.weight;
                sum.erase(sum.begin()+j);
            }
            else ++j;
        }

        for(auto it = sum.begin(); it != sum.end();)
        {
            if(fabs(it->weight) < 1e-14)
                it = sum.erase(it);
            else ++it;
        }
    }
};

struct polynomial
{
    std::string coefficient_name(unsigned i) const
    {
        const coefficient& coef = coefficients[i];
        // Pinned!
        std::string form;
        if(coef.sum.size() > 1)
            form = "(";

        bool first = true;
        for(const var_weight& var: coef.sum)
        {
            float weight = var.weight;
            if(weight == 0 && coef.sum.size() > 1)
                continue; // '0'
            if(!first)
            {
                form += " ";
                if(weight >= 0)
                    form += "+ ";
                else
                {
                    form += "- ";
                    weight = fabs(var.weight);
                }
            }
            else if(weight < 0)
            {
                form += "-";
                weight = fabs(var.weight);
            }

            if(weight != 1 || !var.name)
            {
                form += is_integer(weight) ? 
                    std::to_string(int64_t(weight)) : 
                    std::to_string(weight);
                if(var.name)
                    form += " * ";
            }

            if(var.name)
                form += std::string(1, var.name);

            first = false;
        }

        if(coef.sum.size() > 1)
            form += ")";
        return form;
    }

    void reset(unsigned degree, unsigned dimension)
    {
        output_axis = dimension;
        unsigned count = 1;
        for(unsigned i = 0; i < dimension; ++i)
            count *= (degree+1);
        coefficients.resize(count);
        for(unsigned i = 0; i < coefficients.size(); ++i)
        {
            coefficient c;
            c.degrees.resize(dimension);
            unsigned id = i;
            for(int d = 0; d < dimension; ++d)
            {
                c.degrees[d] = id % (degree+1);
                id /= (degree+1);
            }
            c.sum = {{char('a'+i), 1}};
            coefficients[i] = c;
        }
    }

    polynomial regroup() const
    {
        polynomial copy = *this;
        for(unsigned i = 0; i < copy.coefficients.size(); ++i)
        for(unsigned j = i+1; j < copy.coefficients.size();)
        {
            coefficient& parent = copy.coefficients[i];
            coefficient& child = copy.coefficients[j];
            if(parent.degrees == child.degrees)
            {
                parent.sum.insert(parent.sum.end(), child.sum.begin(), child.sum.end());
                copy.coefficients.erase(copy.coefficients.begin()+j);
            }
            else ++j;
        }
        for(auto it = copy.coefficients.begin(); it != copy.coefficients.end();)
        {
            it->regroup();
            if(it->sum.size() == 0)
                it = copy.coefficients.erase(it);
            else ++it;
        }
        return copy;
    }

    polynomial resolve(int axis, double at) const
    {
        polynomial copy = *this;

        for(auto it = copy.coefficients.begin(); it != copy.coefficients.end();)
        {
            coefficient& c = *it;
            double mul = 1;
            for(int i = 0; i < c.degrees[axis]; ++i)
                mul *= at;

            c.degrees[axis] = 0;
            if(mul == 0)
            {
                it = copy.coefficients.erase(it);
                continue;
            }
            else if(mul != 1)
            {
                for(var_weight& w: c.sum)
                {
                    w.weight *= mul;
                }
            }
            ++it;
        }
        
        return copy.regroup();
    }

    polynomial differentiate(int axis) const
    {
        polynomial copy = *this;

        for(auto it = copy.coefficients.begin(); it != copy.coefficients.end();)
        {
            coefficient& c = *it;
            double mul = 0;
            for(int i = 0; i < c.degrees.size(); ++i)
            {
                if(i != axis)
                    c.degrees[i] = 0;
                else if(c.degrees[i] > 0)
                {
                    mul = c.degrees[i];
                    c.degrees[i]--;
                }
            }
            if(mul == 0)
            {
                it = copy.coefficients.erase(it);
                continue;
            }
            else if(mul != 1)
            {
                for(var_weight& w: c.sum)
                    w.weight *= mul;
            }
            ++it;
        }
        
        return copy.regroup();
    }

    polynomial replace(char name, const std::vector<var_weight>& sum) const
    {
        polynomial copy = *this;

        for(coefficient& c: copy.coefficients)
        {
            double weight = 0;
            for(auto it = c.sum.begin(); it != c.sum.end(); )
            {
                if(it->name == name)
                {
                    weight += it->weight;
                    it = c.sum.erase(it);
                }
                else ++it;
            }

            if(weight != 0)
            {
                for(const var_weight& w: sum)
                    c.sum.push_back({w.name, weight * w.weight});
            }
        }

        return copy.regroup();
    }

    std::optional<polynomial> pin(const std::vector<pin_constraint>& constraints) const
    {
        std::vector<pin_constraint> output_constraints;
        std::vector<pin_constraint> input_constraints;
        for(const pin_constraint& pc: constraints)
        {
            if(pc.derivative || pc.axis == output_axis)
                output_constraints.push_back(pc);
            else
                input_constraints.push_back(pc);
        }

        if(output_constraints.size() == 0)
        {
            printf("No output constraint in pin!\n");
            return {};
        }

        if(output_constraints.size() > 1)
        {
            printf("Multiple output constraints in pin!\n");
            return {};
        }
        pin_constraint& output_constraint = output_constraints[0];

        polynomial constrained = *this;
        if(output_constraint.derivative)
            constrained = constrained.differentiate(output_constraint.axis);

        for(const pin_constraint& pc: input_constraints)
            constrained = constrained.resolve(pc.axis, pc.value);

        constrained.coefficients.push_back(
            coefficient{std::vector<unsigned>(output_axis, 0), {{0, -output_constraint.value}}}
        );
        constrained = constrained.regroup();

        polynomial result = *this;

        // Go over all coefficients of the constraint and a way to make them zero.
        while(constrained.coefficients.size() != 0)
        {
            coefficient c = constrained.coefficients.back();
            constrained.coefficients.pop_back();

            char replace_name = 0;
            double replace_weight = 0;
            for(var_weight& w: c.sum)
            {
                if(w.name > replace_name)
                {
                    replace_name = w.name;
                    replace_weight = w.weight;
                }
            }

            if(replace_name == 0 || replace_weight == 0)
            {
                printf("Impossible to pin; no free variables able to achieve this pin are left!\n");
                return {};
            }

            std::vector<var_weight> sum;
            for(var_weight& w: c.sum)
            {
                if(w.name != replace_name)
                {
                    w.weight = -w.weight / replace_weight;
                    sum.push_back(w);
                }
            }

            result = result.replace(replace_name, sum);
            constrained = constrained.replace(replace_name, sum);
        }
        return result;
    }

    void print(bool separate_lines)
    {
        bool first = true;
        for(unsigned i = 0; i < coefficients.size(); ++i)
        {
            std::string name = coefficient_name(i);
            const std::vector<unsigned>& degrees = coefficients[i].degrees;
            bool all_zero = true;
            for(unsigned n: degrees)
                all_zero = all_zero && n == 0;

            if(name == "1" && !all_zero) name = "";
            if(name == "0" || name.size() == 0) continue;
            if(!first)
            {
                if(separate_lines)
                    printf("\n");
                else printf(" ");
                if(name[0] != '-')
                    printf("+ ");
                else
                {
                    printf("- ");
                    name.erase(name.begin());
                }
            }
            printf("%s", name.c_str());

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
                    term += "^"+std::to_string(degrees[d]);
            }

            if(term.size() > 0)
                printf(" * %s", term.c_str());
            first = false;
        }
        if(first)
            printf("0");
        printf("\n");
    }

    unsigned output_axis;
    std::vector<coefficient> coefficients;
};

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
using command_handler = std::function<bool(polynomial&, const std::vector<parameter>& parameters)>;

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

const std::unordered_map<std::string, command_handler> command_handlers = {
    {"printl", [](polynomial& p, const std::vector<parameter>& parameters)->bool
    {
        NO_PARAMS();
        p.print(true);
        printf("\n");
        return true;
    }},
    {"print", [](polynomial& p, const std::vector<parameter>& parameters)->bool
    {
        NO_PARAMS();
        p.print(false);
        return true;
    }},
    {"domain", [](polynomial& p, const std::vector<parameter>& parameters)->bool
    {
        int degree;
        int dimension;
        PARAMS(degree, dimension);
        if(degree < 0)
        {
            printf("Degree must be a positive integer!\n");
            return false;
        }
        if(dimension < 0 || dimension > 4)
        {
            printf("Currently, only dimensions [0, 4] are supported.\n");
            return false;
        }
        p.reset(degree, dimension);
        return true;
    }},
    {"resolve", [](polynomial& p, const std::vector<parameter>& parameters)->bool
    {
        std::string label;
        double value;
        PARAMS(label, value);
        pin_constraint pc;
        if(!pin_constraint::create(label.c_str(), value, pc))
        {
            printf("Unrecognized resolve axis %s\n", label.c_str());
            return false;
        }
        if(pc.derivative)
        {
            printf("Derivatives are not resolvable yet\n");
            return false;
        }
        p = p.resolve(pc.axis, pc.value);
        return true;
    }},
    {"differentiate", [](polynomial& p, const std::vector<parameter>& parameters)->bool
    {
        std::string label;
        PARAMS(label);

        int axis = get_axis(label.c_str());

        if(axis < 0)
        {
            printf("Unrecognized differentiation axis %s\n", label.c_str());
            return false;
        }
        p = p.differentiate(axis);
        return true;
    }},
    {"pin", [](polynomial& p, const std::vector<parameter>& parameters)->bool
    {
        std::vector<pin_constraint> constraints;
        for(int i = 0; i < parameters.size(); i+=2)
        {
            std::string label;
            double value;
            int res = get_params(parameters.data()+i, std::min(parameters.size()-i, 2lu), label, value);
            if(res != 0)
            {
                printf("Pin parameters must be in the format [d]<axis>=<value>.\n");
                return false;
            }

            pin_constraint pc;
            if(!pin_constraint::create(label.c_str(), value, pc))
            {
                printf("Unrecognized pin constraint axis %s\n", label.c_str());
                return false;
            }
            constraints.push_back(pc);
        }
        auto val = p.pin(constraints);
        if(val.has_value())
        {
            p = *val;
            return true;
        }
        return false;
    }}
};

int main(int argc, char** argv)
{
    if(argc != 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "Usage: %s <path-to-command-list>\n%s", argv[0], command_list_format);
        return 1;
    }

    polynomial pol;

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
        if(!it->second(pol, parameters))
        {
            fprintf(stderr, "Command \"%s\" failed, exiting.\n", command_name.c_str());
            return 4;
        }
        skip_whitespace(command_list);
    }
    return 0;
}
