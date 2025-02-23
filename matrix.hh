#ifndef POLYNOMIALPINNER_MATRIX_HH
#define POLYNOMIALPINNER_MATRIX_HH
#include <vector>
#include <optional>
#include <tuple>
#include <string>

struct matrix
{
    unsigned w = 0;
    unsigned h = 0;
    std::vector<double> values;

    double& operator()(unsigned x, unsigned y) { return values[x+y*w]; }
    double operator()(unsigned x, unsigned y) const { return values[x+y*w]; }

    matrix column(unsigned x) const;
    matrix row(unsigned y) const;

    matrix clip(unsigned x, unsigned y, unsigned w, unsigned h) const;
    void insert(unsigned x, unsigned y, const matrix& other);

    static matrix zeroes(unsigned w, unsigned h);
    static matrix eyes(unsigned w, unsigned h, double value = 1.0);
    static matrix vector(const std::vector<double>& values);
    static matrix from_values(unsigned w, unsigned h, const std::vector<double>& values);
};

double length(const matrix& m);
std::optional<matrix> mul(const matrix& a, const matrix& b);
matrix mul(const matrix& a, double b);
std::optional<matrix> add(const matrix& a, const matrix& b);
std::optional<matrix> sub(const matrix& a, const matrix& b);
matrix transpose(const matrix& m);
std::tuple<matrix, matrix> qr_decompose(const matrix& m);

std::string to_string(const matrix& m, const char* line_prefix = "");

#endif
