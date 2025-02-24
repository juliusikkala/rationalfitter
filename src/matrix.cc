#include "matrix.hh"
#include "math.hh"
#include <algorithm>
#include <cmath>

matrix matrix::zeroes(unsigned w, unsigned h)
{
    matrix m;
    m.w = w;
    m.h = h;
    m.values.resize(w*h, 0.0);
    return m;
}

matrix matrix::eyes(unsigned w, unsigned h, double value)
{
    matrix m;
    m.w = w;
    m.h = h;
    m.values.resize(w*h, 0.0);
    for(int i = 0; i < std::min(w, h); ++i)
        m.values[i+i*w] = value;
    return m;
}

matrix matrix::vector(const std::vector<double>& values)
{
    matrix m;
    m.w = 1;
    m.h = values.size();
    m.values = values;
    return m;
}

matrix matrix::from_values(unsigned w, unsigned h, const std::vector<double>& values)
{
    matrix m;
    m.w = w;
    m.h = h;
    m.values = values;
    return m;
}

matrix matrix::column(unsigned x) const
{
    matrix c;
    c.w = 1;
    c.h = h;
    c.values.resize(h);
    for(unsigned y = 0; y < h; ++y)
        c.values[y] = values[x+y*w];
    return c;
}

matrix matrix::row(unsigned y) const
{
    matrix r;
    r.w = w;
    r.h = 1;
    r.values.resize(w);
    for(unsigned x = 0; x < w; ++x)
        r.values[x] = values[x+y*w];
    return r;
}

matrix matrix::clip(unsigned x, unsigned y, unsigned w, unsigned h) const
{
    matrix c;
    c.w = w;
    c.h = h;
    c.values.resize(w * h);

    for(unsigned i = 0; i < w; ++i)
    for(unsigned j = 0; j < h; ++j)
        c.values[i+j*w] = values[x+i+(y+j)*this->w];

    return c;
}

void matrix::insert(unsigned x, unsigned y, const matrix& other)
{
    for(unsigned i = 0; i < other.w; ++i)
    for(unsigned j = 0; j < other.h; ++j)
        values[x+i+(y+j)*this->w] = other(i, j);
}

bool operator==(const matrix& a, const matrix& b)
{
    return a.w == b.w && a.h == b.h && a.values == b.values;
}

double length(const matrix& m)
{
    double sum = 0;
    for(double v: m.values)
        sum += v*v;
    return sqrt(sum);
}

std::optional<matrix> mul(const matrix& a, const matrix& b)
{
    if(a.w != b.h)
        return {};
    matrix result = matrix::zeroes(b.w, a.h);
    for(unsigned x = 0; x < result.w; ++x)
    for(unsigned y = 0; y < result.h; ++y)
    for(unsigned t = 0; t < a.w; ++t)
        result(x, y) += a(t, y) * b(x, t);
    return result;
}

matrix mul(const matrix& a, double b)
{
    matrix result = matrix::zeroes(a.w, a.h);
    for(unsigned x = 0; x < result.w; ++x)
    for(unsigned y = 0; y < result.h; ++y)
        result(x, y) = a(x, y) * b;
    return result;
}

std::optional<matrix> add(const matrix& a, const matrix& b)
{
    if(a.w != b.w || a.h != b.h)
        return {};
    matrix result = matrix::zeroes(a.w, a.h);
    for(unsigned x = 0; x < result.w; ++x)
    for(unsigned y = 0; y < result.h; ++y)
        result(x, y) = a(x, y) + b(x, y);
    return result;
}

std::optional<matrix> sub(const matrix& a, const matrix& b)
{
    if(a.w != b.w || a.h != b.h)
        return {};
    matrix result = matrix::zeroes(a.w, a.h);
    for(unsigned x = 0; x < result.w; ++x)
    for(unsigned y = 0; y < result.h; ++y)
        result(x, y) = a(x, y) - b(x, y);
    return result;
}

matrix transpose(const matrix& m)
{
    matrix t;
    t.w = m.h;
    t.h = m.w;
    t.values.resize(m.values.size());
    for(unsigned x = 0; x < t.w; ++x)
    for(unsigned y = 0; y < t.h; ++y)
    {
        t(x, y) = m(y, x);
    }
    return t;
}

std::tuple<matrix, matrix> qr_decompose(const matrix& A)
{
    if(A.w >= A.h)
        return {{}, {}};
    matrix Qt = matrix::eyes(A.h, A.h);
    matrix R = A;
    for(int i = 0; i < std::min(A.w, A.h-1); ++i)
    {
        matrix x = R.column(i).clip(0, i, 1, R.h-i);
        double alpha = -sign(x(0, 0)) * length(x);
        x(0, 0) -= alpha;

        double u_length = length(x);
        x = mul(x, 1.0 / u_length);

        matrix xx = mul(x, transpose(x)).value();

        matrix Qi_small = sub(matrix::eyes(R.h-i, R.h-i), mul(xx, 2.0)).value();
        matrix Qi_expand = matrix::eyes(R.h, R.h, 1.0);
        Qi_expand.insert(i, i, Qi_small);

        Qt = mul(Qi_expand, Qt).value();
        R = mul(Qi_expand, R).value();
    }
    return {transpose(Qt), mul(Qt, A).value()};
}

std::string to_string(const matrix& m, const char* line_prefix)
{
    std::vector<std::string> entries(m.w*m.h);
    for(unsigned y = 0; y < m.h; ++y)
    for(unsigned x = 0; x < m.w; ++x)
    {
        entries[x+y*m.w] = std::to_string(m(x, y));
    }

    // Ensure matching column widths
    for(unsigned x = 0; x < m.w; ++x)
    {
        size_t w = 0;
        for(unsigned y = 0; y < m.h; ++y)
            w = std::max(w, entries[x+y*m.w].length());

        for(unsigned y = 0; y < m.h; ++y)
        {
            std::string& entry = entries[x+y*m.w];
            if(entry.length() < w)
                entry.insert(0, w-entry.length(), ' ');
        }
    }

    std::string result;
    for(unsigned y = 0; y < m.h; ++y)
    {
        if(y != 0) result += line_prefix;
        result += "[ ";
        for(unsigned x = 0; x < m.w; ++x)
        {
            if(x != 0)
                result += ", ";
            result += entries[x+y*m.w];
        }
        result += " ]\n";
    }
    return result;
}
