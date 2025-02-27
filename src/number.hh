#ifndef NUMBER_HH
#define NUMBER_HH
#include <cmath>
#include <cfloat>
#include <string>

// This class is effectively just a double, but it has a biased rounding towards
// zero. This makes comparisons more reliable. It's not just a dumb flat
// epsilon, it depends on the scale of the sum that would cause the zero.
class number
{
public:
    static constexpr unsigned epsilon_ulps = 100;

    inline number(double d = 0.0): val(d) {}

    // This is the interesting one, it checks if the result is within rounding
    // error to zero.
    inline number operator+(const number& other) const
    {
        double res = val + other.val;
        double aval = fabs(val);
        double bval = fabs(other.val);
        double err = epsilon_ulps *
            ((nextafter(aval, DBL_MAX) - aval) +
            (nextafter(bval, DBL_MAX) - bval));

        return fabs(res) < err ? number(copysign(0.0, res)) : number(res);
    }

    // These do nothing interesting, they just copy the double's functionality.
    inline number operator*(const number& other) const
    { return number(val * other.val); }
    inline number operator/(const number& other) const
    { return number(val / other.val); }
    inline number operator-() const
    { return number(-val); }
    inline number operator-(const number& other) const
    { return *this + (-other); }
    inline number& operator*=(const number& other)
    { return *this = *this * other; }
    inline number& operator/=(const number& other)
    { return *this = *this / other; }
    inline number& operator+=(const number& other)
    { return *this = *this + other; }
    inline number& operator-=(const number& other)
    { return *this = *this - other; }
    inline bool operator==(const number& other) const
    { return (*this-other).val == 0.0; }
    inline bool operator!=(const number& other) const
    { return (*this-other).val != 0.0; }
    inline bool operator<(const number& other) const
    { return *this != other && val < other.val; }
    inline bool operator>(const number& other) const
    { return *this != other && val > other.val; }

    explicit operator double() const
    { return val; }

private:
    double val;
};

inline number operator*(double a, const number& b)
{ return number(a) * b; }
inline number operator/(double a, const number& b)
{ return number(a) / b; }
inline number operator+(double a, const number& b)
{ return number(a) + b; }
inline number operator-(double a, const number& b)
{ return number(a) - b; }

inline std::string to_string(number n)
{
    double val = (double)n;
    return double(int32_t(val)) == val ?
        std::to_string(int32_t(val)) :
        std::to_string(val);
}

inline number pow(number base, number exp)
{
    return pow((double)base, (double)exp);
}

inline number sign(number v) { return (double)v < 0 ? -1 : 1; }

inline number fabs(number v) { return number(fabs((double)v)); }

#endif
