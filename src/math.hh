#ifndef RATIONALFITTER_MATH_HH
#define RATIONALFITTER_MATH_HH

inline double ipow(double base, unsigned exp)
{
    double result = 1;
    while(exp != 0)
    {
        if((exp&1) != 0) result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

inline double sign(double v)
{
    return v < 0 ? -1 : 1;
}

#endif
