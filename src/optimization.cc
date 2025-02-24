#include "optimization.hh"
#include <cstdio>

matrix least_squares(const matrix& A, const matrix& b)
{
    auto[Q, R] = qr_decompose(A);

    printf("A = %s\n", to_string(A, "    ").c_str());
    printf("b = %s\n", to_string(b, "    ").c_str());
    printf("Q = %s\n", to_string(Q, "    ").c_str());
    printf("R = %s\n", to_string(R, "    ").c_str());
    printf("QR = %s\n", to_string(mul(Q, R).value(), "     ").c_str());
    printf("QQ = %s\n", to_string(mul(Q, transpose(Q)).value(), "     ").c_str());

    matrix Qt = transpose(Q);

    matrix v = mul(Qt, b).value();
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
