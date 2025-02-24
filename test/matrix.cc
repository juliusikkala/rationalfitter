#include "matrix.hh"
#include "test.hh"
#include <cmath>

int main()
{
    // Creation: zeroes
    matrix zero_matrix = matrix::zeroes(5, 10);
    CASE(zero_matrix.w == 5);
    CASE(zero_matrix.h == 10);
    CASE(zero_matrix.values.size() == zero_matrix.w * zero_matrix.h);
    for(double v: zero_matrix.values)
    {
        CASE(v == 0);
    }

    // Creation: eyes
    matrix eye_matrix = matrix::eyes(3, 5, 5.5);
    CASE(eye_matrix.w == 3);
    CASE(eye_matrix.h == 5);
    CASE(eye_matrix.values.size() == eye_matrix.w * eye_matrix.h);
    for(unsigned i = 0; i < eye_matrix.w; ++i)
    for(unsigned j = 0; j < eye_matrix.h; ++j)
    {
        if(i == j)
        {
            CASE(eye_matrix(i, j) == 5.5);
        }
        else
        {
            CASE(eye_matrix(i, j) == 0);
        }
    }

    // Creation: vector
    matrix vector_matrix = matrix::vector({0,1,2,3,4});
    CASE(vector_matrix.w == 1);
    CASE(vector_matrix.h == 5);
    CASE(vector_matrix.values.size() == vector_matrix.w * vector_matrix.h);
    for(unsigned i = 0; i < vector_matrix.h; ++i)
    {
        CASE(vector_matrix(0, i) == i);
    }

    // Creation: from values
    matrix value_matrix = matrix::from_values(3,2, {
        0, 1, 2,
        3, 4, 5
    });
    CASE(value_matrix.w == 3);
    CASE(value_matrix.h == 2);
    CASE(value_matrix.values.size() == value_matrix.w * value_matrix.h);
    for(unsigned i = 0; i < value_matrix.w; ++i)
    for(unsigned j = 0; j < value_matrix.h; ++j)
    {
        CASE(value_matrix(i, j) == i+j*value_matrix.w);
    }

    // Indexing has been incidentally also tested above, so that's not repeated.
    matrix column_matrix = value_matrix.column(1);
    CASE(column_matrix.w == 1);
    CASE(column_matrix.h == 2);
    CASE(column_matrix.values.size() == column_matrix.w * column_matrix.h);
    CASE(column_matrix.values[0] == 1.0 && column_matrix.values[1] == 4.0);

    matrix row_matrix = value_matrix.row(1);
    CASE(row_matrix.w == 3);
    CASE(row_matrix.h == 1);
    CASE(row_matrix.values.size() == row_matrix.w * row_matrix.h);
    CASE(row_matrix.values[0] == 3.0 && row_matrix.values[1] == 4.0 && row_matrix.values[2] == 5.0);

    matrix clip_matrix = value_matrix.clip(1,1,2,1);
    CASE(clip_matrix.w == 2);
    CASE(clip_matrix.h == 1);
    CASE(clip_matrix.values.size() == clip_matrix.w * clip_matrix.h);
    CASE(clip_matrix.values[0] == 4.0 && clip_matrix.values[1] == 5.0);

    value_matrix.insert(0, 0, clip_matrix);
    CASE(value_matrix.values[0] == 4.0);
    CASE(value_matrix.values[1] == 5.0);

    // Vector ops
    matrix test_vector = matrix::vector({2,3,-1});
    CASE(length(test_vector) == sqrt(2*2+3*3+1*1));
    test_vector = transpose(test_vector);
    CASE(length(test_vector) == sqrt(2*2+3*3+1*1));

    // Matrix ops
    matrix a = matrix::from_values(4,3, {
        0, 1, 2, 3,
        4, 5, 6, 7,
        8, 9, 10, 11
    });
    matrix b = matrix::from_values(3,4, {
        -5, 1, 35,
        3, 2, 15,
        7, 3, 0,
        105, -52, 14
    });

    // Multiplication!
    CASE(!mul(a,a).has_value());

    std::optional<matrix> c = mul(a, b);
    CASE(c.has_value());
    CASE(
        *c == matrix::from_values(3, 3, {
            332, -148, 57,
            772, -332, 313,
            1212, -516, 569
        })
    );
    c = mul(b, a);
    CASE(c.has_value());
    CASE(
        *c == matrix::from_values(4, 4, {
            284, 315, 346, 377,
            128, 148, 168, 188,
            12, 22, 32, 42,
            -96, -29, 38, 105
        })
    );
    c = mul(a, 5.0);
    CASE(
        *c == matrix::from_values(4, 3, {
            0, 5, 10, 15,
            20, 25, 30, 35,
            40, 45, 50, 55
        })
    );

    // Transpose!
    matrix bt = transpose(b);
    CASE(
        bt == matrix::from_values(4, 3, {
            -5, 3, 7, 105,
            1, 2, 3, -52,
            35, 15, 0, 14
        })
    );

    // add / sub!
    CASE(!add(a,b).has_value());
    CASE(!sub(a,b).has_value());
    c = add(a, bt);
    CASE(c.has_value());
    CASE(
        *c == matrix::from_values(4, 3, {
            -5, 4, 9, 108,
            5, 7, 9, -45,
            43, 24, 10, 25
        })
    );

    c = sub(a, bt);
    CASE(c.has_value());
    CASE(
        *c == matrix::from_values(4, 3, {
            0+5, 1-3, 2-7, 3-105,
            4-1, 5-2, 6-3, 7+52,
            8-35, 9-15, 10-0, 11-14
        })
    );

    // QR decomposition
    auto[fQ, fR] = qr_decompose(bt);
    CASE(fQ.w == 0 && fQ.h == 0 && fR.w == 0 && fR.h == 0);

    auto[Q, R] = qr_decompose(b);

    CASE(Q.w == b.h);
    CASE(Q.h == b.h);
    CASE(R.w == b.w);
    CASE(R.h == b.h);

    // Q * Q^T should be an eyes matrix
    c = sub(mul(Q, transpose(Q)).value(), matrix::eyes(4, 4));
    CASE(c.has_value());
    CASE(length(*c) < 1e-15);

    // R must be an upper triangle matrix
    for(unsigned i = 0; i < R.w; ++i)
    for(unsigned j = 0; j < R.h; ++j)
    {
        if(j > i) CASE(eye_matrix(i, j) == 0);
    }

    // Q * R must be the original matrix
    c = sub(mul(Q, R).value(), b);
    CASE(c.has_value());
    CASE(length(*c) < 1e-13);

    std::string str = to_string(b);
    printf("%s\n", str.c_str());
    CASE(str == 
        "[  -5.000000,   1.000000, 35.000000 ]\n"
        "[   3.000000,   2.000000, 15.000000 ]\n"
        "[   7.000000,   3.000000,  0.000000 ]\n"
        "[ 105.000000, -52.000000, 14.000000 ]\n"
    );

    FINISH;
}
