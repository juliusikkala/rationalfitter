#include "math.hh"
#include "test.hh"

int main()
{
    // ipow()
    CASE(ipow(42.0, 0) == 1.0);
    CASE(ipow(0.5, 2) == 0.25);
    CASE(ipow(-100.1, 1) == -100.1);
    CASE(ipow(2.0, 10) == 1024.0);
    CASE(ipow(0.0, 0) == 1);
    CASE(ipow(0.0, 1) == 0);

    // sign()
    CASE(sign(-1.0) == -1.0);
    CASE(sign(15.0) == 1.0);
    CASE(sign(0) == 1.0);

    FINISH;
}
