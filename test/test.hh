#ifndef TEST_HH
#define TEST_HH
#include <cstdlib>
#include <cstdio>

inline unsigned& fail_counter() { static unsigned c = 0; return c; }
#define FAIL(explanation) { fprintf(stderr, "Test \"%s\" failed\n", explanation); fail_counter()++; }
#define CASE(expr) if(!(expr)) FAIL(#expr)
#define FINISH return fail_counter() > 255 ? 255 : fail_counter();

#endif
