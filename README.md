PolynomialPinner
================

This tool creates N-D polynomials with a set of given limitations. This is
useful e.g. to ensure a least-squares fit of a function must satisfy some given
constraints.

## Building

The only dependency is the C++17 standard library.

```sh
cmake -S . -B build
cmake --build build
```

## Running

```sh
build/polynomialpinner mymath.pin
```

where the contents `mymath.pin` follows [the next chapter](#pin-files).

## Pin files

Pin files are simple text files specifying the polynomial and the constraints placed on it. Commands:

* `polynomial <degree> <dimensions>`: start with a polynomial of `<degree>` and number of `<dimensions>`s. E.g. `polynomial 3 2` creates a 3rd degree polynomial over `x` and `y`.
* `pin <in-axis>=<value> <out-axis>=<value>`: ensure the polynomial has a specific value at a given point. For a 2D polynomial, `<in-axis>` is either `x` or `y`, and `<out-axis>` can be `z`, `dx` or `dy`. You may specify multiple in-axes, `pin x=0 y=0 z=2` is valid for a 2D polynomial.
* `print`: prints the current state of the polynomial. `print multiline` splits each term on a new line, and `print lc` factors variables first. `print numpy` prints fitting code for use with Python & Numpy. `print c` prints C code that computes the fit.
* `reassign-names`: renames existing variables with successive letters.
* `differentiate <axis>`: differentiates the polynomial over given axis.
* `let <axis>=<value>`: assigns a given value to an `<axis>` (e.g. `x=1`).

There's also comment support, '#' starts a single-line comment.
See examples below on how to use this tool in practice.

### 3rd degree polynomial

```
# Specify a 3rd degree 1D polynomial (y = a + b*x + c*x^2)
polynomial 3 1

# At x = 0, ensure dy/dx = 0.
pin x=0 dx=0

# At x = 0, ensure y = 0.
pin x=0 y=0

# At x = 1, ensure dy/dx = 0.
pin x=1 dx=0

# At x = 1, ensure y = 1.
pin x=1 y=1

# Print the resulting polynomial.
print
```

The resulting polynomial is then `3 * x^2 - 2 * x^3`.
There are no free variables left, so the polynomial cannot be adjusted further.

If you set `polynomial 4 1`, the output is instead:

```
c * x^2 + (-2 * c + 4) * x^3 + (c - 3) * x^4
```

This means that no matter the value of `c`, the polynomial will always pass through (0, 0) and (1,1), and have a slope of 0 at x=0 and x=1.

Further, by changing `print` to `print lc`:

```
(4 * x^3 - 3 * x^4) + c * (x^2 - 2 * x^3 + x^4)
```

This form avoids repeating `c` in the output.
