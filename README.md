RationalFitter
==============

This tool creates polynomial & rational functions (multivariate, Nth degree)
with a given set of constraints, and allows fitting them to data.

Some functions, especially multivariate ones, can be difficult to fit
effectively in some circumstances:
* You have constraints (e.g. "Function's second derivative must be exactly 0 at x=0 y=1")
* You want a [rational function fit](https://en.wikipedia.org/wiki/Rational_function) to better deal with difficult curves
    - E.g. polynomial fits to sigmoids often require lots of coefficients, but rational functions fare much better
* You want to eliminate the least important coefficients for performance (when there are dozens or hundreds)

This project aims to make it much easier and less manual to find effective
polynomial & rational approximations in those circumstances.

While you could probably do all of this with Mathematica, this project costs
100% less to use ;)

## Features

[x] Easily generate polynomials & rational functions
[x] Place constraints on values & derivatives
[x] Output final C code for the fit
[x] Load datasets from raw binary formats & CSV
[x] Fit polynomials with linear least-squares
[x] Fit rationals with non-linear least-squares (Gauss-Newton)
[x] Automated elimination of low-importance coefficients
[ ] Tests (WIP)

(Because tests are WIP, please check the output in e.g. Desmos to make sure it
does what it should, before trusting this program. There's no warranty in any
case, see [license.](LICENSE))

## Building

The only dependency is the C++17 standard library.

```sh
cmake -S . -B build
cmake --build build
```

## Running

```sh
build/rationalfitter mymath.fit
```

where the contents `mymath.fit` follows [the next chapter](#fit-files).

## Fit files

Fit files are simple text files specifying the polynomial and the constraints placed on it.

The "language" (if you want to call it that) is very simple:
* One line per command
* All commands implicitly operate on the function you're working on
* `#` starts a single-line comment
* All numbers are doubles (though they may be interpreted as integers when necessary)
* All text parameters are strings, quotes are optional but useful if you need spaces

Commands format is `command-name unnamed-parameter named-parameter=value`, i.e.
parameters are delimited by spaces and some commands expect named parameters.

Commands are listed below.

### `polynomial`

Start with a polynomial function. The first parameter is the degree of the
polynomial, and the rest specify the axes, where the last axis is the result
axis.

Examples:
* `polynomial 2 x y` creates a 2nd degree polynomial `a + b * x + c * x^2 = y`
* `polynomial 3 x y z` creates a 3rd degree polynomial over `x` and `y`, resulting in `z`.

### `rational`

Start with a rational function. The first parameter is the degree of the
numerator polynomial, and the second parameter is the degree of the denominator
polynomial. The rest specify the axes, where the last axis is the result axis.

Examples:
* `rational 2 1 x y` creates `(a + b * x + c * x^2)/(1 + d * x) = y`

### `print`

Prints the current state of the function. Optionally, takes one of the following
format parameters: `lc`, which factors coefficients first; `c`, which prints C
code; and `multiline`, which splits the output into lines by term for readability.

Examples:
* `print` prints the current state of the function, e.g. `a + b * x + c * x^2`
* `print c`:
```c
float func(float a, float b, float c, float x)
{
    float x2 = x * x;
    float num = 0;
    num += a;
    num += x * b;
    num += x2 * c;
    return num;
}
```

### `let`

`let <variable>=<value>` sets one of the variables from the polynomial to the
given value, erasing the variable in the process.

Examples:
* `let b=1` makes `a + b * x` become `a + x`

### `pin`

Sets a constraint on the function. The input axes can be given as named
parameters to choose a point to constrain, and the output axis can be used
to specify the value that must be met. Additionally, derivatives can be
used as constraints as well, where `'x` means `dy/dx` (assuming `y` is the
output axis). Higher-order derivatives can be specified with multiple ticks,
e.g. `'''x`.

For a multivariate function, you don't need to specify all input axes. E.g.
`pin x=0 z=1` is valid, and applies for all values of `y`.

Examples:
* `pin x=0 y=-1` makes `a + b * x + c * x^2` become `- 1 + b * x + c * x^2`
* `pin x=1 'x=0` makes `a + b * x + c * x^2` become `a + b * x - 0.500000 * b * x^2`

### `differentiate`

`differentiate <axis>` differentiates the polynomial over the given axis.

Examples:
* `differentiate x` makes `a + b * x + c * x^2` become `b + 2 * c * x`

### `reassign-names`

Reassigns variable names in alphabetic order. Without this, the original
variable names from the `polynomial` or `rational` command are always preserved.

Examples:
* `reassign-names` makes `b + 2 * c * x` become `a + 2 * b * x`

### `load-float`, `load-double`, `load-int32`, `load-int64`

These commands load unformatted binary datasets into memory. They take 4 named parameters:
* `dataset=<name>` gives a name to the loaded dataset
* `path=<path>` tells where to load the file
* (optional) `offset=<offset>` gives that starting offset in bytes. Default = `0`.
* (optional) `stride=<stride>` gives the distance between entries in bytes. Default = `sizeof(type)`.

Examples:
* `load-float dataset=y_data path=zdata.bin offset=4 stride=8` loads a dataset
  called `y_data` from `zdata.bin`, and reads 32-bit floating point numbers in
  native endian starting from byte 4, and the entries are placed 8 bytes apart.

### `load-csv`

Loads a dataset in CSV format into memory. Parameters:
* `dataset=<name>` gives a name to the dataset
* `path=<path>` tells the CSV file location
* (optional) `delimiter=<delim-char>` sets the delimiter character. Default = `","`.
* (optional) `column=<column-name or column-index>` specifies the column. If the
  value is a string, the CSV file is assumed to have a header. If it's a number,
  it's just the column index directly and there must be no header. Default = `0`.

Examples:
* `load-csv dataset=zvalues path="dataset.csv" column="Height"`

### `dump-dataset`

`dump-dataset <name>` prints out the named dataset. You can use this to check
that it got loaded properly, which may be especially useful with `load-float` & co.,
where you need to get the offset and stride right to not get garbage dataj

### `fit`

`fit` searches for optimal values for the coefficients, minimizing L2 loss.
It replaces all coefficients with the found values. It takes dataset names for
each axis as parameters, e.g. `x=x_data y=y_data z=z_data`, where `z_data` and
friends refer to datasets loaded with e.g. `load-csv`.

Internally, `fit` uses linear least-squares. For polynomials, that's all. For
rational functions, it further refines the fit with the Gauss-Newton algorithm.

Named parameters, all are optional:
* `eliminate=<count>` eliminates `count` least-contributing coefficients. Default = 0.
* `maxloss=<limit>` eliminates coefficients until further removals would cause
  loss to get worse than `limit`. Default = 0 (no eliminations).
* `step=<step-size>` adjusts the step size for Gauss-Newton. Default = 0.1.
* `convergence=<limit>` adjusts the convergence heuristic for Gauss-Newton. Lower
  values let the algorithm run longer even if there is low to no improvement. Default = 0.01.
* `maxiterations=<count>` adjusts the maximum number of Gauss-Newton iterations
  before giving up. Default = 100.

Examples:
* `fit x=x_data y=y_data z=t1_data maxloss=0.001` fits a 2D function to the
  given data, and sets least important coefficients to 0 until loss would go
  worse than `0.001`.

## Examples

### 3rd degree polynomial

```
# Specify a 3rd degree 1D polynomial (y = a + b*x + c*x^2)
polynomial 3 x y

# At x = 0, ensure dy/dx = 0.
pin x=0 'x=0

# At x = 0, ensure y = 0.
pin x=0 y=0

# At x = 1, ensure dy/dx = 0.
pin x=1 'x=0

# At x = 1, ensure y = 1.
pin x=1 y=1

# Print the resulting polynomial.
print
```

The resulting polynomial is then `3 * x^2 - 2 * x^3`.
There are no free variables left, so the polynomial cannot be adjusted further.

If you set `polynomial 4 x y`, the output is instead:

```
c * x^2 + (-2 * c + 4) * x^3 + (c - 3) * x^4
```

This means that no matter the value of `c`, the polynomial will always pass through (0, 0) and (1,1), and have a slope of 0 at x=0 and x=1.

Further, by changing `print` to `print lc`:

```
(4 * x^3 - 3 * x^4) + c * (x^2 - 2 * x^3 + x^4)
```

This form avoids repeating `c` in the output.
