# Null default operator

Provides a convenient infix operator for replacing `NULL` values.
`x \%||\% y` is equivalent to `if (is.null(x)) y else x`.

## Usage

``` r
x %||% y
```

## Arguments

- x, y:

  Any R objects.

## Examples

``` r
NULL %||% "fallback"
#> [1] "fallback"
list(a = 1) %||% list(b = 2)
#> $a
#> [1] 1
#> 
```
