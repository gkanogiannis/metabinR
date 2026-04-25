# Convert a [MetabinResult](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md) to a `data.frame`

Preserves backwards compatibility with the `data.frame` return of
metabinR \<= 1.x.

## Usage

``` r
# S4 method for class 'MetabinResult'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A
  [MetabinResult](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md).

- row.names, optional, ...:

  Unused; retained for S3 signature compatibility.

## Value

A base `data.frame` of the assignments.
