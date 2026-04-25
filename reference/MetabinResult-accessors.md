# Accessors for [MetabinResult](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md)

Accessors for
[MetabinResult](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md)

## Usage

``` r
assignments(x)

nClusters(x)

parameters(x)

algorithm(x)

# S4 method for class 'MetabinResult'
assignments(x)

# S4 method for class 'MetabinResult'
parameters(x)

# S4 method for class 'MetabinResult'
algorithm(x)

# S4 method for class 'MetabinResult'
nClusters(x)
```

## Arguments

- x:

  A
  [MetabinResult](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md)
  object.

## Value

`assignments()` returns a
[`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with the per-read cluster assignments and distances. `nClusters()`
returns an integer scalar: the number of clusters inferred by the
algorithm. `parameters()` returns the list of arguments passed to the
binning function. `algorithm()` returns the algorithm tag (`"AB"`,
`"CB"`, or `"ABxCB"`).

## Examples

``` r
res <- abundance_based_binning(
    system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
    dryRun = TRUE, kMerSizeAB = 4, numOfClustersAB = 2
)
assignments(res)
#> DataFrame with 26664 rows and 4 columns
#>           read_id        AB      AB.1      AB.2
#>       <character> <integer> <numeric> <numeric>
#> 1          S0R0/1         1  0.569785         0
#> 2          S0R0/2         1  0.622971         0
#> 3          S0R1/1         1  0.635513         0
#> 4          S0R1/2         1  0.698513         0
#> 5          S0R2/1         1  0.728616         0
#> ...           ...       ...       ...       ...
#> 26660  S0R13329/2         1  0.733242         0
#> 26661  S0R13330/1         1  0.647232         0
#> 26662  S0R13330/2         1  0.673660         0
#> 26663  S0R13331/2         1  0.701561         0
#> 26664  S0R13331/1         1  0.715246         0
nClusters(res)
#> [1] 1
algorithm(res)
#> [1] "AB"
parameters(res)$kMerSizeAB
#> [1] 4
```
