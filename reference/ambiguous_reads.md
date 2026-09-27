# Find reads with close competing bins

Returns reads for which the second-smallest finite distance is within
`margin` of the smallest distance. In hierarchical results, distances to
bins outside the read's parent abundance bin are `NA` and are ignored.
Reads with fewer than two finite distances are omitted. The margin is an
absolute difference on the algorithm's distance scale, not a
probability.

## Usage

``` r
ambiguous_reads(x, margin = 0.05)
```

## Arguments

- x:

  A
  [MetabinResult](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md)
  object.

- margin:

  Nonnegative maximum difference between the two smallest distances.
  Defaults to `0.05`.

## Value

A
[S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with `read_id`, `bin`, `best_distance`, `second_distance`, and
`distance_margin`, in input order.

## Examples

``` r
res <- composition_based_binning(
    system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
    dryRun = TRUE, kMerSizeCB = 2, numOfClustersCB = 2
)
ambiguous_reads(res, margin = 0.05)
#> DataFrame with 142 rows and 5 columns
#>         read_id         bin best_distance second_distance distance_margin
#>     <character> <character>     <numeric>       <numeric>       <numeric>
#> 1      S0R235/1           1            60              60               0
#> 2      S0R350/1           1            78              78               0
#> 3      S0R451/2           1            60              60               0
#> 4      S0R623/1           1            68              68               0
#> 5      S0R747/1           1            76              76               0
#> ...         ...         ...           ...             ...             ...
#> 138  S0R12993/2           1            60              60               0
#> 139  S0R12996/1           1            60              60               0
#> 140  S0R13037/2           1            64              64               0
#> 141  S0R13043/1           1            68              68               0
#> 142  S0R13268/2           1            62              62               0
```
