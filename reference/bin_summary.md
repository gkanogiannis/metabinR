# Summarize the bins in a metabinR result

Counts reads in each observed bin and summarizes their distance to the
assigned bin. Distances retain the scale of the binning algorithm; they
are not probabilities or comparable across algorithms.

## Usage

``` r
bin_summary(x)
```

## Arguments

- x:

  A
  [MetabinResult](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md)
  object.

## Value

A
[S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with `bin`, `n_reads`, `proportion`, `mean_distance`, and
`median_distance`. Only observed bins are included.

## Examples

``` r
res <- composition_based_binning(
    system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
    dryRun = TRUE, kMerSizeCB = 2, numOfClustersCB = 2
)
bin_summary(res)
#> DataFrame with 2 rows and 5 columns
#>           bin   n_reads proportion mean_distance median_distance
#>   <character> <integer>  <numeric>     <numeric>       <numeric>
#> 1           2     19647   0.736836       36.0165              34
#> 2           1      7017   0.263164       47.3955              46
```
