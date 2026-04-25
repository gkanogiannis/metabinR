# MetabinResult: binning result container

An S4 class returned by \[abundance_based_binning()\],
\[composition_based_binning()\] and \[hierarchical_binning()\].

## Usage

``` r
# S4 method for class 'MetabinResult'
show(object)
```

## Arguments

- object:

  A MetabinResult.

## Value

Objects of this class are returned by the binning functions. Use
[`assignments`](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-accessors.md),
[`nClusters`](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-accessors.md),
[`parameters`](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-accessors.md),
[`algorithm`](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-accessors.md),
or [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to
access the results.

## Slots

- `assignments`:

  A
  [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  of cluster assignments. The first column is `read_id`; subsequent
  columns are algorithm-specific (see the corresponding binning
  function).

- `parameters`:

  Named list of the parameters passed to the algorithm.

- `inputs`:

  Character vector of input file paths that were processed.

- `algorithm`:

  Character scalar: one of `"AB"`, `"CB"`, `"ABxCB"`.

## Examples

``` r
res <- abundance_based_binning(
    system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
    dryRun = TRUE, kMerSizeAB = 4, numOfClustersAB = 2
)
res
#> MetabinResult (AB)
#>   reads:     26,664
#>   clusters:  1
#>   inputs:    1 file(s)
#>     - /home/runner/work/_temp/Library/metabinR/extdata/reads.metagenome.fasta.gz
is(res, "MetabinResult")
#> [1] TRUE
```
