# Evaluate read-level bin assignments against known origins

Joins assignments to a truth table by read identifier, then calculates
per-bin purity, per-origin best-bin recovery, and the adjusted Rand
index (ARI). These are read-count metrics. They do not measure genome or
metagenome-assembled genome completeness, and they are not CAMI/AMBER
base-pair-weighted scores. Bin labels are treated as arbitrary
identifiers.

## Usage

``` r
evaluate_bins(x, truth, id = "read_id", label = "genome_id")
```

## Arguments

- x:

  A
  [MetabinResult](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md)
  object.

- truth:

  A `data.frame` or
  [S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  containing read IDs and known origins or classes.

- id:

  Name of the read-ID column in `truth`. The result's `read_id` column
  is matched to this column.

- label:

  Name of the truth-label column in `truth`.

## Value

A list with `confusion` (bins as rows, truth groups as columns),
`per_bin` (read count, dominant truth group, purity), `per_origin` (read
count, dominant bin, best-bin recovery), and `overall` (read count,
numbers of bins and truth groups, weighted purity, weighted recovery,
and ARI). Purity is the dominant truth count divided by bin size;
best-bin recovery is the largest count for an origin in any one bin
divided by that origin's read count. Both weighted measures sum the
relevant dominant counts and divide by the total number of evaluated
reads. ARI is `NA` for fewer than two reads.

## Details

Every result read must have exactly one matching truth row. Extra truth
rows are ignored. Duplicate or missing identifiers are errors. If
several truth groups or bins tie for the maximum count, the first
encountered is reported as dominant.

## Examples

``` r
x <- composition_based_binning(
    system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
    dryRun = TRUE, kMerSizeCB = 2, numOfClustersCB = 2
)
truth <- read.delim(system.file(
    "extdata", "reads_mapping.tsv.gz", package = "metabinR"
))
evaluation <- evaluate_bins(x, truth, id = "anonymous_read_id",
                            label = "genome_id")
evaluation$overall
#> DataFrame with 1 row and 6 columns
#>     n_reads    n_bins n_origins weighted_purity weighted_recovery
#>   <integer> <integer> <integer>       <numeric>         <numeric>
#> 1     26664         2        10        0.598822          0.907478
#>   adjusted_rand_index
#>             <numeric>
#> 1            0.255206
```
