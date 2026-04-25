# Composition based binning on metagenomic samples

This function performs composition based binning on metagenomic samples,
directly from FASTA or FASTQ files, by short kmer analysis (k\<8). See
[doi:10.1186/s12859-016-1186-3](https://doi.org/10.1186/s12859-016-1186-3)
for more details.

## Usage

``` r
composition_based_binning(
  ...,
  kMerSizeCB = 4,
  numOfClustersCB = 5,
  outputCB = "CB.cluster",
  keepQuality = FALSE,
  dryRun = FALSE,
  gzip = FALSE,
  numOfThreads = BiocParallel::bpworkers()
)
```

## Arguments

- ...:

  Input sequences. Either character paths to FASTA/FASTQ files
  (uncompressed or gzip compressed), a
  [`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  /
  [`QualityScaledDNAStringSet`](https://rdrr.io/pkg/Biostrings/man/QualityScaledXStringSet-class.html),
  or a
  [`ShortReadQ`](https://rdrr.io/pkg/ShortRead/man/ShortReadQ-class.html)
  object. Non-file inputs are staged to a temporary FASTA/FASTQ file for
  the Java backend.

- kMerSizeCB:

  kmer length for Composition based Binning.

- numOfClustersCB:

  Number of Clusters for Composition based Binning.

- outputCB:

  Output Composition based Binning Clusters files location and prefix.

- keepQuality:

  Keep fastq qualities on the output files. (will produce .fastq)

- dryRun:

  Don't write any output files.

- gzip:

  Gzip output files.

- numOfThreads:

  Number of threads to use. Defaults to
  [`bpworkers()`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html).

## Value

A
[MetabinResult](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md)
object. Its `assignments` slot is a
[`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with `numOfClustersCB + 2` columns:

- `read_id` : read identifier from fasta header

- `CB` : read was assigned to this CB cluster index

- `CB.n` : read to cluster CB.n distance

For backwards-compatible `data.frame` output use
`as.data.frame(result)`.

## References

<https://github.com/gkanogiannis/metabinR>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
res <- composition_based_binning(
    system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
    dryRun = TRUE, kMerSizeCB = 2
)
res
#> MetabinResult (CB)
#>   reads:     26,664
#>   clusters:  5
#>   inputs:    1 file(s)
#>     - /home/runner/work/_temp/Library/metabinR/extdata/reads.metagenome.fasta.gz
```
