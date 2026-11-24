# Abundance based binning on metagenomic samples

This function performs abundance based binning on metagenomic samples,
directly from FASTA or FASTQ files, by long kmer analysis (k\>8). See
[doi:10.1186/s12859-016-1186-3](https://doi.org/10.1186/s12859-016-1186-3)
for more details.

## Usage

``` r
abundance_based_binning(
  ...,
  eMin = 1,
  eMax = 0,
  kMerSizeAB = 10,
  numOfClustersAB = 3,
  outputAB = "AB.cluster",
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

- eMin:

  Exclude kmers of less or equal count.

- eMax:

  Exclude kmers of more or equal count.

- kMerSizeAB:

  kmer length for Abundance based Binning.

- numOfClustersAB:

  Number of Clusters for Abundance based Binning.

- outputAB:

  Output Abundance based Binning Clusters files location and prefix.

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
with `numOfClustersAB + 2` columns:

- `read_id` : read identifier from fasta header

- `AB` : read was assigned to this AB cluster index

- `AB.n` : read to cluster AB.n distance

For backwards-compatible `data.frame` output use
`as.data.frame(result)`.

## References

<https://github.com/gkanogiannis/metabinR>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
res <- abundance_based_binning(
    system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
    dryRun = TRUE, kMerSizeAB = 8
)
res
#> MetabinResult (AB)
#>   reads:     26,664
#>   clusters:  3
#>   inputs:    1 file(s)
#>     - /home/runner/work/_temp/Library/metabinR/extdata/reads.metagenome.fasta.gz
head(as.data.frame(res))
#>   read_id AB       AB.1       AB.2        AB.3
#> 1  S0R0/1  1 0.09945798 0.04015390 0.014257520
#> 2  S0R0/2  1 0.09309136 0.04982843 0.009528324
#> 3  S0R1/1  1 0.07955002 0.03811272 0.025665234
#> 4  S0R1/2  1 0.06132334 0.04118644 0.029910034
#> 5  S0R2/2  2 0.01154256 0.06789413 0.029774552
#> 6  S0R2/1  1 0.04883219 0.03985131 0.037270951
```
