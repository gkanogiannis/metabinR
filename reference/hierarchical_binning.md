# Hierarchical (ABxCB) binning on metagenomic samples

This function performs hierarchical binning on metagenomic samples,
directly from FASTA or FASTQ files. First it analyzes sequences by long
kmer analysis (k\>8), as in
[`abundance_based_binning`](https://gkanogiannis.github.io/metabinR/reference/abundance_based_binning.md).
Then for each AB bin, it guesses the number of composition bins in it
and performs composition based binning by short kmer analysis (k\<8), as
in
[`composition_based_binning`](https://gkanogiannis.github.io/metabinR/reference/composition_based_binning.md).
See
[doi:10.1186/s12859-016-1186-3](https://doi.org/10.1186/s12859-016-1186-3)
for more details.

## Usage

``` r
hierarchical_binning(
  ...,
  eMin = 1,
  eMax = 0,
  kMerSizeAB = 10,
  kMerSizeCB = 4,
  genomeSize = 3e+06,
  numOfClustersAB = 3,
  outputC = "ABxCB.cluster",
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

- kMerSizeCB:

  kmer length for Composition based Binning.

- genomeSize:

  Average genome size of taxa in the metagenome data.

- numOfClustersAB:

  Number of Clusters for Abundance based Binning.

- outputC:

  Output Hierarchical Binning (ABxCB) Clusters files location and
  prefix.

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
[`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html):

- `read_id` : read identifier from fasta header

- `ABxCB` : read was assigned to this ABxCB cluster index

- `ABxCB.n` : read to cluster ABxCB.n distance

For backwards-compatible `data.frame` output use
`as.data.frame(result)`.

## References

<https://github.com/gkanogiannis/metabinR>

## Author

Anestis Gkanogiannis, <anestis@gkanogiannis.com>

## Examples

``` r
res <- hierarchical_binning(
    system.file("extdata", "reads.metagenome.fasta.gz", package = "metabinR"),
    dryRun = TRUE, kMerSizeAB = 4, kMerSizeCB = 2
)
res
#> MetabinResult (ABxCB)
#>   reads:     26,664
#>   clusters:  1
#>   inputs:    1 file(s)
#>     - /home/runner/work/_temp/Library/metabinR/extdata/reads.metagenome.fasta.gz
```
