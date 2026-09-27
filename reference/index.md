# Package index

## Binning functions

High-level entry points for metagenome binning.

- [`abundance_based_binning()`](https://gkanogiannis.github.io/metabinR/reference/abundance_based_binning.md)
  : Abundance based binning on metagenomic samples
- [`composition_based_binning()`](https://gkanogiannis.github.io/metabinR/reference/composition_based_binning.md)
  : Composition based binning on metagenomic samples
- [`hierarchical_binning()`](https://gkanogiannis.github.io/metabinR/reference/hierarchical_binning.md)
  : Hierarchical (ABxCB) binning on metagenomic samples

## Result analysis

Summarize assignments and compare them with known read origins.

- [`bin_summary()`](https://gkanogiannis.github.io/metabinR/reference/bin_summary.md)
  : Summarize the bins in a metabinR result
- [`ambiguous_reads()`](https://gkanogiannis.github.io/metabinR/reference/ambiguous_reads.md)
  : Find reads with close competing bins
- [`evaluate_bins()`](https://gkanogiannis.github.io/metabinR/reference/evaluate_bins.md)
  : Evaluate read-level bin assignments against known origins

## MetabinResult

S4 container for binning assignments and parameters.

- [`show(`*`<MetabinResult>`*`)`](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-class.md)
  : MetabinResult: binning result container

- [`assignments()`](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-accessors.md)
  [`nClusters()`](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-accessors.md)
  [`parameters()`](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-accessors.md)
  [`algorithm()`](https://gkanogiannis.github.io/metabinR/reference/MetabinResult-accessors.md)
  : Accessors for MetabinResult

- [`as.data.frame(`*`<MetabinResult>`*`)`](https://gkanogiannis.github.io/metabinR/reference/as.data.frame-MetabinResult-method.md)
  :

  Convert a MetabinResult to a `data.frame`

## JVM configuration

- [`metabinR_jvm_options()`](https://gkanogiannis.github.io/metabinR/reference/metabinR_jvm_options.md)
  : JVM options used when metabinR loads
