# metabinR: Abundance and Compositional Based Binning of Metagenomes

Provide functions for performing abundance and compositional based
binning on metagenomic samples, directly from FASTA or FASTQ files.
Functions are implemented in Java and called via rJava. Parallel
implementation that operates directly on input FASTA/FASTQ files for
fast execution. Inputs may be file paths or Biostrings/ShortRead
sequence objects; results are returned as a MetabinResult S4 object
wrapping cluster assignments, algorithm parameters, and input metadata.

## See also

Useful links:

- <https://github.com/gkanogiannis/metabinR>

- Report bugs at <https://github.com/gkanogiannis/metabinR/issues>

## Author

**Maintainer**: Anestis Gkanogiannis <anestis@gkanogiannis.com>
([ORCID](https://orcid.org/0000-0002-6441-0688))
