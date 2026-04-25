# JVM options used when metabinR loads

Returns the vector of JVM flags passed to
[`.jpackage`](https://rdrr.io/pkg/rJava/man/jpackage.html) on package
load. Set `options(metabinR.jvm.flags = c(...))` before loading the
package to override; set `options(java.parameters = ...)` to prepend
heap-size flags (e.g. `"-Xmx4g"`) in the usual rJava way.

## Usage

``` r
metabinR_jvm_options()
```

## Value

A character vector of JVM flags.

## Examples

``` r
metabinR_jvm_options()
#> [1] "-Djava.awt.headless=true"    "-XX:+UseG1GC"               
#> [3] "-XX:+UseStringDeduplication"
```
