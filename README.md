
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RecPD

<!-- badges: start -->

<!-- badges: end -->

The RecPD package was developed to calculate the recombination-adjusted
phylogenetic diversity of a given feature (phenotype, variant, etc. in
binary presence/absence format) found across a given species
phylogenetic tree utilizing ancestral state reconstruction.

RecPD also provides functions for:

  - Visualizing feature ancestral lineage reconstructions.
  - Derived measures such as: **nRecPD** (a measure of the degree of
    feature recombination), **Span**, **Clustering**, **Longevity**, and
    **Lability**.
  - **RecPDcor** - a measure of pairwise feature lineage correlation.

For a more detailed description of the RecPD methodology, and its
advantage over using simple Prevalence to measure feature diversity in
microbial population genetics studies, see
[Bundalovic-Torma C. & Guttman D. PLoS CompBiol.
(2022)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009899).

## Installation

You can install the development version of RecPD from
[GitHub](https://github.com/)

``` r
# install.packages("devtools")
devtools::install_github("cedatorma/recpd")
```

## Tutorial and Usage Examples

See the ./vignettes directory for a tutorial providing a general overview of RecPD and the various functionalities offered in the package.
