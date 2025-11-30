# generate_commonness_matrix_from_gdm

Creates a pairwise site by site commonness matrix from estimates of
species richness and Bray-Curtis dissimilarity.

## Usage

``` r
generate_commonness_matrix_from_gdm(gdm_predictions, alpha_list)
```

## Arguments

- gdm_predictions:

  a square pairwise `matrix` of Bray-Curtis dissimilarity estimates
  between site pairs. We recommend using the gdm-package (Fitzpatrick et
  al. 2020) to generate this matrix

- alpha_list:

  a `vector` of species richness for every site in the study area. The
  length of this vector must be equivalent to one of the dimensions of
  the `gdm_predictions`

## Value

A pairwise site by site `matrix` of the number of species in common
between each site pair, with dimensions equal to that of the provided
dissimilarity matrix.

## Details

`generate_commonness_matrix_from_gdm` uses a vector of estimated species
richness per site and a pairwise matrix of site by site Bray-Curtis
dissimilarity (we recommend using the gdm-package (Fitzpatrick et al.
2020) to generate this matrix) to produce a matrix of the estimated
species in common between site pairs (referred to as a commonness
matrix). The commonness between sites is calculated using
\$\$C\_{ij}=(1-\beta\_{ij})(S\_{i} + S\_{j})/2\$\$ Where \\\beta\_{ij}\\
is the dissimilarity between sites, \\C\_{ij}\\ is the species in common
between sites, and S is the number of species in each site. For more
details see Mokany et al 2011.

## References

Mokany, K., Harwood, T.D., Overton, J.M., Barker, G.M., & Ferrier, S.
(2011). Combining \\\alpha\\ and \\\beta\\ diversity models to fill gaps
in our knowledge of biodiversity. Ecology Letters, 14(10), 1043-1051.
