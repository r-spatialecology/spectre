# calc_commonness_error

Calculate commonness error

## Usage

``` r
calc_commonness_error(x, objective_matrix)
```

## Arguments

- x:

  Results object from run_optimization_min_conf.

- objective_matrix:

  Matrix from (modeled) alpha-diversity and Bray-Curtis dissimilarity

## Value

vector

## Details

Calculate mean absolute commonness error (MAE_c) and relative commonness
error in percentage (RCE).
