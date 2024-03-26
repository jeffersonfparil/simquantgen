# simquantgen

Simulate genotype and phenotype data for quantitative genetics analyses in R.

|**Build Status**|**License**|
|:--------------:|:---------:|
| <a href="https://github.com/jeffersonfparil/simquantgen/actions"><img src="https://github.com/jeffersonfparil/simquantgen/actions/workflows/r.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

## Installation

```R
devtools::install_github("jeffersonfparil/simquantgen")
```

## Usage

```R
library(simquantgen)
?fn_simulate_genotypes
?fn_simulate_phenotypes
?fn_simulate_gxe
?fn_simulate_sparse_mixed_model_data
```

## Examples

```R
library(simquantgen)
G = fn_simulate_genotypes(n=120, l=512, ploidy=42, n_alleles=4, verbose=TRUE)
list_Y_b_E_b_epi = fn_simulate_phenotypes(G=G, n_alleles=4, dist_effects="chi2", n_effects=25, h2=0.75, pheno_reps=5, verbose=TRUE)
list_df_CORR = fn_simulate_gxe(G=G, n_effects=50, purely_additive=FALSE, n_networks=10, n_effects_per_network=50, h2=0.5, env_factor_levels=c(5, 3), env_factor_effects_sd=0.2, frac_additional_QTL_per_env=0.15, n_reps=5, verbose=TRUE)
list_y_complete_y_X_Z_D = fn_simulate_sparse_mixed_model_data(G=G, df_gxe=list_df_CORR$df, frac_gen_missing=0.25)
```


## Developer Notes

Clone and navigate into the repository. Open R (>=4.3.1) and run:

```R
install.packages("devtools")
devtools::load_all()
devtools::document()
devtools::test()
```
