# simquantgen
Simulate genotype and phenotype data for quantitative genetics analyses

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

## Tests

```R
devtools::load_all()
devtools::document()
devtools::test()
```
