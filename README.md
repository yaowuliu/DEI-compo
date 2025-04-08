# DEIcompo
This package concerns the direct estimation of the type I error methods for composite null testing. In particular, it implemented the DEI_B, a fast and powerful method designed for large-scale composite null testing, with application to mediation analysis.

## Installation
```
library(devtools)
devtools::install_github("yaowuliu/DEI-compo")
```

## Usage
The code below provides an example about the usage of the DEI_B method.

```
set.seed(123)

# number of items (e.g., mediators)
N = 100000  

### Generate the individual p-values from z scores. The first 100 items are under the composite alternative
x = rnorm(N); y = rnorm(N)
x[1:100] = x[1:100] + runif(100,min=4,max = 4)
y[1:100] = y[1:100] + runif(100,min=4,max = 4)
x[101:200]  = x[101:200] + runif(100,min=4,max = 4)
y[201:300]  = y[201:300] + runif(100,min=4,max = 4)
px = 2*pnorm(-abs(x)); py = 2*pnorm(-abs(y))

### Compute the p-values for the composite null.
pvals_compo = DEI_B(px,py)

### Apply the BH procedure to control FDR. Note there are many NA's in pvals_compo
pvals_compo[is.na(pvals_compo)] = 1
which(p.adjust(pvals_compo, method = "BH") < 0.05)
```
