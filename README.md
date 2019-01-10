# Parainflammation


# Install

devtools::install_github('dviraran/Parainflammation')

# Usage

library(Parainflammation)

To calculate the PI scores:

```R
pi = calculatePIscore(expr,tissue)
```

tissue options can be found in fit_params$tissues

To adjust gene expression to CD45 expression:

```R
adjustExpressionToCD45(expr,params)
```

params is a the column in fit_params$coef corresponding to the the associated tissue. For example, if the tissue is Kidney, then: 

```R
adjustExpressionToCD45(expr,fit_params$coef[,fit_params$tissue=='Kidney'])
```

