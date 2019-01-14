# Parainflammation score


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

# Citation

Dvir Aran, Audrey Lasry, Adar Zinger, Moshe Biton, Eli Pikarsky, Asaf Hellman, Atul J. Butte and Yinon Ben-Neriah. Widespread parainflammation in human cancer. Genome Biology 2016 17:145

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0995-z
