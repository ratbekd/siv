# siv
Synthetic Instrumental Variable Method R-package

The package estimates a linear regressions with one endogeneous variable using three different techniques of the SIV method.The essense of the mthod is in construction of a synthetic IV using  $s = x + k \delta r$, where $x$ is the endogenous variable in the reduced form, and $r$ is a vector determined in the plane spanned by the outcome variable $y$ and endogenous variable $x$. Vector $r$ is constructed as an orthogonal vector to $x$. 
The first method is a simple SIV method that assumes homoscedasdicity of the error term. This method implies that if we synthesize such an  SIV,  s*, that satisfies $E(s*'e'e) = 0$, then, $E(s*'u) = 0$ also must hold. 
Thus, in this case, 
the SIV method determines a valid SIV such that $E(s|u) = 0$   by $s*=x+k\delta_0 r$ where $\delta_0=arg[E(ee'| s*)=0]$, and $e$ is the first-stage error term. Then, $\beta$, the  parameter in  the regression equaition, is identified by an IV estimator: 
$\hat{\beta}_{IV}=(x's^*)^{-1} x'y.$

In the robust to heteroscedasticity approach, we use the difference in the degree of the heteroscedasticity, $\Delta$, is estimated by  prametrically or non-parametrically and their locus over $\delta \in (0, \bar{\delta})$ as given by  function $D_{\Delta}$.
Then,  a valid SIV such that $E(  u| s*)=0$ is identified by $s*= x+k\delta_0  r$  where $\delta_0 =argmin_{\delta}(  D_{\Delta})$.   The details can be found in the related paper at  https://github.com/ratbekd/SIV/blob/main/SIV_DT_R7.pdf.

The package can be installed in Rstudio paltform using this command:
remotes::install_git("https://github.com/ratbekd/siv.git")
library(siv)
## Example based on Mroz data
data <- wooldridge::mroz  # Use sample data set

data <- data[complete.cases(data), ]  # Remove missing values

attach(data)
# Run regression
#Y="hours" # outcome variable
#X="lwage"# endogenous variable
#H=c("educ", "age", "kidslt6", "kidsge6", "nwifeinc")# exogenous variables
result <- siv_regression(data, "hours", "lwage", c("educ", "age", "kidslt6", "kidsge6", "nwifeinc"), reps=5)

iv1 <- (result$IV1)

iv2 <-(result$IV2)

iv3 <-(result$IV3)

summ.iv1 <- summary(iv1, diagnostics=T)

summ.iv2 <- summary(iv2, diagnostics=T)

summ.iv3 <- summary(iv3, diagnostics=T)
One can review the $delta_0$ values found using the different appoaches:

$result\$$ $delta_0$
 The first0satge equations can be obtained by using result$FS3.

