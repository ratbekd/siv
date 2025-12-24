# siv
Synthetic Instrumental Variable Method R-package

The package estimates a linear regression with one endogenous variable using three different techniques of the SIV method. The essence of the process is in the construction of a synthetic IV using  $s = x + k \delta r$, where $x$ is the endogenous variable in the reduced form, and $r$ is a vector determined in the plane spanned by the outcome variable $y$ and endogenous variable $x$. Vector $r$ is constructed as an orthogonal vector to $x$. 
The first method is a simple SIV method that assumes homoscedasticity of the error term. This method implies that if we synthesise such an  SIV,  s*, that satisfies $cov(s*, e^2) = 0$, then, $E(s*'u) = 0$ also must hold. 
Thus, in this case, 
the SIV method determines a valid SIV such that $E(s|u) = 0$   by $s*=x+k\delta_0 r$ where $\delta_0=arg[M(\delta)=0]$. Here $M(\delta)$ is a covariance function of $e$, the first-stage error term. This condition is used for the method "GMM". To save time, one can use  a simpler method, where one uses a $\delta_0=arg[cov(e^2(\delta), s(\delta))=0]$. 
Then, $\beta$, the  parameter in  the regression equation, is identified by an IV estimator: 
$\hat{\beta}_{IV}=(x's^*)^{-1} x'y.$

In the robust to heteroscedasticity approach, we use the difference in the degree of the heteroscedasticity, $D$, which is estimated by  parametrically or non-parametrically, and its locus over $\delta \in (0, \bar{\delta})$ as given by  the function $D_{E}$.
Then,  a valid SIV such that $E( u| s*)=0$ is identified by $s*= x+k\delta_0  r$  where $\delta_0 =argmin_{\delta}(  D_{E})$.   The details can be found in the related paper at  [SIV_DT.pdf](https://github.com/ratbekd/SIV-method/blob/main/SIV_DT_R13.pdf).

The package can be installed in the RStudio platform using this command:
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
result <- siv_reg(data, "hours", "lwage", c("educ", "age", "kidslt6", "kidsge6", "nwifeinc"), method="simple", reps=5)

iv1 <- (result$IV1)

iv2 <-(result$IV2)

iv3 <-(result$IV3)

summ.iv1 <- summary(iv1, diagnostics=T)

summ.iv2 <- summary(iv2, diagnostics=T)

summ.iv3 <- summary(iv3, diagnostics=T)
One can review the $delta_0$ values found using the different approaches:

$result\$$ $delta_0$
 The first-stage equations can be obtained by using result$FS3.

