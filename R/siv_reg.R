#' Perform Instrumental Variable Regression using SIV
#'
#' @param data A data frame.
#' @param Y Name of the dependent variable.
#' @param X Name of the endogenous variable.
#' @param H Vector of exogenous variable names.
#' @param reps Number of bootstrap loops.
#' @return A list containing SIV regression results.IV1- a simple SIV, IV2- a parametric SIV, IV3- a non-parametric SIV, citable- a table of confidence intervals.
#' @export
#' @importFrom stats  predict
#' @importFrom stats  complete.cases
#' @importFrom stats  fitted
#' @importFrom stats pf
#' @importFrom stats pchisq
#' @importFrom stats sum
#' @importFrom stats sumsq
#' @importFrom stats var
#' @importFrom stats as.formula
#' @importFrom stats lm
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats residuals
#' @importFrom stats df.residual
#' @importFrom stats qt
#' @importFrom stats sign
#' @importFrom stats rnrom
#' @importFrom stats sd
#' @importFrom stats mean
#' @importFrom stats abs#'
#' @examples
#' df <- wooldridge::mroz  # Use sample data set
#' data <- df[complete.cases(df), ]  # Remove missing values
#' attach(data)
#' H0 <- data.frame(hours, lwage, educ, age, kidslt6, kidsge6, nwifeinc)
# Vector of variables used in the regression. Firts var is outcome, the second var is endogenous regressor
#' Y <- as.character(colnames(H0))[1] ###OUTCOME variable
#' X <-as.character(colnames(H0))[2] ###ENDOEGENOUS variable
#' H<- as.character(colnames(H0))[-(1:2)] #
#' result <- siv_reg(data, Y, X, H, reps=5)
#' iv1 <-(result$IV1)# a simple SIV
#' iv2 <-(result$IV2)# a robust parametric SIV (RSIV-p)
#' iv3 <-(result$IV3)# a robust non-parametric SIV (RSIV-n)
#' summ.iv1 <- summary(iv1, diagnostics=TRUE)
#' summ.iv2<- summary(iv2, diagnostics=TRUE)
#' summ.iv3 <- summary(iv3, diagnostics=TRUE)
#'result$citable ###  renders the CI table of beta estimates
#'result$delta0 ###  renders the delta0 estimates:SIV, RSIV-p, RSIV-n.
#'SIv <- result$sivs[[2]] ## renders the robust SIV estimation
siv_reg <- function(data, Y, X, H, reps) {
  library(dplyr)
  library(ivreg)
  library(lmtest)
  library(sandwich)
  library(rMR)
  library(AER)
  library()
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  deg2rad <- function(deg) {(deg * pi) / (180)}
  # Ensure variables exist in data
  if (!all(c(Y, X, H) %in% names(data))) {
    stop("Some variables are missing in the dataset.")
  }
  # Y <- as.character(colnames(H0))[1] ###OUTCOME variable
  # X <-as.character(colnames(H0))[2] ###ENDOEGENOUS variable
  # H<- as.character(colnames(H0))[-(1:2)] ##### EXOGENOUS variables
  #print(H)
  formula_str <- paste(paste0(Y," ~ ", X,"+"), paste(H, collapse = " + "))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object
  #print(formula)

  # determine the number of rows
  N<-nrow(data)
  #####SIV Method
  ###Basic vectors for  SIV calculation###
  # ### the outcome variable
  ## Factoring out the effects of other exogenous variables
  formula_str <- paste(paste0(Y," ~ "), paste(H, collapse = " + "))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object
  fity<-lm(formula, data=data)
  y<-resid(fity)
  # the endogenous variable
  ## Factoring out the effects of other exogenous variables
  formula_str <- paste(paste0(X," ~ "), paste(H, collapse = " + "))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object
  fitx<-lm(formula, data=data)
  x<-resid(fitx)
  xx <- min(length((x-mean(x))),length((y-mean(y))))
  data <- data[1:xx,]
  #saving the transformed x and y
  data$x<-(x-mean(x))
  data$y<-(y-mean(y))
  y0<-y
  x0<-x
  V=0
  ### Generating a vector orthogonal to x
  fity<-lm(y0~(x0), data=data)
  V<-resid(fity)
  V<-(V-mean(V))/sd(V)
  V<-V*sd(x0)
  data$R<-V
  theta=0
  theta1=0
  d<-0.01# starting value for delta
  delt<-0.01# step to change delta
  while(theta1<70){
    data$siv<-(data$x-d*data$R) ### SIV
    theta <- acos( sum(data$x*data$siv) / ( sqrt(sum(data$x^2)) * sqrt(sum(data$siv^2)) ) )
    theta1<-rad2deg(theta)
    d=d+delt
  }
  dd <- d
  dl=round(dd/delt)
  #####################################################################
  ############ Determining the sign of cor(x,u)
  k=0
  j=1
  data$siv=0
  signc=matrix(ncol = 5, nrow = 2)
  signc[1,1] <-1
  signc[2,1] <--1
  #

  for (j in 1:2) {
    if(j<2){k=1}else{k=-1}#the assumed sign for cor(x,u)

    # IV regression
    #dd<-3#end value for delta
    d<-0.01# starting value for delta
    delt<-0.01# step to change delta
    i<-1 ### starting value for a counter
    ## placeholders for variables
    # m1=0
    # m2=0
    dl=round(dd/delt)-1
    m1 <- numeric(length = dl)  # Pre-allocate memory
    m2 <- numeric(length = dl)
    while (d<dd){# we compute m1 until siv is close to become perpendicular to x
      data$siv<-(data$x-k*d*data$R) ### SIV
      ####OLS estimates
      rls<-(lm(data$x~data$siv, data=data))
      s.rls<-summary(rls,vcov. = function(x) vcovHC(x, type="HC1"), diagnostics=T)
      data$ev21<-resid(s.rls)
      m1[i]<-(cov(data$ev21^2,data$siv))
      m2[i] <- (cor(data$ev21^2,data$siv))
      theta <- acos( sum(data$x*data$siv) / ( sqrt(sum(data$x^2)) * sqrt(sum(data$siv^2)) ) )
      theta1<-rad2deg(theta)
      d<-d+delt
      i<-i+1
    }

    signc[j,3] <- check_initial_abs_increase(m1)#
    signc[j,5] <- check_initial_abs_increase(m2)
    if(signc[j,3]!=1){
      index <- check_sign_change(m1)
      m=m1[index:i]
    }else{m=m1}
    signc[j,2] <-check_sign_change(m1)
    signc[j,3] <- check_initial_abs_increase(m1)#
    signc[j,4] <- check_sign_change(m2)
  }

  ch=0
  for(j in 1:2){
    ch[j] <- signc[j,2]+signc[j,3]+signc[j,4]+signc[j,5]
  }

  k <- signc[which.max(ch),1]
 # cat("The true sign for cor(xu) is", k )
  #####################################

  if(k!=0){
    ############# Initial settings for SIV
    vvar=0
    d0i=0
    d0ri=0
    d0rni=0
    b2=0
    b2r=0
    b2rp=0
    b2t=0
    N<-nrow(data)
    reps=reps
    S=round(N*.999)
    l=1
    fitc <- matrix(ncol = 1, nrow = reps)
    sumb2=matrix(ncol = 1, nrow = reps)
    fitcr <- matrix(ncol = 1, nrow = reps)
    sumb2r=matrix(ncol = 1, nrow = reps)
    fitcrn <- matrix(ncol = 1, nrow = reps)
    sumb2rn=matrix(ncol = 1, nrow = reps)
    fitct <- matrix(ncol = 1, nrow = reps)
    sumb2t=matrix(ncol = 1, nrow = reps)
    lowbp=0
    upbp=0

    ##### Bootstrap sampling loop. You may use data <- mydata instead of data <- mydata[sample(1:N, S),  TRUE]
    #if you just want see how it works for the original sample data.
    while (l<reps){
      # IV regression
      set.seed(3*(l))   # a different seed for each sub sample
      # data <- data[sample(1:N, S),  TRUE]
      data_sample <- data[sample(1:N, S, replace = TRUE), ]
      ####Computation of m1

      #dd<-3#end value for delta
      d<-0.01# starting value for delta
      delt<-0.01# step to change delta
      i<-1 ### starting value for a counter
      ## placeholders for variables
      m1=numeric(length = dl)  # Pre-allocate memory
      st=0
      ev22=0
      dv=0
      dv2=numeric(length = dl)
      x4=numeric(length = dl)
      l1=0
      l2=0
      while (d<dd){# we compute m1 until siv is close to become perpendicular to x
        data_sample$siv<-(data_sample$x-k*d*data_sample$R)
        #
        rls<-(lm(x~siv, data=data_sample))
        s.rls<-summary(rls, diagnostics=T)
        data_sample$siv<-(data_sample$x-k*d*data_sample$R)
        s.rls<-summary(rls,vcov. = function(x) vcovHC(x, type="HC1"), diagnostics=T)
        data_sample$ev21<-resid(s.rls)
        ##FGLS
        ehatsq <- resid(rls)^2
        sighatsq.rls  <- lm(log(ehatsq)~siv, data=data_sample)
        data_sample$vari <- sqrt(exp(fitted(sighatsq.rls)))
        vvar[i] <- var(data_sample$vari)
        fgls <- lm(x~siv, weights=1/vari, data=data_sample)
        data_sample$ev22<-resid(fgls)
        m1[i]<-cor(data_sample$ev21^2,data_sample$siv)#
        l1 <- summary(lm((data_sample$ev21^2)~data_sample$siv,data=data_sample))
        l2 <-summary( lm((data_sample$ev22^2)~data_sample$siv,data=data_sample))
        n=length(data_sample$ev21)
        ssr1 <- sumsq(predict(lm((ev21^2)~siv,data=data_sample))-mean(data_sample$ev21^2))
        sse1=sumsq(data_sample$ev21)
        x1= (ssr1/2)/(sse1/n^2)^2
        ssr2 <- sumsq(predict(lm((ev22^2)~siv,data=data_sample))-mean(data_sample$ev22^2))
        sse2=sumsq(data_sample$ev22)
        x2= (ssr2/2)/(sse2/n^2)^2
        #
        x3 <- x1/x2#
        dv2[i] <- pf(x3, df = 1, df2 = 1, lower.tail = TRUE)
        st[i]<-d
        # non-parametric
        samp1 <- (predict(lm((data_sample$ev21^2)~data_sample$siv,data=data_sample)))^2
        samp2 <- (predict(lm((data_sample$ev22^2)~data_sample$siv,data=data_sample)))^2
        xx0 <- samp1
        yy0 <- samp2
        ad0 <- ad2_stat(x = xx0, y = yy0)
        x4[i] <- 1-ad0
        d<-d+delt
        i<-i+1
      }
      ### updating the formula for regressions
      formula_str <- paste(paste0(Y," ~ ", X,"+"), paste(H, collapse = " + "))## Construct formula as a string
      formula <- as.formula(formula_str)# Convert to formula object
      ### update with your own instruments: instruments<-~ siv+all exogenous variables
      iv_str <- paste(paste0(" ~ ", "siv","+"), paste(H, collapse = " + "))## Construct formula as a string
      instruments <- as.formula(iv_str)# Convert to formula object
      #instruments<-~ siv+all exogenous variables
      # formula <-hours~lwage+educ+ age+kidslt6+ kidsge6+ nwifeinc
      ###the First-Stage equation formula
      fs_formula_str <- paste(paste0(X," ~ ","siv","+"), paste(H, collapse = " + "))#
      #### DT condition of homoscedatic case
      d0 <- (which.min(abs(m1)))*delt
      d0i[l] <- d0
      data_sample$siv<-(data_sample$x-k*d0*data_sample$R)
      iv2<-ivreg(formula, instruments, data=data_sample)
      summ.iv2 <- summary(iv2, diagnostics=T)#
      ## saving the estimation parameters for each sample
      fitc[l] <- iv2$coefficients[2]
      sumb2[l] <-  summ.iv2$coefficients[2,2]

      ### DT point for heteroscedastic case- parametric approach
      d0r <-  which.min(dv2)*delt
      d0ri[l] <- d0r
      data_sample$siv<-(data_sample$x-k*d0r*data_sample$R)
      iv3<-ivreg(formula, instruments, data=data_sample)
      summ.iv3 <- summary(iv3, diagnostics=T)#,
      ## saving the estimation paramters for each sample
      fitcr[l] <- iv3$coefficients[2]
      sumb2r[l] <-  summ.iv3$coefficients[2,2]

      #### DT point for heteroscedastic case- non-parametric approach
      d0rn <- which.min(x4)*delt
      d0rni[l] <- d0rn
      data_sample$siv<-(data_sample$x-k*d0rn*data_sample$R)
      iv4<-ivreg(formula, instruments, data=data_sample)
      summ.iv4 <- summary(iv4, diagnostics=T)#,
      ## saving the estimation parameters for each sample
      fitcrn[l] <- iv4$coefficients[2]
      sumb2rn[l] <-  summ.iv4$coefficients[2,2]
      l <- l+1
    }
    sign_cor_ux <- k
    ## Distribution of sample paramters
    # the simple homogenous case
    fitc <- fitc[complete.cases(fitc)]
    sumb2<- sumb2[complete.cases(sumb2)]
    alpha <- 0.05 # chosen significance level
    b2[1] <- mean(fitc)# the endogenous parameter beta
    df <- df.residual(iv2)# degrees of freedom for t-stat
    seb2 <-mean(sumb2)# the average standard error of the parameter beta
    tc <- qt(1-alpha/2, df) ## t-statistic
    lowbp[1] <- b2[1]-tc*seb2  # lower bound for beta
    upbp[1] <- b2[1]+tc*seb2   # upper bound for beta

    #### The parametric heterogenous case
    fitcr <- fitcr[complete.cases(fitcr)]
    sumb2r<- sumb2r[complete.cases(sumb2r)]
    alpha <- 0.05 # chosen significance level
    b2[2] <- mean(fitcr)
    df <- df.residual(iv3)
    seb2r <-mean(sumb2r)
    tc <- qt(1-alpha/2, df)
    lowbp[2] <- b2[2]-tc*seb2r  # lower bound
    upbp[2] <- b2[2]+tc*seb2r   # upper bound

    #################The non-parametric heterogenous case
    fitcrn <- fitcrn[complete.cases(fitcrn)]
    sumb2rn<- sumb2rn[complete.cases(sumb2rn)]
    alpha <- 0.05 # chosen significance level
    b2[3] <- mean(fitcrn)
    df <- df.residual(iv4)
    seb2rn <-mean(sumb2rn)
    tc <- qt(1-alpha/2, df)
    lowbp[3] <- b2[3]-tc*seb2rn  # lower bound
    upbp[3] <- b2[3]+tc*seb2rn   # upper bound

    ### Table for CI of beta
    mv<-data.frame(lowbp,b2,upbp)
    colnames(mv)<-c("low beta","mean b2", "high beta")
    rownames(mv)<- c("SIV","SIVRr","SIVRn")#, "nearc4")
    ###   (mv)

    ### final satge estimations
    ###### Simple homogenous assumption case
    d0i <-  d0i[complete.cases(d0i)]
    d0m <- mean(d0i)
    data$siv<-(data$x-k*d0m*data$R)
    siv<-(data$x-k*d0m*data$R)
    fs1<-lm(as.formula(fs_formula_str), data=data)
    iv1<-ivreg(formula, instruments, data=data)
    #summ.iv1 <- summary(iv1, diagnostics=T)#
    ############Paramteric heterogenous case
    d0ri <-  d0ri[complete.cases(d0ri)]
    d0rm <- mean(d0ri)
    data$siv<-(data$x-k*d0rm*data$R)
    sivr<-(data$x-k*d0rm*data$R)
    fs2<-lm(as.formula(fs_formula_str), data=data)
    iv2<-ivreg(formula, instruments, data=data)
    #

    #### Non-parameteric heterogenous case
    d0rni <-  d0rni[complete.cases(d0rni)]
    d0rnm <- mean(d0rni)
    data$siv<-(data$x-k*d0rnm*data$R)
    sivrn<-(data$x-k*d0rnm*data$R)
    fs3<-lm(as.formula(fs_formula_str), data=data)
    iv3<-ivreg(formula, instruments, data=data)
    #
      }else{
    v=rnorm(N,0,1)
    print("NO endogeneity problem. All SIV estimates are the same as the OLS")
    d0m=0.001
    data$siv<-(data$x-k*d0m*data$R)+v
    iv1<-ivreg(formula, instruments, data=data)
    #
    data$siv<-(data$x-k*d0rm*data$R)+v
    iv2<-ivreg(formula, instruments, data=data)
    #
    data$siv<-(data$x-k*d0rnm*data$R)+v
    iv3<-ivreg(formula, instruments, data=data)

    siv<-(data$x-k*d0m*data$R)
    ############Parametric heterogeneous case
    sivr <-siv
    #### Non-parametric heterogeneous case
    sivrn<-siv
    #
      }
  ### Table for CI of beta
  siv_tb<-data.frame(siv,sivr,sivrn)
  colnames(siv_tb)<-c("SIV","SIV_r", "SIV_rn")
  #rownames(mv)<- c("SIV","SIVRr","SIVRn")#, "nearc4")
  ###   (mv)
 # cat("Use result$IV1, result$IV3, and result$IV3 for SIV estimations.  Use result$delta0 to  obtain delta0 to compute the SIV" )
  return(list(IV1 = (iv1), IV2 = (iv2), IV3 = (iv3), FS1=(fs1), FS2=(fs2), FS3=(fs3), sign=(sign_cor_ux), sivs=(siv_tb), citable=mv, delta0=c(d0m, d0rm, d0rnm)))
}
