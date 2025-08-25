#' Perform Instrumental Variable Regression using multiple SIVs
#'
#' @param data A data frame.
#' @param Y Name of the dependent variable.
#' @param E Name of the endogenous variables.
#' @param H Vector of exogenous variable names
#' @param reps Number of bootsrap loops.
#' @return A list containing OLS and IV regression results.
#' @export
#' @examples
#' df <- wooldridge::mroz  # Use sample data set
#' data <- df[complete.cases(df), ]  # Remove missing values
#' attach(data)
#' result <- msiv_reg(data, "hours", c("lwage", "educ"),c( "age", "kidslt6", "kidsge6", "nwifeinc"), reps=5)
msiv_reg <- function(data, Y, E, H, reps) {
  library(dplyr)
  library(ivreg)
  library(lmtest)
  library(sandwich)
  library(AER)
  library()
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  deg2rad <- function(deg) {(deg * pi) / (180)}
  #Two-sample Anderson-Darling statistic
  ad2_stat <- function(x, y) {
    # Sample sizes
    n <- length(x)
    m <- length(y)

    # Pooled sample and pooled ecdf
    z <- c(x, y)
    z <- z[-which.max(z)] # Exclude the largest point
    H <- rank(z) / (n + m)

    # Statistic computation via ecdf()
    (n * m / (n + m)^2) * sum((ecdf(x)(z) - ecdf(y)(z))^2 / ((1 - H) * H))

  }

  check_sign_change <- function(x) {
    # Step 2: Check if there is a sign change
    signs <- sign(x)  # Get sign (-1, 0, or 1)
    sign_changes <- any(diff(signs) != 0, na.rm = TRUE)  # Detect if sign changes

    # Step 3: Return 1 if both conditions are met, otherwise 0
    return(as.integer(sign_changes))
  }


  check_initial_abs_increase <- function(x) {
    # Step 1: Check if the absolute values start by increasing (non-decreasing trend)
    abs_x <- abs(x)
    initial_increasing <- all(diff(abs_x[1:min(20, length(abs_x))]) >= 0)  # Check first 3 points or as many as available

    # Step 3: Return 1 if both conditions are met, otherwise 0
    return(as.integer(initial_increasing))
  }
  find_first_sign_change <- function(x) {
    sign_changes <- which(diff(sign(x)) != 0)  # Find indices where sign changes
    if (length(sign_changes) > 0) {
      return(sign_changes[1] + 1)  # Return the first occurrence (adjust for diff)
    } else {
      return(NA)  # Return NA if no sign change
    }
  }
  #Ensure variables exist in data
   if (!all(c(Y, E, H) %in% names(data))) {
     stop("Some variables are missing in the dataset.")
   }
  # Y <- as.character(colnames(H0))[1] ###OUTCOME variable
  # X <-as.character(colnames(H0))[2] ###ENDOEGENOUS variable
  # H<- as.character(colnames(H0))[-(1:2)] ##### EXOGENOUS variables
  if (reps>10){reps=10} else{reps=reps}
  # determine the number of rows
  N<-nrow(data)
  zk=length(E)
  N=nrow(data)
  siv_list <- list()
  signk=0
  #####loop for each endogenous variable
for (g in 1:length(E)) {
  X <- E[g]

  #####SIV Method

  ###Basic vectors for  SIV calculation###
  ## Factoring out the effects of other exogenous variables
  if(length(E)>1){
  formula_str <- paste(paste0(Y," ~ "), paste(H, collapse = " + "),paste("+",E[-g],collapse = "+"))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object

  fity<-lm(formula, data=data)
  y<-resid(fity)
  ## Factoring out the effects of other exogenous variables
  formula_str <- paste(paste0(X," ~ "), paste(H, collapse = " + "),paste("+",E[-g],collapse = "+"))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object
  fitx<-lm(formula, data=data)
  x<-resid(fitx)
  }else{
    formula_str <- paste(paste0(Y," ~ "), paste(H, collapse = " + "))## Construct formula as a string
    formula <- as.formula(formula_str)# Convert to formula object
    fity<-lm(formula, data=data)
    y<-resid(fity)
    #x1<-lwage### the endogenous variable
    ## Factoring out the effects of other exogenous variables
    formula_str <- paste(paste0(X," ~ "), paste(H, collapse = " + "))## Construct formula as a string
    formula <- as.formula(formula_str)# Convert to formula object
    fitx<-lm(formula, data=data)
    x<-resid(fitx)
  }
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

  #####################################################################
  ############ Determining the sign of cor(x,u)
  k=0
  j=1

  signc=matrix(ncol = 5, nrow = 2)
  signc[1,1] <-1
  signc[2,1] <--1
  #

  for (j in 1:2) {
    if(j<2){k=1}else{k=-1}#the assumed sign for cor(x,u)

    # IV regression
    theta=0
    theta1=0
    dd<-3#end value for delta
    d<-0.01# starting value for delta
    delt<-0.01# step to change delta
    i<-1 ### starting value for a counter
    ## placeholders for variables

    dl=round(dd/delt)
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
    par(mfrow = c(1, 2))

    plot(m1[0:300])
    signc[j,3] <- check_initial_abs_increase(m1)#result
    signc[j,5] <- check_initial_abs_increase(m2)
    if(signc[j,3]!=1){
      index <- find_first_sign_change(m1)
      m=m1[index:i]
    }else{m=m1}
    signc[j,2] <-check_sign_change(m1)
    signc[j,3] <- check_initial_abs_increase(m1)#result
    signc[j,4] <- check_sign_change(m2)
  }

  ch=0
  for(j in 1:2){
    ch[j] <- signc[j,2]+signc[j,3]+signc[j,4]+signc[j,5]
  }


  k <- signc[which.max(ch),1]
  signk[g] <- k
  #cat("The true sign for cor(xu) is", k )
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
   #
    if (reps<2){reps=2} else{reps=reps}
    S=round(N*.999)
    l=1

    ##### Bootstrap sampling loop. You may use data <- mydata instead of data <- mydata[sample(1:N, S),  TRUE]
    #if you just want see how it works for the original sample data.
    while (l<reps){
      # IV regression
      set.seed(3*(l))   # a different seed for each sub sample
      # data <- data[sample(1:N, S),  TRUE]
      data_sample <- data[sample(1:N, S, replace = TRUE), ]
      ####Computation of m1

      dd<-3#end value for delta
      d<-0.01# starting value for delta
      delt<-0.01# step to change delta
      i<-1 ### starting value for a counter
      ## placeholders for variables
      m1=0
      st=0
      ev22=0
      dv=0
      dv2=0
      x4=0
      l1=0
      l2=0
      while (d<dd){# we compute m1 until siv is close to become perpendicular to x
        data_sample$siv<-(data_sample$x-k*d*data_sample$R)
        #rls<-(lm(p401k~siv+inc+incsq+age+agesq+marr+fsize, data=data))
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
        ssr1 <- sum((predict(lm((ev21^2)~siv,data=data_sample))-mean(data_sample$ev21^2))^2)
        sse1=sum(data_sample$ev21^2)
        x1= (ssr1/2)/(sse1/n^2)^2
        ssr2 <- sum((predict(lm((ev22^2)~siv,data=data_sample))-mean(data_sample$ev22^2))^2)
        sse2=sum(data_sample$ev22^2)
        x2= (ssr2/2)/(sse2/n^2)^2
        #dv[i] <-pchisq(x2, df =1,lower.tail=FALSE)-pchisq(x1, df =1,lower.tail=FALSE)#
        x3 <- x1/x2#sumsq(predict(lm((ev22^2)~siv,data=data))/predict(lm((ev21^2)~siv,data=data))/predict(lm((ev21^2)~siv,data=data)))
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
            # #### DT condition of homoscedatic case
      d0 <- (which.min(abs(m1)))*delt
      d0i[l] <- d0

      ### DT point for heteroscedastic case- parametric approach
      d0r <-  which.min(dv2)*delt
      d0ri[l] <- d0r

      #### DT point for heteroscedastic case- non-parametric approach
      d0rn <- which.min(x4)*delt
      d0rni[l] <- d0rn

        l <- l+1
    }

   ### random shocks for SIVs
  v1 <- rnorm(N, 0, sd(x))
  v2 <- rnorm(N, 0, sd(x))
  v3 <- rnorm(N, 0, sd(x))
    ### final stage estimations
    ###### Simple homogenous assumption case
    d0i <-  d0i[complete.cases(d0i)]
    d0m <- mean(d0i)
    siv_list[[paste0("siv1", g)]]<-(data$x-k*d0m*data$R)
    siv_list[[paste0("siv1b",g)]]<-(data$x-k*d0m*data$R+v1)
    ############Paramteric heterogenous case
    d0ri <-  d0ri[complete.cases(d0ri)]
    d0rm <- mean(d0ri)
    siv_list[[paste0("siv2", g)]]<-(data$x-k*d0rm*data$R)
    siv_list[[paste0("siv2b",g)]]<-(data$x-k*d0rm*data$R+v2)
      #### Non-parameteric heterogenous case
    d0rni <-  d0rni[complete.cases(d0rni)]
    d0rnm <- mean(d0rni)
    siv_list[[paste0("siv3", g)]]<-(data$x-k*d0rnm*data$R)
    siv_list[[paste0("siv3b",g)]]<-(data$x-k*d0rnm*data$R+v3)

  }else{
    v=rnorm(N,0,1)
    print("NO endogeneity problem. All SIV estimates are the same as the OLS")
    d0m=0.001
    siv_list[[paste0("siv1", g)]]<-(data$x-k*d0m*data$R)
      ############Parametric heterogeneous case
    d0rm <- 0.001
    siv_list[[paste0("siv2", g)]]<-(data$x-k*d0rm*data$R)
    #### Non-parameteric heterogenous case
    d0rnm <- 0.001
    siv_list[[paste0("siv3", g)]]<-(data$x-k*d0rnm*data$R)

  }


}
  my_df <- data.frame(matrix(nrow = N, ncol = 0))  # empty data frame with 5 rows
  for (name in names(siv_list)) {
  my_df[[name]] <- siv_list[[name]]
}

  for (name in names(siv_list)) {
    data[[name]] <- siv_list[[name]]
  }

  vars1 <- names(my_df)[grepl("^siv1", names(my_df))]
  vars2 <- names(my_df)[grepl("^siv2", names(my_df))]
  vars3 <- names(my_df)[grepl("^siv3", names(my_df))]
  ###### Simple homogeneous assumption case
  formula_str <- paste(paste0(Y," ~ "), paste(paste(E,collapse ="+"), "+"), paste(H, collapse = " + "))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object

  iv_str <- paste(paste0(" ~ "), paste(vars1, collapse = " + "), paste0("+"),paste( H, collapse = " + "))## Construct formula as a string
  instruments <- as.formula(iv_str)# Convert to formula object
  iv1<-ivreg(formula, instruments, data=data)
   ############Parametric heterogeneous case
  iv_str <- paste(paste0(" ~ "), paste(vars2, collapse = " + "), paste0("+"),paste( H, collapse = " + "))## Construct formula as a string
  instruments <- as.formula(iv_str)# Convert to formula object
  iv2<-ivreg(formula, instruments, data=data)

  #### Non-parametric heterogeneous case
  v_str <- paste(paste0(" ~ "), paste(vars3, collapse = " + "), paste0("+"),paste( H, collapse = " + "))## Construct formula as a string
  instruments <- as.formula(iv_str)# Convert to formula object
  iv3<-ivreg(formula, instruments, data=data)

  return(list(IV1 = (iv1), IV2 = (iv2), IV3 = (iv3),siv_list = (siv_list), signk = (signk)))

              }
