library(ggplot2)
library(ggpubr)
library(readxl)
dat <- as.matrix(read_excel("math_test.xlsx"))

# Check data distribution
plot(density(dat), main = "Density Plot of Your Variable", xlab = "Variable Values", ylab = "Density")


n <- nrow(dat)
p <- ncol(dat)

##### Generate ability values referring to the start.val.rasch function in the Itm package
rs <- as.vector(rowSums(dat, na.rm = TRUE)) # Calculate the sum of each row in the dataset dat, convert the result to a one-dimensional vector and assign it to rs.
len.uni <- length(unique(rs)) # Calculate the number of unique values in rs, and assign the result to len.uni.
rs <- factor(rs, labels = 1:len.uni) # Convert rs into a factor variable and set the labels to range from 1 to len.uni.
rs <- as.numeric(levels(rs))[as.integer(rs)] # Convert the factor variable rs into a numeric vector. (rs can be seen as scores)
# Calculate Z-score normalization
z_score <- (rs - mean(rs)) / sd(rs)
z <- cbind(1, z_score)

##### Find maximum likelihood
getMLE <- function(x, y, w) {
  d <- ncol(x)  # Number of columns in x
  beta <- rep(0, d)
  loop  <- 1
  Loop  <- 100  # Set the maximum number of iterations to 100
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(1 - 1 / (1 + exp(x %*% beta))) # Probability vector pr
    H <- t(x) %*% (pr * (1 - pr) * w * x) # Solve (X^T * X)
    S <- colSums((y - pr) * w * x) # Solve X^T*Y
    tryCatch(
      {shs <- NA
      shs <- solve(H, S) }, # Solve beta = (X^T * X)^(-1)X^T*Y
      error=function(e){
        cat("\n ERROR :", loop, conditionMessage(e), "\n")}) # Used to catch errors when solving linear equations
    if (is.na(shs[1])) {
      msg <- "Not converge"
      beta <- loop <- NA
      break
    }
    beta.new <- beta + shs
    tlr  <- sum((beta.new - beta)^2)
    beta  <- beta.new
    if(tlr < 0.000001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop)
      warning("Maximum iteration reached")
    loop  <- loop + 1
  }
  list(par= beta/beta[length(beta)], message=msg, iter=loop)
  # The estimated difficulty parameter remains 1 by default, consistent with model assumptions
}

##### Two-step method
twostep <- function(X, Y, r1, r2, method=c("mvc", "mmse", "uni")) {
  call <- match.call()
  method <- match.arg(method)
  n <- length(Y)
  if (method == "uni") {
    idx.simp <- sample(1:n, r1+r2 , T)   # Select subsamples
    x.simp <- X[idx.simp,]
    y.simp <- Y[idx.simp]
    pinv.simp <- n
    fit.simp <- getMLE(x=x.simp, y=y.simp, w=pinv.simp)   # Initially obtain maximum likelihood estimation with weight as n
    beta.simp <- fit.simp$par
    msg <- fit.simp$message
    
      p.simp  <- 1 - 1 / (1 + exp(c(x.simp %*% beta.simp)))   # Directly approximate the total sample after convergence, find results
      w.simp <- p.simp * (1 - p.simp)
      W.simp <- solve(t(x.simp) %*% (x.simp * (w.simp * n))) * (r1+r2) * n
      Vc.simp <- t(x.simp) %*% (x.simp * (y.simp-p.simp)^2 * n^2) / (r1+r2)^2 / n^2
      V.simp <- W.simp %*% Vc.simp %*% W.simp
      amse.simp <- sqrt(diag(V.simp))
    return(list(par=beta.simp, amse=amse.simp, msg=msg, method="uni")) # amse is the asymptotic mean square error of the estimator
  }else{

  # In step 1, different subsampling probabilities can be specified for sets S0 and S1,
  # Each subsampling probability equals half the inverse of the set size. Its purpose is to balance the number of 0s and 1s in subsample responses.
  # If the full data is highly unbalanced, the MLE obtained from this method's subsamples has a higher chance of existing than those obtained through uniform subsampling.
  # This procedure is known as case-control sampling. If the proportion of 1s is close to 0.5, uniform SSP is preferred in step 1 due to its simplicity.
  
  n1 <- sum(Y)
  n0 <- n - n1
  PI.prop <- rep(1/(2*n0), n)
  PI.prop[Y==1] <- 1/(2*n1)
  idx.prop <- sample(1:n, r1, T, PI.prop)   # Case-control sampling selects r1 subsamples, balancing the numbers of 0s and 1s, beneficial for calculating MLE
  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop]
  pinv.prop <- n
  pinv.prop <- 1/PI.prop[idx.prop]   # For samples not drawn, the weight is n; for drawn samples, replace the weights with the inverse of their probabilities
  fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
  beta.prop <- fit.prop$par
  if (is.na(beta.prop[1]))
    return(list(opt=NA, msg="first stage not converge"))
  
  
  if (method == "mmse") {
    P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
    p.prop <- P.prop[idx.prop]
    w.prop <- p.prop * (1 - p.prop)
    W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))
    PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE / sum(PI.mMSE)   # Use formula (10) to find sampling probabilities
    idx.mMSE <- sample(1:n, r2, T, PI.mMSE)
    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]   # Combine r1 and r2 samples
    pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)    # Similarly, change the sampling probability weights for r2 samples on the original basis, replacing them with inverses
    fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)
    
    ru <- length(y.mMSE)
    beta.mMSE <- fit.mMSE$par
    p.mMSE  <- 1 - 1 / (1 + exp(c(x.mMSE %*% beta.mMSE)))
    w.mMSE <- p.mMSE * (1 - p.mMSE)
    W.mMSE <- solve(t(x.mMSE) %*% (x.mMSE * (w.mMSE * pinv.mMSE))) * ru * n
    Vc.mMSE <- t(x.mMSE) %*% (x.mMSE * (y.mMSE-p.mMSE)^2 * pinv.mMSE^2) / ru^2 / n^2
    V.mMSE <- W.mMSE %*% Vc.mMSE %*% W.mMSE
    
    amse <- sqrt(diag(V.mMSE))
    msg <- c(fit.prop$message, fit.mMSE$message)
    return(list(par=beta.mMSE, amse=amse, msg=msg, method="mmse",beta.prop=beta.prop))
  }
  
  if (method == "mvc") {
    P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
    PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
    PI.mVc <- PI.mVc / sum(PI.mVc)
    idx.mVc <- sample(1:n, r2, T, PI.mVc)
    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
    fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
    
    ru <- length(y.mVc)
    beta.mVc <- fit.mVc$par
    p.mVc  <- 1 - 1 / (1 + exp(c(x.mVc %*% beta.mVc)))
    w.mVc <- p.mVc * (1 - p.mVc)
    W.mVc <- solve(t(x.mVc) %*% (x.mVc * (w.mVc * pinv.mVc))) * ru * n
    Vc.mVc <- t(x.mVc) %*% (x.mVc * (y.mVc-p.mVc)^2 * pinv.mVc^2) / ru^2 / n^2
    V.mVc <- W.mVc %*% Vc.mVc %*% W.mVc
    
    amse <- sqrt(diag(V.mVc))
    msg <- c(fit.prop$message, fit.mVc$message)
    return(list(par=beta.mVc, amse=amse, msg=msg, method="mvc",beta.prop=beta.prop))
  }
  }
}


##### Calculate parameter values
getbeta <- function(x,y,r1,r2,p,method=c("mvc", "mmse", "uni")){
  beta <- c() 
  amse <- c() 
  for (i in 1:p) {
    result <- twostep(x, y[, i], r1, r2, method)
    beta <- c(beta, result$par[1])
    amse <- c(amse, result$amse[1]) 
  }
  out <- cbind(beta, amse)
  return(out)
}
# Combine estimates and standard errors into a matrix for output

# Calculate full data parameter estimation
getfullbeta <- function(x, y, p) {
  full.fit <- c()
  for (i in 1:p) {
    result <- getMLE(x, y[, i], 1)
    full.fit <- c(full.fit, result$par[1])
  }
  return(full.fit)
}

# Summarize empirical and estimated MSE
summ <- function(method, n.s, xdist = NA, crate = NA, sigma = NA,
                 real = FALSE) {
  if(real){
    full.coe <- full.fit
    res <- eval(parse(text = paste0("real", ".", method, ".", n.s)))
  }else{
    full <- coe.full[[which(dat.xdist == xdist & dat.crate == crate 
                            & dat.sig == sigma)]]
    full.coe <- full[, 1]
    res <- eval(parse(text = paste0(xdist, crate, ".", sigma,
                                    ".", method, ".", n.s)))
  }
  tmp <- array(unlist(res), dim = c(nrow(res), ncol(res), length(res)))
  emp.mse <- sum(rowMeans((tmp[, 1, ] - full.coe) ^ 2))
  est.mse <- sum(rowMeans(tmp[, 2, ] ^ 2))
  atime <- sum(tmp[1, 3, ])
  if(real){
    return(data.frame(value = c(est.mse, emp.mse, atime), 
                      type = c("Estimated", "Empirical", "Time"),
                      method = rep(method, 3), n.s = rep(n.s, 3)))
  }else{
    return(data.frame(value = c(est.mse, emp.mse, atime), 
                      type = c("est.mse", "emp.mse", "Time"),
                      method = rep(method, 3), xdist = rep(xdist, 3),
                      crate = rep(crate, 3), sigma = rep(sigma, 3),
                      n.s = rep(n.s, 3)))
  }
}


##### Calculate estimates under different samplings
library(tcltk)
pb <- tkProgressBar("Progress","Completed %", 0, 100) 
# Using parallel processing
library(parallel)

RNGkind("L'Ecuyer-CMRG")
set.seed(999)
seeds <- list(.Random.seed)
for (i in 2:1000) {
  seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
}
# Create a list containing 1000 different random number seeds.

for (i in seq(200, 1000, 200)) {
  for (k in c("mmse", "mvc", "uni")) {
    
    # Start cluster
    clust <- makeCluster(10)
    
    # Export variables to each cluster
    clusterExport(clust, c("seeds", "dat", "z","p","i", "k","getMLE","twostep","getbeta"))
    
    eval(parse(text = paste0("real", ".", k, ".", i,'<- 
    
     parLapply(clust, seeds, function(seed, dat, n.s, meds){
        assign(".Random.seed", seed, envir = .GlobalEnv)
        t <- system.time(res <- getbeta(z, dat, 200, n.s, p , meds))
        return(cbind(res, "Time" = c(t[1], rep(NA, (length(res) - 1)))))},
      dat = dat, meds = k, n.s = i)')))
    
    # Stop cluster
    stopCluster(clust) 
  }
  info <- sprintf("Completed %d%%", round(i/10))
  setTkProgressBar(pb, i/10, sprintf("Progress (%s)", info),info)
}

# Close progress bar
close(pb)

##### Time and parameter estimation for full data
time.full <- system.time({
  full.fit <-  getfullbeta(z, dat, p)
  })[1]

##### Derive empirical standard deviation of maximum likelihood estimation for full data using Bootstrap method

library(parallel)
clust <- makeCluster(10)
clusterExport(clust, c("seeds", "dat", "z","p", "n","getMLE","twostep","getbeta"))

boots <- parLapply(clust, seeds, function(seed, dat){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  getbeta(z, dat, 0, n, p , "uni")},
  dat = dat)

stopCluster(clust)
