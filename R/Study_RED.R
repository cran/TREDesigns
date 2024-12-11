#' Studying Properties of Ternary Residual Effect Designs
#'
#' @param design Provide a ternary residual effect design
#'@description
#'To study the properties of any given ternary residual effect design.
#'
#' @return It returns Information matrix (C), Average Variance Factor (AVF), and Canonical Efficiency Factor (CEF) for both treatment and residual effects for a given ternary residual design.
#' @export
#'
#' @examples
#' library(TREDesigns)
#' design=PBtRED3(v = 5)$PBTRED
#' Study_tRED(design)
Study_tRED<-function(design){
  a=as.matrix(design)
  session<-nrow(a)
  #####ginv
  # ginv<-function(matrix){
  #   N<-as.matrix(matrix)
  #   NtN<-t(N)%*%N
  #   NNt<-N%*%t(N)
  #   eig_val<-eigen(NtN)$values
  #   positive<-NULL
  #   for(i in eig_val){
  #     if(i>10^-7){
  #       positive<-c(positive,TRUE)
  #     }else{
  #       positive<-c(positive,FALSE)
  #     }
  #   }
  #   sigma<-1/sqrt(eig_val[eig_val>10^-7])
  #   u<-eigen(NtN)$vectors[,positive,drop=FALSE]
  #   v<-eigen(NNt)$vectors[,positive,drop=FALSE]
  #   return(u%*%diag(sigma,nrow(t(v)))%*%t(v))
  # }
  # Step 1: Initialize matrices
  m <- matrix(1, nrow = nrow(a) * ncol(a), ncol = 1)  # Mean vector
  c <- matrix(0, nrow = nrow(a), ncol = ncol(a))

  # Step 2: Fill c matrix based on sessions
  k <- 2
  kk <- 1
  for (i in seq_along(session)) {
    for (j in 1:(session[i] - 1)) {
      c[k, ] <- a[kk, ]
      k <- k + 1
      kk <- kk + 1
    }
    k <- k + 1
    kk <- kk + 1
  }

  # Step 3: Construct trt matrix (design matrix for direct treatment)
  trt <- matrix(0, nrow = nrow(a) * ncol(a), ncol = max(a))
  k <- 1
  for (i in 1:nrow(a)) {
    for (j in 1:ncol(a)) {
      if (a[i, j] > 0) {
        trt[k, a[i, j]] <- 1
        k <- k + 1
      }
    }
  }

  # Step 4: Construct residual matrix (design matrix for residual treatment)
  residual <- matrix(0, nrow = nrow(c) * ncol(c), ncol = max(a))
  k <- 1
  for (i in 1:nrow(c)) {
    for (j in 1:ncol(c)) {
      if (c[i, j] > 0) {
        residual[k, c[i, j]] <- 1
      }
      k <- k + 1
    }
  }

  # Step 5: Construct periods matrix
  periods <- matrix(0, nrow = nrow(a) * ncol(a), ncol = nrow(a))
  k <- 1
  for (i in 1:nrow(a)) {
    for (j in 1:ncol(a)) {
      if (a[i, j] > 0) {
        periods[k, i] <- 1
        k <- k + 1
      }
    }
  }
  ########panellist
  # Initialize the panelist matrix (design matrix - obs VS column)
  panelist <- matrix(0, nrow = nrow(a) * ncol(a), ncol = ncol(a))

  # Fill the panelist matrix based on observations and columns
  k <- 1
  for (i in 1:nrow(a)) {
    for (j in 1:ncol(a)) {
      if (a[i, j] > 0) {
        panelist[k, j] <- 1
      }
      k <- k + 1
    }
  }

  # Print the panelist matrix to verify
  #print(panelist)
  ###session
  # Initialize the sessions matrix (design matrix - obs VS session)
  sessions <- matrix(0, nrow = nrow(a) * ncol(a), ncol = length(session))

  # Variables to track rows in `a` and observations in `sessions`
  kk <- 1
  l <- 1

  # Fill the sessions matrix based on `session` and `a`
  for (k in 1:length(session)) {
    for (i in 1:session[k]) {
      for (j in 1:ncol(a)) {
        if (a[l, j] > 0) {
          sessions[kk, k] <- 1
        }
        kk <- kk + 1
      }
      l <- l + 1
    }
  }

  # Print the sessions matrix to verify
  #print(sessions)

  # Step 6: Compute joint C matrix
  x1 <- cbind(trt, residual)
  x2 <- cbind(m, periods,panelist,sessions)
  c_mat_joint <- t(x1)%*%x1 - t(x1)%*%x2 %*% MASS::ginv(t(x2)%*%x2) %*% (t(x2)%*%x1)

  c11 <- c_mat_joint[1:max(a), 1:max(a)]
  c22 <- c_mat_joint[(max(a) + 1):(2 * max(a)), (max(a) + 1):(2 * max(a))]
  c12 <- c_mat_joint[1:max(a), (max(a) + 1):(2 * max(a))]

  c_trt <- c11 - c12 %*% MASS::ginv(c22) %*% t(c12)
  c_resid <- c22 - t(c12) %*% MASS::ginv(c11) %*% c12
  ##AVF

  # Define the maximum value of `a`
  max_a <- max(a)

  # Compute the total number of combinations
  t1 <- choose(max_a, 2)

  # Initialize the matrix `cot` with zeros
  cot <- matrix(0, nrow = t1, ncol = max_a)

  # Fill the `cot` matrix
  k <- 1
  for (i in 1:(max_a - 1)) {
    for (j in (i + 1):max_a) {
      cot[k, i] <- 1
      cot[k, j] <- -1
      k <- k + 1
    }
  }

  # Compute covt
  covt <- cot %*% MASS::ginv(c_trt) %*% t(cot)
  vart1 <- diag(covt)
  onet <- matrix(1, nrow = t1, ncol = 1)
  variance <- vart1 %*% onet

  # Compute covr
  covr <- cot %*% MASS::ginv(c_resid) %*% t(cot)
  vart1r <- diag(covr)
  variance_r <- vart1r %*% onet

  # Calculate average variances
  av_var <- sum(variance) / length(variance)
  av_var_r <- sum(variance_r) / length(variance_r)

  # Step 7: Eigenvalues and Canonical Efficiency Factor
  eig_trt <- (eigen(c_trt)$values)
  eig_resid <-(eigen(c_resid)$values)

  # Filter positive eigenvalues
  eig_trt <- eig_trt[(eig_trt) > (10^(-7))]
  eig_resid <- eig_resid[(eig_resid) > (10^(-7))]

  # Compute canonical efficiency factors
  rep_trt <-(t(trt)%*%(trt))[1, 1]
  rep_resid <- (t(residual)%*%residual)[1, 1]

  can_eff_factor_trt <- (length(eig_trt) / sum(1 / (eig_trt / rep_trt)))
  can_eff_factor_resid <- (length(eig_resid) / sum(1 / (eig_resid / rep_resid)))

  # Output results
  lm<-list("Cmatrix_trt"=round(c_trt,3),"Cmatrix_residual"=round(c_resid,3),"AVF_trt"=round(av_var,3),"AVF_residual"=round(av_var_r,3),"CEF_trt"=round(can_eff_factor_trt,3),"CEF_residual"=round(can_eff_factor_resid,3))
  return(lm)
}

