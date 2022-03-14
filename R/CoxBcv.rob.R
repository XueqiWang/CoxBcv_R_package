#' Uncorrected robust sandwich variance estimator
#' 
#' Calculate the uncorrected robust sandwich variance estimator for marginal Cox analysis of cluster randomized trials (Spiekerman and Lin, 1998).
#' 
#' @param Y vector of observed time-to-event data.
#' @param Delta vector of censoring indicators.
#' @param X matrix of marginal mean covariates with one column for one covariate (design matrix excluding intercept).
#' @param ID vector of cluster identifiers.
#' 
#' @return
#' \itemize{
#'   \item coef - estimate of coefficients.
#'   \item exp(coef) - estimate of hazard ratio.
#'   \item ROB-var - uncorrected robust sandwich variance estimate of coef.
#' }
#' 
#' @export
#' 
#' @references 
#' Spiekerman, C. F., & Lin, D. Y. (1998). 
#' Marginal regression models for multivariate failure time data. 
#' Journal of the American Statistical Association, 93(443), 1164-1175.
#' 
#' @examples 
#' Y <- c(11,19,43,100,7,100,100,62,52,1,7,6)
#' Delta <- c(1,1,1,0,1,0,0,1,1,1,1,1)
#' X1 <- c(0,0,0,0,0,0,1,1,1,1,1,1)
#' X2 <- c(-19,6,-25,48,10,-25,15,22,17,-9,45,12)
#' ID <- c(1,1,2,2,3,3,4,4,5,5,6,6)
#' 
#' X <- X1
#' CoxBcv.rob(Y,Delta,X,ID)
#' 
#' X <- cbind(X1,X2)
#' CoxBcv.rob(Y,Delta,X,ID)
#' 
#' @importFrom stats coef
#' @import pracma
#' @import survival

CoxBcv.rob <- function(Y,Delta,X,ID){

  ######################################
  # Step 1: prepare data elements
  ######################################
  # point estimate
  test.cox_cluster <- coxph(Surv(Y,Delta)~X+cluster(ID))
  beta <- as.matrix(coef(test.cox_cluster))
  
  # sort observations by time
  b <- order(Y)
  Y <- sort(Y)
  X <- as.matrix(X)[b,,drop=FALSE]
  ID <- ID[b]
  Delta <- Delta[b]
  ny <- length(Y)
  nbeta <- dim(as.matrix(X))[2]
  UID <- sort(unique(ID))
  n <- length(UID)
  
  IDind <- zeros(length(UID), ny)
  for (i in 1:length(UID)){
    IDind[i, ID==UID[i]] <- 1
  }
  
  ######################################################
  # Step 2: Model-based variance estimates
  ######################################################
  # the rate of counting process of event, or dN(t)
  NN <- diag(Delta)
  # use the following trick to obtain the at-risk process Y(t)
  # each row is an individual
  # each column is a specific time point (recall the counting process notation)
  IndYY <- (t(repmat(t(Y),ny,1))>=repmat(t(Y),ny,1))
  Xbeta <- c(X%*%beta)
  
  # Three S matrices for variance calculation
  S0beta <- colSums(IndYY*exp(Xbeta))
  S1beta <- zeros(nbeta,ny)
  for(k in 1:nbeta){
    S1beta[k,] <- colSums(IndYY*exp(Xbeta)*X[,k])
  }
  S2beta <- array(0, c(ny,1,nbeta,nbeta))
  for(k in 1:nbeta){
    for(s in 1:nbeta){
      S2beta[,,k,s] <- colSums(IndYY*exp(Xbeta)*X[,k]*X[,s])
    }
  }
  
  Omega <- array(0,c(nbeta,nbeta,n))
  
  # obtain cluster-specific matrices
  for(i in UID){
    # subset observations from each cluster
    S0beta_c <- S0beta[IDind[i,]==1]
    S1beta_c <- S1beta[,IDind[i,]==1,drop=FALSE]
    Delta_c <- Delta[IDind[i,]==1]
    
    # components for A matrix
    for (k in 1:nbeta){
      for (s in 1:nbeta){
        Omega[k,s,i] <- sum(Delta_c*(S2beta[IDind[i,]==1,,k,s]/S0beta_c-
                                       S1beta_c[k,]*S1beta_c[s,]/S0beta_c^2))
      }
    }
  }
  
  Ustar <- apply(Omega,c(1,2),sum)
  naive <- solve(Ustar)
  
  ######################################################
  # Step 3: Robust variance estimate
  ######################################################
  # Breslow estimator of baseline hazard
  dHY <- colSums(NN)/c(S0beta)
  HY <- cumsum(dHY)
  
  # obtain martingale increment: 
  # recall that the martingale is the "residual" in survival context
  epsilon <- NN-IndYY*repmat(t(dHY),ny,1)*exp(Xbeta)
  nom <- zeros(nbeta,n)
  
  # obtain cluster-specific matrices
  for(i in UID){
    # subset observations from each cluster
    X_c <- X[IDind[i,]==1,,drop=FALSE]
    epsilon_c <- epsilon[IDind[i,]==1,,drop=FALSE]
    epsilon_c_all <- colSums(epsilon_c)
    S0beta_c <- S0beta[IDind[i,]==1]
    S1beta_c <- S1beta[,IDind[i,]==1,drop=FALSE]
    ny_c <- sum(IDind[i,]==1)
    Delta_c <- Delta[IDind[i,]==1]
    IndYY_c <- IndYY[IDind[i,]==1,,drop=FALSE]
    Xbeta_c <- Xbeta[IDind[i,]==1]
    ylxb_c <- IndYY_c*exp(Xbeta_c)*repmat(dHY,ny_c,1)
    
    # the trick is to loop through the dimension of the coefficients
    # otherwise need to deal with multi-dimensional array, very complex
    for (k in 1:nbeta){
      # components for B matrix
      tempk <- repmat(X_c[,k,drop=FALSE],1,ny)-repmat(S1beta[k,,drop=FALSE]/t(S0beta),ny_c,1)
      nom[k,i] <- sum(as.matrix(tempk*epsilon_c)%*%repmat(1,ny,1))
    }
  }
  
  # robust variance estimator
  UUtran <- tcrossprod(nom)
  robust <- naive%*%UUtran%*%naive
  
  #############################################
  # Output
  # robust: robust sandwich var
  #############################################
  bvarROB <- diag(robust)
  outbeta <- cbind(matrix(summary(test.cox_cluster)$coefficients[,1:2],ncol=2),bvarROB)
  colnames(outbeta) <- c("coef","exp(coef)","ROB-var")
  
  return(list(outbeta=outbeta))
}


