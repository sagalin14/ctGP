library(AdequacyModel)
library(globpso)

library(RcppArmadillo)


# Create Armadillo function
Rcpp::cppFunction(depends = "RcppArmadillo", code = '
                  Rcpp::NumericMatrix solveCpp(SEXP A) {
                  arma::mat A_c = (arma::mat)Rcpp::as<arma::mat>(A);
                  arma::mat B_c = inv_sympd(A_c);
                  return wrap(B_c);
                  }')

QQGP.ori.model.QQ5level <- function(design, data, initParam, upperBound = NULL, lowerBound = NULL){
  
  #cat("QQGP.ori.model is fitted. \n" )
  
  # transform the data and design to matrix to prevent troubles
  design <- as.matrix(design)
  #data <- as.matrix(data)
  
  # define some dimensions about data
  contiDim <- ifelse(is.null(dim(design)), 1, dim(design)[2])
  levelVec <- unique(data[, contiDim+1])
  levelVecTotal <- data[, contiDim+1]
  levelNum <- length(levelVec)
  contiPosNum <- ifelse(is.null(dim(data)), length(design), dim(data)[1])
  designTmp <- data[, 1:contiDim, drop=FALSE]
  totalPtsNum <- contiPosNum
  
  
  # seperate resonse
  yVec <- data[, contiDim+2]
  
  
  # matrix index
  rcoord <- cbind(rep(seq_len(totalPtsNum-1L), times = rev(seq_len(totalPtsNum-1L))), 
                  unlist(lapply(X=rev(seq_len(totalPtsNum-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = totalPtsNum)))
  # design points differences
  xdiff <- (designTmp[rcoord[, 2L], ] - designTmp[rcoord[, 1L], ])**2
  ###############################################################################
  ####                                                                       ####
  ####                                                                       ####
  ####             STEP 1 : objective function and optimization              ####
  ####                      with respect to the parameters                   ####
  ####                                                                       ####
  ####                                                                       ####
  ###############################################################################
  
  # ========================== Logliklihood Function ========================== #
  
  QQGPobjFun <- function(par, x, design, data){
    
    # transform the data and design to matrix to prevent troubles
    designTmp <- as.matrix(design)
    #data <- as.matrix(data)
    
    # define some dimensions about data
    contiDim <- ifelse(is.null(dim(designTmp)), 1, dim(designTmp)[2])
    contiPosNum <- ifelse(is.null(dim(data)), length(design), dim(data)[1])
    design <- data[, 1:contiDim, drop=FALSE]
    levelVec <- unique(data[, contiDim+1])
    levelVecTotal <- data[, contiDim+1]
    levelNum <- length(levelVec)
    totalPtsNum <- contiPosNum
    
    # seperate resonse
    yVec <- data[, contiDim+2]
    
    # seperate parameters vector
    fi <- par[1:contiDim]
    theta <- par[(contiDim+1):(contiDim+choose(levelNum, 2))]
    
    # correlation matrix of continuous positions (Hmatrix)
    Hmatrix <- matrix(0, ncol = totalPtsNum, nrow = totalPtsNum)
    # matrix index
    rcoord <- cbind(rep(seq_len(totalPtsNum-1L), times = rev(seq_len(totalPtsNum-1L))), 
                    unlist(lapply(X=rev(seq_len(totalPtsNum-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = totalPtsNum)))
    # design points differences
    xdiff <- (design[rcoord[, 2L], ] - design[rcoord[, 1L], ])**2
    Beta <- matrix(fi, nrow = choose(totalPtsNum, 2), ncol = contiDim, byrow = TRUE)
    Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
    Hmatrix[rcoord] <- Htemp; Hmatrix <- Hmatrix + t(Hmatrix); Hmatrix <- exp(-Hmatrix)
    
    # to prevent computational sigular problem, if the condition number of Hmatrix is too small (<= 1e-6)
    # we add a small nugget value (1e-10) to the diagonal line of Hmatrix
    #if(rcond(Hmatrix) <= 1e-6){
    #  Hmatrix <- Hmatrix + diag(rep(1e-10, totalPtsNum))
    #}
    
    #-------------------------------------------------------#
    
    # correlation matrix of categorical combinations (Tmatrix) (by using sphere decomposition)
    Lmatrix <- diag(rep(0, levelNum)); Lmatrix[1, 1] <- 1
    sphereDecLowRowAndLen <- c(2:levelNum)
    LowerTemp <- diag(rep(0, levelNum)); LowerTemp[lower.tri(LowerTemp)] <- theta#; LowerTemp <- t(LowerTemp)
    for (r in sphereDecLowRowAndLen){
      Lmatrix[r,1] = cos(LowerTemp[r,1])
      for(s in 2:r){
        if(s<r){
          Lmatrix[r,s] <- prod(sapply(LowerTemp[r, 1:(s-1)], FUN = sin))*cos(LowerTemp[r, s])
        }else{
          Lmatrix[r,s] <- prod(sapply(LowerTemp[r, 1:(s-1)], FUN = sin))
        }
      }
    }
    Tmatrix <- Lmatrix %*% t(Lmatrix)
    #if(rcond(Tmatrix) <= 1e-6){
    #  Tmatrix <- Tmatrix + diag(rep(1e-10, levelNum))
    #}
    
    #-------------------------------------------------------#
    
    # correlation matrix of whole data (Rmatrix) (by using kronecker product)
    # CAREFULLY using kronecker product function, since the order of input should be inverse
    # Rmatrix <- kronecker(Hmatrix, Tmatrix, FUN = "*", make.dimnames = FALSE)
    Rmatrix <- matrix(0, ncol = totalPtsNum, nrow = totalPtsNum)
    for(row in 1:totalPtsNum){
      for(col in 1:totalPtsNum){
        Rmatrix[row, col] <- Hmatrix[row, col] * Tmatrix[levelVecTotal[row], levelVecTotal[col]]
      }
    }
    
    # inverse of H.matrix T.matrix and R.matrix (should not have singular problem now)
    # TinvMatrix <- solve(Tmatrix)
    # HinvMatrix <- solve(Hmatrix)
    # RinvMatrix <- kronecker(HinvMatrix, TinvMatrix, FUN = "*", make.dimnames = FALSE)
    if(rcond(Rmatrix) <= 1e-6){
      Rmatrix <- Rmatrix + diag(rep(1e-10, totalPtsNum))
    }
    RinvMatrix <- solveCpp(Rmatrix)
    # RmatrixLogDet <- contiPosNum*log(det(Tmatrix)) + levelNum*log(det(Hmatrix))
    
    #-------------------------------------------------------#
    
    # estimators
    FVec <- matrix(0, ncol = levelNum, nrow = totalPtsNum)
    for(i in 1:totalPtsNum){
      FVec[i,which(data[i,(contiDim+1)] == (1:levelNum))] <- 1
    }
    FVec[,1] <- rep(1, totalPtsNum)
    
    toBeInvMat <- t(FVec)%*%RinvMatrix%*%FVec
    for(row in 1:dim(toBeInvMat)[1]){
      for(col in row:dim(toBeInvMat)[2]){
        toBeInvMat[row, col] <- toBeInvMat[col, row]
      }
    }
    
    betaVecHat <- solveCpp(toBeInvMat)%*%t(FVec)%*%RinvMatrix%*%yVec
    sigmaSqHat <- (t(yVec-FVec%*%betaVecHat)%*%RinvMatrix%*%(yVec-FVec%*%betaVecHat)) / totalPtsNum
    
    #-------------------------------------------------------#
    
    # return objective
    # subject <- totalPtsNum*log(sigmaSqHat) + RmatrixLogDet
    subject <- totalPtsNum*log(sigmaSqHat) + log(det(Rmatrix))
    
    return(subject)
  }
  
  
  # make official negative log-likelihood function by pre-setting data and design matrix
  # to simplify the input object to only "par" and "x", thus we can use PSO to solve the optimization process
  QQGPnegLogLik <- function(par, x){
    return(QQGPobjFun(par, x, data = data, design = design))
  }
  
  QQGPnegLogLik2 <- function(x){
    return(QQGPobjFun(par = x, x=NULL, data = data, design = design))
  }
  
  
  # ============================== optimization =============================== #
  # optimize the negLogLik function to find the parameters
  # set.seed(10000)
  # QQGPoptResult <- pso(func = QQGPnegLogLik, S = 32, lim_inf = c(rep(-4, contiDim), rep(0, choose(levelNum, 2))),
  #                      lim_sup = c(rep(4, contiDim), rep(pi, choose(levelNum, 2))), e = 0.05, N=5)
  # The search domain is [-5, 5]^3
  if(is.null(upperBound)){
    if(is.null(dim(xdiff))){
      upp_bound <- c(rep(log10(-log(0.00001)/max(xdiff)), contiDim), rep(1, choose(levelNum, 2)))
      #upp_bound <- c(rep(3, contiDim), rep(3, choose(levelNum, 2)))
    }else{
      upp_bound <- c(rep(log10(-log(0.00001)/max(apply(xdiff, 1,sum))), contiDim), rep(1, choose(levelNum, 2)))
      #upp_bound <- c(rep(3, contiDim), rep(3, choose(levelNum, 2)))
    }
  }else{
    upp_bound <- upperBound
  }
  
  if(is.null(lowerBound)){
    if(is.null(dim(xdiff))){
      low_bound <- c(rep(log10(-log(0.99999)/max(xdiff)), contiDim), rep(-1, choose(levelNum, 2)))
      #low_bound <- c(rep(1, contiDim), rep(-3, choose(levelNum, 2)))
    }else{
      low_bound <- c(rep(log10(-log(0.99999)/max(apply(xdiff, 1,sum))), contiDim), rep(-1, choose(levelNum, 2))) # c(rep(-2, contiDim), rep(-0.5, choose(levelNum, 2)))
      #low_bound <- c(rep(1, contiDim), rep(-3, choose(levelNum, 2)))
    }
  }else{
    low_bound <- lowerBound
  }
  
  
  ######################
  ######################
  #initParam <- (upp_bound+low_bound)/2
  ######################
  ######################
  
  
  # Use getPSOInfo() to change the PSO options
  alg_setting <- getPSOInfo(nSwarm = 32, maxIter = 10, psoType = "quantum")
  res_c_large <- globpso(objFunc = QQGPnegLogLik2, lower = low_bound, upper = upp_bound, PSO_INFO = alg_setting, seed = 1000000,
                         init = c(initParam, rep(0, choose(levelNum, 2))), verbose = F)#, loc = loc_shift)
  # res_c_large$par
  # res_c_large$val
  
  # record the result parameters
  estimatedParameters <- #c(initParam, res_c_large$par[-c(1:length(initParam))])
    res_c_large$par
  # estimatedParameters <- QQGPoptResult$par
  
  
  HmatrixHat <- matrix(0, ncol = totalPtsNum, nrow = totalPtsNum)
  # set the result parameters above as the parameters we will use later
  fiVecHat <- estimatedParameters[1:contiDim]
  thetaVecHat <- estimatedParameters[(contiDim+1):(contiDim+choose(levelNum, 2))]
  Beta <- matrix(fiVecHat, nrow = choose(totalPtsNum, 2), ncol = contiDim, byrow = TRUE)
  Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
  HmatrixHat[rcoord] <- Htemp; HmatrixHat <- HmatrixHat + t(HmatrixHat); HmatrixHat <- exp(-HmatrixHat)
  
  # make estimation of correlation matrix of categorical combinations (TmatrixHat) (by using sphere decomposition)
  Lmatrix <- diag(rep(0, levelNum)); Lmatrix[1, 1] <- 1
  sphereDecLowRowAndLen <- c(2:levelNum)
  LowerTemp <- diag(rep(0, levelNum)); LowerTemp[lower.tri(LowerTemp)] <- thetaVecHat
  for (r in sphereDecLowRowAndLen){
    Lmatrix[r,1] = cos(LowerTemp[r,1])
    for(s in 2:r){
      if(s<r){
        Lmatrix[r,s] <- prod(sapply(LowerTemp[r, 1:(s-1)], FUN = sin))*cos(LowerTemp[r, s])
      }else{
        Lmatrix[r,s] <- prod(sapply(LowerTemp[r, 1:(s-1)], FUN = sin))
      }
    }
  }
  TmatrixHat <- Lmatrix %*% t(Lmatrix)
  
  # make estimation of inverse correlation matrix of whole data (RmatrixHat) (by using kronecker product)
  # CAREFULLY using kronecker product function, since the order of input should be inverse
  
  #############
  #levels(levelVecTotal) <- 1:length(levels(levelVecTotal))
  #############
  
  RmatrixHat <- matrix(0, ncol = totalPtsNum, nrow = totalPtsNum)
  for(row in 1:totalPtsNum){
    for(col in 1:totalPtsNum){
      #print(row); print(col)
      RmatrixHat[row, col] <- HmatrixHat[row, col] * TmatrixHat[levelVecTotal[row], levelVecTotal[col]]
    }
  }
  
  if(rcond(RmatrixHat) <= 1e-6){
    RmatrixHat <- RmatrixHat + diag(rep(1e-10, totalPtsNum))
  }
  RinvMatrixHat <- solveCpp(RmatrixHat)
  
  
  # estimators
  FVec <- matrix(0, ncol = levelNum, nrow = totalPtsNum)
  for(i in 1:totalPtsNum){
    FVec[i,which(data[i,(contiDim+1)] == (1:levelNum))] <- 1
  }
  FVec[,1] <- rep(1, totalPtsNum)
  
  toBeInvMat <- t(FVec)%*%RinvMatrixHat%*%FVec
  for(row in 1:dim(toBeInvMat)[1]){
    for(col in row:dim(toBeInvMat)[2]){
      toBeInvMat[row, col] <- toBeInvMat[col, row]
    }
  }
  
  betaVecHat <- solveCpp(toBeInvMat)%*%t(FVec)%*%RinvMatrixHat%*%yVec
  sigmaSqHat <- (t(yVec-FVec%*%betaVecHat)%*%RinvMatrixHat%*%(yVec-FVec%*%betaVecHat)) / totalPtsNum
  
  
  
  
  
  #FVec <- rep(1L,totalPtsNum)
  #betaVecHat <- solveCpp(t(FVec)%*%RinvMatrixHat%*%FVec)%*%t(FVec)%*%RinvMatrixHat%*%yVec
  #sigmaSqHat <- (t(yVec-FVec%*%betaVecHat)%*%RinvMatrixHat%*%(yVec-FVec%*%betaVecHat)) / totalPtsNum
  
  # combine all estimations
  estimatedList <- list(fi = fiVecHat, theta = thetaVecHat, RinvMat = RinvMatrixHat, Tmat = TmatrixHat, design = design, data = data,
                        FVec = FVec, betaVec = betaVecHat, sigmaSq = sigmaSqHat, nugget = sqrt(.Machine$double.eps))
  
  return(estimatedList)
}




QQGP.ori.prediction <- function(QQGP.model, newPoint){
  
  ###############################################################################
  ####                                                                       ####
  ####                                                                       ####
  ####                         STEP 2 : estimation                           ####
  ####                                                                       ####
  ####                                                                       ####
  ###############################################################################
  
  # set the result parameters above as the parameters we will use later
  fiVecHat <- QQGP.model$fi
  thetaVecHat <- QQGP.model$theta
  RinvMatrixHat <- QQGP.model$RinvMat
  TmatrixHat <- QQGP.model$Tmat
  
  # transform the data and design to matrix to prevent troubles
  designTmp <- data.matrix(as.matrix(QQGP.model$design))
  #data <- coerceFun.ori(data, design)#data.matrix(data)
  
  # define some dimensions about data
  contiDim <- ifelse(is.null(dim(designTmp)), 1, dim(designTmp)[2])
  contiPosNum <- ifelse(is.null(dim(QQGP.model$data)), length(design), dim(QQGP.model$data)[1])
  design <- data[, 1:contiDim, drop=FALSE]
  levelVec <- unique(QQGP.model$data[, contiDim+1])
  levelVecTotal <- QQGP.model$data[, contiDim+1]
  levelNum <- length(levelVec)
  totalPtsNum <- contiPosNum
  
  betaVecHat <- QQGP.model$betaVec
  FVec <- QQGP.model$FVec
  
  
  # seperate resonse
  yVec <- QQGP.model$data[, contiDim+2]
  
  # estimation of rediduals vector, this will be use to find the prediction value of newPoint. 
  if(is.null(dim(newPoint))){
    newPointNum <- 1
    rVecHat <- rep(0, totalPtsNum)
    
    predictFVec <- rep(0, levelNum)
    predictFVec[1] <- 1
    predictFVec[which(newPoint[(contiDim+1)] == (1:levelNum))] <- 1
    
    
    for (k in (1:totalPtsNum)){
      tau <- TmatrixHat[as.numeric(QQGP.model$data[k, (dim(QQGP.model$data)[2]-1)]), as.numeric(newPoint[length(newPoint)-1])]
      corr_K <- exp(-sum(sapply(c(1:contiDim), FUN = function(x) return((as.numeric(10^fiVecHat[x])%*%(as.numeric(QQGP.model$data[k,x]-newPoint[x])^2))))))
      rVecHat[k] <- tau*corr_K
    }
    
    yPredict <- as.numeric(predictFVec%*%betaVecHat + (((t(rVecHat)) %*% (RinvMatrixHat)) %*% (yVec-FVec%*%betaVecHat)))
    
  }else{
    newPointNum <- dim(newPoint)[1]
    rVecHat <- matrix(0, ncol = totalPtsNum, nrow = newPointNum)
    
    predictFVec <- matrix(0, ncol = levelNum, nrow = newPointNum)
    predictFVec[,1] <- rep(1, newPointNum)
    for(i in 1:newPointNum){
      predictFVec[i, which(newPoint[i,(contiDim+1)] == (1:levelNum))] <- 1
    }
    
    
    
    for (j in 1:newPointNum){
      for (k in (1:totalPtsNum)){
        tau <- TmatrixHat[as.numeric(QQGP.model$data[k, (dim(QQGP.model$data)[2]-1)]), as.numeric(newPoint[j, dim(newPoint)[2]-1])]
        corr_K <- exp(-sum(sapply(c(1:contiDim), FUN = function(x) return((as.numeric(10^fiVecHat[x])%*%(as.numeric(QQGP.model$data[k,x]-newPoint[j, x])^2))))))
        rVecHat[j, k] <- tau*corr_K
      }
    }
    
    yPredict <- as.numeric(predictFVec%*%betaVecHat + ((rVecHat %*% (RinvMatrixHat)) %*% (yVec-FVec%*%betaVecHat)))
    
  }
  
  # return estimated list and prediction function
  return(yPredict)
}


#startTime <- proc.time()
#####################################################################################
# setwd("~/Desktop/RB/TGP")
# source("generateQQGPoriginTestData.R")
# a <- QQGP.ori.model(design = design[1:19,], data = data[1:19,], initParam =  c(0, 0))
# QQGP.ori.prediction(QQGP.model = a, newPoint = data[20,])
# data[20,4]
#####################################################################################
#endTime <- proc.time()
#(endTime - startTime)[3]



