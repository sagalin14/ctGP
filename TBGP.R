library(AdequacyModel)
library(globpso)
library(RcppArmadillo)
library(plgp)
library(foreach)
library(doParallel)
library(mvtnorm)
library(parallel)
library(dplyr)
library(ggplot2)


#########################################
########## some tool functions ##########
#########################################

# Armadillo inverse function
Rcpp::cppFunction(depends = "RcppArmadillo", code = '
                  Rcpp::NumericMatrix solveCpp(SEXP A) {
                  arma::mat A_c = (arma::mat)Rcpp::as<arma::mat>(A);
                  arma::mat B_c = inv_sympd(A_c);
                  return wrap(B_c);
                  }')

# cov matrix to corr matrix
cov2cor <- function (V)
{
  p <- (d <- dim(V))[1L]
  if (!is.numeric(V) || length(d) != 2L || p != d[2L])
    stop("'V' is not a square numeric matrix")
  Is <- sqrt(1/diag(V))
  if (any(!is.finite(Is)))
    warning("diag(.) had 0 or NA entries; non-finite result is doubtful")
  r <- V
  r[] <- Is * V * rep(Is, each = p)
  r[cbind(1L:p, 1L:p)] <- 1
  r
}


# compute nugget function
computeNug <- function(mat)
{
  ev <- eigen(mat)$values
  checkCond <- abs(max(ev)) - 1e8 * abs(min(ev))
  if(checkCond >= 0){
    nug <- checkCond/(1e8 - 1)
  }else{
    nug <- sqrt(.Machine$double.eps)
  }
  return(nug)
}



# coerce function
coerceFun.ori <- function(data, design, includeSubCate = T){
  
  if(includeSubCate){
    # transform the data and design to matrix to prevent troubles
    design <- as.matrix(design)
    
    # define some parameters about data
    contiDim <- ifelse(is.null(dim(design)), 1, dim(design)[2])
    contiPosNum <- ifelse(is.null(dim(design)), length(design), dim(design)[1])
    cateNum <- dim(data)[2] - contiDim - 1
    cateMatrix <- data[, (contiDim+1):(contiDim+cateNum)]
    levelVec <- character(dim(data)[1])
    for(i in 1:cateNum){
      levelVec <- paste0(levelVec, data[,(contiDim+i)])
    }
    levelVec.unique <- unique(levelVec)
    
    for(i in 1:length(levelVec.unique)){
      levelVec[levelVec == levelVec.unique[i]] <- i 
    }
    
    levelVec <- as.numeric(levelVec)
    
    returnData <- cbind(data.frame(data[, 1:contiDim]), c.T = levelVec, Y = as.numeric(data[, (contiDim+cateNum+1)]))
    
    # return data
    return(returnData)  
  }else{
    
    # transform the data and design to matrix to prevent troubles
    design <- as.matrix(design)
    
    # define some parameters about data
    contiDim <- ifelse(is.null(dim(design)), 1, dim(design)[2])
    contiPosNum <- ifelse(is.null(dim(design)), length(design), dim(design)[1])
    cateNum <- dim(data)[2] - contiDim - 1
    cateMatrix <- data[, (contiDim+1):(contiDim+cateNum)]
    levelVec <- character(dim(data)[1])
    for(i in 1:cateNum){
      levelVec <- paste0(levelVec, data[,(contiDim+i)])
    }
    levelVec.unique <- unique(levelVec)
    
    for(i in 1:length(levelVec.unique)){
      levelVec[levelVec == levelVec.unique[i]] <- i 
    }
    
    levelVec <- as.numeric(levelVec)
    returnData <- cbind(data.frame(data[, 1:contiDim]), Y = as.numeric(data[, (contiDim+cateNum+1)]))
    returnCate <- levelVec
    
    # return data
    return(list(returnData = returnData, levelVec = levelVec))  
  }
}






############################################
########## Gaussian Process model ##########
############################################

GP.model <- function(data, initParam = NULL, range = NULL, cpp = F, nugget = 0){
  
  if(!cpp){
    solveCpp = solve
  }
  
  # transform the data to matrix to prevent troubles
  data <- as.matrix(data)
  
  # define some dimensions about data
  contiDim <- dim(data)[2]-1
  ptsNum <- dim(data)[1]
  
  # seperate resonse
  design <- data[, c(1:contiDim), drop=FALSE]
  yVec <- data[, contiDim+1]
  rcoord <- cbind(rep(seq_len(ptsNum-1L), times = rev(seq_len(ptsNum-1L))), 
                  unlist(lapply(X=rev(seq_len(ptsNum-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = ptsNum)))
  # design points differences
  xdiff <- (design[rcoord[, 2L], ] - design[rcoord[, 1L], ])**2
  
  # ========================== Logliklihood Function ========================== #
  
  GPobjFun <- function(par, data){
    
    # transform the data to matrix to prevent troubles
    data <- as.matrix(data)
    
    # define some dimensions about data
    contiDim <- dim(data)[2]-1
    ptsNum <- dim(data)[1]
    
    # seperate resonse
    design <- data[, c(1:contiDim), drop=FALSE]
    yVec <- data[, contiDim+1]
    
    # seperate parameters vector
    fi <- par
    
    # correlation matrix of continuous positions (Hmatrix)
    Hmatrix <- matrix(0, ncol = ptsNum, nrow = ptsNum)
    # matrix index
    rcoord <- cbind(rep(seq_len(ptsNum-1L), times = rev(seq_len(ptsNum-1L))), 
                    unlist(lapply(X=rev(seq_len(ptsNum-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = ptsNum)))
    # design points differences
    xdiff <- (design[rcoord[, 2L], ] - design[rcoord[, 1L], ])**2
    Beta <- matrix(fi, nrow = choose(ptsNum, 2), ncol = contiDim, byrow = TRUE)
    Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
    Hmatrix[rcoord] <- Htemp; Hmatrix <- Hmatrix + t(Hmatrix); Hmatrix <- exp(-Hmatrix)
    
    # to prevent computational sigular problem, if the condition number of Hmatrix is too small (<= 1e-6)
    # we add a small nugget value to the diagonal line of Hmatrix
    if(rcond(Hmatrix) <= 1e-6){
      Hmatrix <- Hmatrix + diag(rep(nugget, ptsNum))
    }
    
    #-------------------------------------------------------#
    
    # inverse of H.matrix T.matrix and R.matrix (should not have singular problem now)
    HinvMatrix <- solveCpp(Hmatrix)
    
    #-------------------------------------------------------#
    
    # estimators
    FVec <- rep(1L, ptsNum)
    betaVecHat <- solveCpp(t(FVec)%*%HinvMatrix%*%FVec)%*%t(FVec)%*%HinvMatrix%*%yVec
    sigmaSqHat <- (t(yVec-FVec%*%betaVecHat)%*%HinvMatrix%*%(yVec-FVec%*%betaVecHat)) / ptsNum
    
    #-------------------------------------------------------#
    
    # return objective
    subject <- ptsNum*log(sigmaSqHat) + log(det(Hmatrix))
    
    return(subject)
  }
  
  
  # make official negative log-likelihood function by pre-setting data and design matrix
  # to simplify the input object to only "par" and "x", thus we can use PSO to solve the optimization process
  
  GPnegLogLik <- function(x){
    return(as.numeric(GPobjFun(par = x, data = data)))
  }
  
  # ============================== optimization =============================== #
  # optimize the negLogLik function to find the parameters
  
  # define the optimize boundary
  if(is.null(range)){
    if(is.null(dim(xdiff))){
      upp_bound <- rep(log10(-log(0.1)/max(xdiff)), contiDim)
      low_bound <- rep(log10(-log(0.9)/max(xdiff)), contiDim)
    }else{
      upp_bound <- rep(log10(-log(0.1)/max(apply(xdiff, 1,sum))), contiDim)
      low_bound <- rep(log10(-log(0.9)/max(apply(xdiff, 1,sum))), contiDim)
    }
  }else{
    upp_bound <- range[(length(range)/2+1):length(range)]
    low_bound <- range[1:(length(range)/2)]
  }
  
  if(is.null(initParam)){
    initParam <- (0.8*upp_bound+0.2*low_bound)
  }
  
  # Use getPSOInfo() to change the PSO options
  alg_setting <- getPSOInfo(nSwarm = 64, maxIter = 10, psoType = "quantum")
  res_c_large <- globpso(objFunc = GPnegLogLik, lower = low_bound, upper = upp_bound,
                         init = initParam, PSO_INFO = alg_setting, seed = 10000,
                         verbose = F)
  
  # record the result parameters
  estimatedParameters <- res_c_large$par
  
  # correlation matrix of continuous positions (Hmatrix)
  Hmatrix <- matrix(0, ncol = ptsNum, nrow = ptsNum)
  Beta <- matrix(estimatedParameters, nrow = choose(ptsNum, 2), ncol = contiDim, byrow = TRUE)
  Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
  Hmatrix[rcoord] <- Htemp; Hmatrix <- Hmatrix + t(Hmatrix); Hmatrix <- exp(-Hmatrix)
  
  # to prevent computational sigular problem, if the condition number of Hmatrix is too small (<= 1e-6)
  # we add a small nugget value to the diagonal line of Hmatrix
  if(rcond(Hmatrix) <= 1e-6){
    Hmatrix <- Hmatrix + diag(nugget, ptsNum)
  }
  
  # make sure Hmatrix will be symmetric
  Hmatrix <- (Hmatrix + t(Hmatrix))/2
  
  #-------------------------------------------------------#
  
  # inverse of H.matrix T.matrix and R.matrix (should not have singular problem now)
  HinvMatrix <- solveCpp(Hmatrix)
  #-------------------------------------------------------#
  
  # estimators
  FVec <- rep(1L, ptsNum)
  betaVecHat <- solveCpp(t(FVec)%*%HinvMatrix%*%FVec)%*%t(FVec)%*%HinvMatrix%*%yVec
  sigmaSqHat <- (t(yVec-FVec%*%betaVecHat)%*%HinvMatrix%*%(yVec-FVec%*%betaVecHat)) / ptsNum
  
  
  # combine all estimations
  estimatedList <- list(fi = estimatedParameters, design = design, data = data, HinvMat = HinvMatrix,
                        FVec = FVec, betaVec = betaVecHat, sigmaSq =  sigmaSqHat, nugget = nugget,
                        negLogLik = as.numeric(GPnegLogLik(estimatedParameters)))
  
  return(estimatedList)
}




##########################################################
########## Gaussian Process prediction function ##########
##########################################################

GP.prediction <- function(GP.model, newPoint, cpp = F, indep = F){
  
  if(!cpp){
    solveCpp = solve
  }
  
  # transform the data to matrix to prevent troubles
  data <- as.matrix(GP.model$data)
  
  # define some dimensions about data
  contiDim <- dim(data)[2]-1
  ptsNum <- dim(data)[1]
  
  # seperate resonse
  design <- data[, c(1:contiDim), drop=FALSE]
  yVec <- data[, contiDim+1]
  
  # set the result parameters above as the parameters we will use later
  fiVecHat <- GP.model$fi
  
  # make estimation of inverse correlation matrix of whole data (HmatrixHat)
  HinvMatrixHat <- GP.model$HinvMat
  betaVecHat <- GP.model$betaVec
  sigmaSqHat <- GP.model$sigmaSq
  FVec <- GP.model$FVec
  nugget <- GP.model$nugget
  
  if(is.vector(newPoint)){
    newPoint <- matrix(newPoint, nrow = 1)
    np <- newPoint[,1:contiDim, drop=FALSE]
    RX <- distance(t(t(np)*sqrt(as.numeric(10^fiVecHat))), t(t(design)*sqrt(as.numeric(10^fiVecHat))))
    rVecHat <- exp(-RX)
    predictFVec <- rep(1, 1)
    yPredict <- as.numeric(predictFVec%*%betaVecHat + (((rVecHat) %*% (HinvMatrixHat)) %*% (yVec-FVec%*%betaVecHat)))
  }else{
    newPoint <- as.matrix(newPoint)
    np <- newPoint[,1:contiDim, drop=FALSE]
    RX <- distance(t(t(np)*sqrt(as.numeric(10^fiVecHat))), t(t(design)*sqrt(as.numeric(10^fiVecHat))))
    rVecHat <- exp(-RX)
    predictFVec <- rep(1, dim(newPoint)[1])
    yPredict <- as.numeric(predictFVec%*%betaVecHat + (((rVecHat) %*% (HinvMatrixHat)) %*% (yVec-FVec%*%betaVecHat)))
  }
  
  DXX <- distance(t(t(np)*sqrt(as.numeric(10^fiVecHat))))
  SXX <- exp(-DXX) + diag(nugget, dim(DXX))
  predSigmaSq <- (as.numeric(sigmaSqHat) * (SXX - rVecHat %*% HinvMatrixHat %*% t(rVecHat)))
  diag(predSigmaSq) <- pmax(0, diag(predSigmaSq))
  
  if(indep){
    predSigmaSq <- diag(predSigmaSq)
  }
  
  # combine all estimations
  estimatedList <- list(beta = as.numeric(betaVecHat), sigmaSq = as.numeric(sigmaSqHat),
                        predSigmaSq = predSigmaSq)
  
  # return estimated list and prediction function
  return(list(predictValue = yPredict, predSigmaSq = predSigmaSq))
}





################################
########## QQGP model ##########
################################
QQGP.ori.model <- function(design, data, initParam, fiRange = NULL, thetaRange = NULL, Tmatrix = NULL, nugget = 0){
  
  # transform the data and design to matrix to prevent troubles
  design <- as.matrix(design)
  
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
  
  if(is.null(Tmatrix)){
    
    # ========================== Logliklihood Function ========================== #
    
    QQGPobjFun <- function(par, design, data){
      
      # transform the data and design to matrix to prevent troubles
      designTmp <- as.matrix(design)
      
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
      
      
      #-------------------------------------------------------#
      
      # correlation matrix of categorical combinations (Tmatrix) (by using sphere decomposition)
      Lmatrix <- diag(rep(0, levelNum)); Lmatrix[1, 1] <- 1
      sphereDecLowRowAndLen <- c(2:levelNum)
      LowerTemp <- diag(rep(0, levelNum)); LowerTemp[lower.tri(LowerTemp)] <- theta
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
      if(rcond(Rmatrix) <= 1e-6){
        Rmatrix <- Rmatrix + diag(rep(nugget, totalPtsNum))
      }
      RinvMatrix <- solveCpp(Rmatrix)
      
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
    QQGPnegLogLik <- function(x){
      return(QQGPobjFun(par = x, data = data, design = design))
    }
    
    
    if(is.null(thetaRange)){
      thetaLow <- NULL
      thetaUpp <- NULL
    }else{
      thetaLow <- thetaRange[1]
      thetaUpp <- thetaRange[2]
    }
    
    if(is.null(fiRange)){
      fiLow <- NULL
      fiUpp <- NULL
    }else{
      fiLow <- fiRange[1:(length(fiRange)/2)]
      fiUpp <- fiRange[(length(fiRange)/2+1):length(fiRange)]
    }
    
    
    
    
    # ============================== optimization =============================== #
    # optimize the negLogLik function to find the parameters
    if(is.null(thetaUpp) & is.null(fiUpp)){
      if(is.null(dim(xdiff))){
        upp_bound <- c(rep(log10(-log(0.01)/max(xdiff)), contiDim), rep(0.5, choose(levelNum, 2)))
      }else{
        upp_bound <- c(rep(log10(-log(0.01)/max(apply(xdiff, 1,sum))), contiDim), rep(0.5, choose(levelNum, 2)))
      }
    }else if(is.null(thetaUpp) & !is.null(fiUpp)){
      fiUpp <- fiUpp
      thetaUpp <- 0.5
      upp_bound <- c(fiUpp, rep(thetaUpp, choose(levelNum, 2)))
    }else if(!is.null(thetaUpp) & is.null(fiUpp)){
      fiUpp <- rep(log10(-log(0.01)/max(xdiff)), contiDim)
      upp_bound <- c(fiUpp, rep(thetaUpp, choose(levelNum, 2)))
    }else{
      upp_bound <- c(fiUpp, rep(thetaUpp, choose(levelNum, 2)))
    }
    
    
    if(is.null(thetaLow) & is.null(fiLow)){
      if(is.null(dim(xdiff))){
        low_bound <- c(rep(log10(-log(0.99)/max(xdiff)), contiDim), rep(-0.5, choose(levelNum, 2)))
      }else{
        low_bound <- c(rep(log10(-log(0.99)/max(apply(xdiff, 1,sum))), contiDim), rep(-0.5, choose(levelNum, 2)))
      }
    }else if(is.null(thetaLow) & !is.null(fiLow)){
      fiLow <- fiLow
      thetaLow <- -0.5
      low_bound <- c(fiLow, rep(thetaLow, choose(levelNum, 2)))
    }else if(!is.null(thetaUpp) & is.null(fiUpp)){
      fiLow <- rep(log10(-log(0.99)/max(xdiff)), contiDim)
      low_bound <- c(fiLow, rep(thetaLow, choose(levelNum, 2)))
    }else{
      low_bound <- c(fiLow, rep(thetaLow, choose(levelNum, 2)))
    }
    
    thetaInitParam <- (0.8*upp_bound + 0.2*low_bound)[-c(1:contiDim)]
    
    
    # Use getPSOInfo() to change the PSO options
    alg_setting <- getPSOInfo(nSwarm = 64, maxIter = 10, psoType = "quantum")
    res_c_large <- globpso(objFunc = QQGPnegLogLik, lower = low_bound, upper = upp_bound, PSO_INFO = alg_setting, seed = 1000000,
                           init = c(initParam, thetaInitParam), verbose = F)
    # record the result parameters
    estimatedParameters <- res_c_large$par
    
    
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
    RmatrixHat <- matrix(0, ncol = totalPtsNum, nrow = totalPtsNum)
    for(row in 1:totalPtsNum){
      for(col in 1:totalPtsNum){
        RmatrixHat[row, col] <- HmatrixHat[row, col] * TmatrixHat[levelVecTotal[row], levelVecTotal[col]]
      }
    }
    
    if(rcond(RmatrixHat) <= 1e-6){
      RmatrixHat <- RmatrixHat + diag(rep(nugget, totalPtsNum))
    }
    RinvMatrixHat <- solveCpp(RmatrixHat)
    
  }else{
    
    HmatrixHat <- matrix(0, ncol = totalPtsNum, nrow = totalPtsNum)
    # set the result parameters above as the parameters we will use later
    fiVecHat <- initParam
    Beta <- matrix(fiVecHat, nrow = choose(totalPtsNum, 2), ncol = contiDim, byrow = TRUE)
    Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
    HmatrixHat[rcoord] <- Htemp; HmatrixHat <- HmatrixHat + t(HmatrixHat); HmatrixHat <- exp(-HmatrixHat)
    
    TmatrixHat <- Tmatrix
    # make estimation of inverse correlation matrix of whole data (RmatrixHat) (by using kronecker product)
    # CAREFULLY using kronecker product function, since the order of input should be inverse
    
    RmatrixHat <- matrix(0, ncol = totalPtsNum, nrow = totalPtsNum)
    for(row in 1:totalPtsNum){
      for(col in 1:totalPtsNum){
        RmatrixHat[row, col] <- HmatrixHat[row, col] * TmatrixHat[levelVecTotal[row], levelVecTotal[col]]
      }
    }
    if(rcond(RmatrixHat) <= 1e-6){
      RmatrixHat <- RmatrixHat + diag(rep(nugget, totalPtsNum))
    }
    RinvMatrixHat <- solveCpp(RmatrixHat)
    thetaVecHat <- NULL
  }
  
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
  
  # combine all estimations
  estimatedList <- list(fi = fiVecHat, theta = thetaVecHat, RinvMat = RinvMatrixHat, Tmat = TmatrixHat,
                        design = design, data = data, FVec = FVec, betaVec = betaVecHat, sigmaSq = sigmaSqHat, nugget = nugget,
                        negLogLik = as.numeric(QQGPnegLogLik(estimatedParameters)))
  
  return(estimatedList)
}





###########################################
########## QQGP predict function ##########
###########################################
QQGP.ori.prediction <- function(QQGP.model, newPoint){
  
  # set the result parameters above as the parameters we will use later
  fiVecHat <- QQGP.model$fi
  thetaVecHat <- QQGP.model$theta
  RinvMatrixHat <- QQGP.model$RinvMat
  TmatrixHat <- QQGP.model$Tmat
  nugget <- QQGP.model$nugget
  
  # transform the data and design to matrix to prevent troubles
  designTmp <- data.matrix(as.matrix(QQGP.model$design))

  # define some dimensions about data
  contiDim <- ifelse(is.null(dim(designTmp)), 1, dim(designTmp)[2])
  contiPosNum <- ifelse(is.null(dim(QQGP.model$data)), length(design), dim(QQGP.model$data)[1])
  design <- QQGP.model$data[, 1:contiDim, drop=FALSE]
  levelVec <- unique(QQGP.model$data[, contiDim+1])
  levelVecTotal <- QQGP.model$data[, contiDim+1]
  levelNum <- length(levelVec)
  totalPtsNum <- contiPosNum
  
  betaVecHat <- QQGP.model$betaVec
  FVec <- QQGP.model$FVec
  sigmaSqHat <- QQGP.model$sigmaSq
  
  
  # seperate resonse
  yVec <- QQGP.model$data[, contiDim+2]
  
  # estimation of rediduals vector, this will be use to find the prediction value of newPoint. 
  if(is.null(dim(newPoint))){
    newPointNum <- 1
    rVecHat <- rep(0, totalPtsNum)
    
    predictFVec <- rep(0, levelNum)
    predictFVec[1] <- 1
    predictFVec[which(newPoint[(contiDim+1)] == (1:levelNum))] <- 1
    
    newPoint <- matrix(newPoint, nrow = 1)
    tauMat <- t(matrix(TmatrixHat[as.matrix(expand.grid(QQGP.model$data[,(dim(QQGP.model$data)[2]-1)], newPoint[,dim(newPoint)[2]-1]))], ncol=dim(newPoint)))
    np <- newPoint[,1:contiDim, drop=FALSE]
    RX <- distance(t(t(np)*sqrt(as.numeric(10^fiVecHat))), t(t(design)*sqrt(as.numeric(10^fiVecHat))))
    rVecHat <- exp(-RX)*tauMat
    
    SXX <- 1 + nugget
    yPredict <- as.numeric(predictFVec%*%betaVecHat + (((t(rVecHat)) %*% (RinvMatrixHat)) %*% (yVec-FVec%*%betaVecHat)))
    
  }else{
    newPointNum <- dim(newPoint)[1]
    rVecHat <- matrix(0, ncol = totalPtsNum, nrow = newPointNum)
    
    predictFVec <- matrix(0, ncol = levelNum, nrow = newPointNum)
    predictFVec[,1] <- rep(1, newPointNum)
    for(i in 1:newPointNum){
      predictFVec[i, which(newPoint[i,(contiDim+1)] == (1:levelNum))] <- 1
    }
    
    tauMat <- t(matrix(TmatrixHat[as.matrix(expand.grid(QQGP.model$data[,(dim(QQGP.model$data)[2]-1)], newPoint[,dim(newPoint)[2]-1]))], ncol=dim(newPoint)))
    np <- newPoint[,1:contiDim, drop=FALSE]
    RX <- distance(t(t(np)*sqrt(as.numeric(10^fiVecHat))), t(t(design)*sqrt(as.numeric(10^fiVecHat))))
    rVecHat <- exp(-RX)*tauMat
    
    DXX <- distance(t(t(np)*sqrt(as.numeric(10^fiVecHat))))
    SXX <- exp(-DXX) + diag(nugget, dim(DXX))
    
    yPredict <- as.numeric(predictFVec%*%betaVecHat + ((rVecHat %*% (RinvMatrixHat)) %*% (yVec-FVec%*%betaVecHat)))
    
  }
  
  predSigmaSq <- pmax(0,diag(as.numeric(sigmaSqHat) * (SXX - rVecHat %*% RinvMatrixHat %*% t(rVecHat))))
  
  # return estimated list and prediction function
  return(list(predictValue = yPredict, predSigmaSq = predSigmaSq))
}


##############################################
########## node splitting function  ##########
##############################################
splitNode.general <- function(data, logicVector = NULL, currentPos, leafNumVec, crossCorrMat){
  # print(currentPos)
  options(warn=-1)
  # load current node data with logic vector
  if(is.null(logicVector)){
    logicVec <- rep(TRUE, dim(data)[1])
  }else{
    logicVec <- logicVector
  }
  nodeData <- data
  Y <- nodeData$Y
  cateNum <- sum(unlist(lapply(nodeData, FUN=is.factor)))
  contiVarNum <- ncol(nodeData)-(cateNum+1)
  contiVarVec <- paste0(c("x."), 1:contiVarNum)
  X.conti <- nodeData[,which(names(nodeData) %in% contiVarVec), drop=FALSE]
  # create category vector which is in current data
  cateVarVec <- paste0(c("c."), 1:cateNum)
  cateLevels <- lapply(cateVarVec, FUN = function(c) assign(c, levels(factor(nodeData[logicVec,which(names(nodeData) %in% c)]))))
  names(cateLevels) <- paste0(cateVarVec, ".levels")
  X.cate <- nodeData[logicVec,which(names(nodeData) %in% cateVarVec), drop=FALSE]
  leafNumberVec <- leafNumVec
  # create total possible partitions of current node
  cateNumVec <- sapply(cateLevels, FUN = function(c) return(length(c)))
  cateCorrVec <- lapply(cateLevels, FUN = function(c) return(rep(0, ((2^(length(c)-1)-1)))))
  for(c in paste0(cateVarVec, ".combList")){
    assign(c, list()) 
  }
  for(i in 1:cateNum){
    for(j in 1:(cateNumVec[i]-1)){
      if(j != 0){
        assign(paste0(cateVarVec, ".combList")[i], c(get(paste0(cateVarVec, ".combList")[i]), combn(cateLevels[[i]], j, simplify=FALSE))) 
      }
    }
    assign(paste0(cateVarVec, ".combList")[i], get(paste0(cateVarVec, ".combList")[i])[1:(2^(cateNumVec[i]-1)-1)])
  }
  
  # fit pure GPs
  leafNumberVec <- coerceFun.ori(data = data, design = data[,1:contiVarNum], includeSubCate = T)$c.T
  # note subgroup in this node
  GVec <- unique(leafNumberVec[logicVec])[order(unique(leafNumberVec[logicVec]))]
  
  # compute max cross-correlation
  for(i in 1:cateNum){
    for(j in 1:(2^(cateNumVec[i]-1)-1)){
      if(j != 0){
        # setting somethings will be used in later splitting
        rightNodeLogic <- logicVec & (nodeData[,which(names(nodeData) %in% cateVarVec[i])] %in% unlist(get(paste0(cateVarVec, ".combList")[i])[j]))
        leftNodeLogic <- logicVec & (!rightNodeLogic)
        nodeSplitVar <- cateVarVec[i]
        nodeCateTotalLevels <- cateLevels[[i]]
        rightNodeC.levels <- unlist(get(paste0(cateVarVec, ".combList")[i])[j])
        leftNodeC.levels <- setdiff(nodeCateTotalLevels, rightNodeC.levels)
        # total split logic
        subSplitLogic <- lapply(nodeCateTotalLevels, FUN=function(c) return((nodeData[,(names(data)==nodeSplitVar)]==c) & logicVec))
        names(subSplitLogic) <- nodeCateTotalLevels
        # largest corss correlation
        if((sum(rightNodeLogic) == 0) | (sum(leftNodeLogic) == 0)){
          largestCrossCorr <- 2
        }else{
          crossCorrVec <- numeric(0)
          for(row in unique(leafNumVec[leftNodeLogic])){
            for(col in unique(leafNumVec[rightNodeLogic])){
              crossCorrVec <- c(crossCorrVec, crossCorrMat[row, col])
            }
          }
          largestCrossCorr <- max(crossCorrVec)
        }
        # record the largest corss correlation
        cateCorrVec[[i]][j] <- largestCrossCorr
      }
    }
  }
  minMaxCorr <- min(unlist(lapply(cateCorrVec, FUN = function(x) if(length(x)>0){return(min(x))}else{return(0)})))
  if(length(unique(leafNumVec[logicVec])) > 1){
    splitCateNum <- which(unlist(lapply(cateCorrVec, FUN = function(x) return(any(x == minMaxCorr)))))[1]
    minCate <- cateVarVec[splitCateNum]
    minPos <- which(cateCorrVec[[splitCateNum]] == minMaxCorr)
    minCombRight <- unlist(get(paste0(cateVarVec, ".combList")[which(unlist(lapply(cateCorrVec, FUN = function(x) return(any(x == minMaxCorr)))))][1])[minPos[1]])
    minCombLeft <- setdiff(cateLevels[[splitCateNum]], minCombRight)
    # return node position
    leftNodePosition <- 2*currentPos
    rightNodePosition <- 2*currentPos+1
    # return logic
    returnRightNodeLogic <- logicVec & (nodeData[,which(names(nodeData) %in% cateVarVec[splitCateNum])] %in% get(paste0(cateVarVec, ".combList")[splitCateNum])[[minPos[1]]])
    returnLeftNodeLogic <- logicVec & (!returnRightNodeLogic)
    # return message, which can be useful when prunning the tree
    splitCorrMax <- minMaxCorr
    # return if this is a leaf node logic value
    is.leaf <- FALSE
    # return list
    returnList <- list(splitCateVar = minCate,
                       rightSubCombination = minCombRight, leftSubCombination = minCombLeft,
                       splitCorrMax = splitCorrMax, is.leaf = is.leaf,
                       leftNodePosition = leftNodePosition, rightNodePosition = rightNodePosition,
                       rightSubLogic = returnRightNodeLogic, leftSubLogic = returnLeftNodeLogic)
    return(returnList)
  }else{
    # return if this is a leaf node logic value
    is.leaf <- TRUE
    returnList <- list(is.leaf = is.leaf)
    return(returnList)
  }
}








####################################
########## BCD likelihood ##########
####################################
BCDlikelihood <- function(taus, xdiff, rcoord, covMat, yMat, nugget){
  
  y <- as.vector(yMat)
  n <- length(y) / dim(covMat)[1]
  # make estimation of correlation matrix of continuous positions
  commonRmatrix <- matrix(0, ncol = n, nrow = n)
  # design points differences
  Beta <- matrix(taus, nrow = choose(n, 2), ncol = length(taus), byrow = TRUE)
  Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
  commonRmatrix[rcoord] <- Htemp; commonRmatrix <- commonRmatrix + t(commonRmatrix); commonRmatrix <- exp(-commonRmatrix)
  if(is.null(nugget)){
    commonRmatrix <- commonRmatrix + diag(rep(computeNug(commonRmatrix), n))
  }else{
    commonRmatrix <- commonRmatrix + diag(rep(nugget, n))
  }
  invCommonRmatrix <- solveCpp(commonRmatrix)
  # B matrix (precision matrix) and mu matrix (mean surface matrix)
  oneVec <- matrix(rep(1, n), ncol = 1)
  muVec <- matrix((t(oneVec) %*% invCommonRmatrix %*% yMat) / as.numeric(t(oneVec) %*% invCommonRmatrix %*% oneVec), ncol = 1)
  # here the muMatrix is not identical with BCD algorithm's muMatrix, since this muMatrix is an nk by 1 matrix.
  covMat <- covMat + diag(computeNug(covMat), dim(covMat)[1])
  firstTerm <- n*determinant(covMat, logarithm = TRUE)$modulus
  secondTerm <- dim(covMat)[1]*determinant(commonRmatrix, logarithm = TRUE)$modulus
  thirdTerm <- t(y - kronecker(muVec, oneVec)) %*% kronecker(solveCpp(covMat), invCommonRmatrix) %*% (y - kronecker(muVec, oneVec))
  
  negLikelihood <- as.numeric(firstTerm + secondTerm + thirdTerm)/2
  
  return(list(negLik = negLikelihood, muVec = muVec, commonRmatrix = commonRmatrix))
}






################################
########## initial GP ##########
################################
global.GP.model <- function(design, data, range = NULL, nugget = 0){
  
  # define some dimensions about data
  contiDim <- dim(data)[2]-2
  
  # transform the data to matrix to prevent troubles
  data <- as.matrix(data)
  ptsNum <- dim(data)[1]
  
  # seperate resonse
  design <- data[, c(1:contiDim), drop=FALSE]
  
  # matrix index
  rcoord <- cbind(rep(seq_len(ptsNum-1L), times = rev(seq_len(ptsNum-1L))), 
                  unlist(lapply(X=rev(seq_len(ptsNum-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = ptsNum)))
  # design points differences
  xdiff <- (design[rcoord[, 2L], ] - design[rcoord[, 1L], ])**2
  
  
  GPnegLogLik <- function(x){
    
    levelVec <- unique(data[,(dim(data)[2]-1)])
    
    negLogLikVec <- numeric(length(levelVec))
    
    
    for(levelIndex in levelVec){
      
      dataTmp <- data[(data[,(dim(data)[2]-1)] == levelIndex), -c(dim(data)[2]-1)]
      
      # transform the data to matrix to prevent troubles
      dataTmp <- as.matrix(dataTmp)
      
      # define some dimensions about data
      contiDim <- dim(dataTmp)[2]-1
      
      # ========================== Logliklihood Function ========================== #
      
      GPobjFun <- function(par, data){
        
        # transform the data to matrix to prevent troubles
        data <- as.matrix(data)
        
        # define some dimensions about data
        contiDim <- dim(data)[2]-1
        ptsNum <- dim(data)[1]
        
        # seperate resonse
        design <- data[, c(1:contiDim), drop=FALSE]
        yVec <- data[, contiDim+1]
        
        # seperate parameters vector
        fi <- par
        
        # correlation matrix of continuous positions (Hmatrix)
        Hmatrix <- matrix(0, ncol = ptsNum, nrow = ptsNum)
        # matrix index
        rcoord <- cbind(rep(seq_len(ptsNum-1L), times = rev(seq_len(ptsNum-1L))), 
                        unlist(lapply(X=rev(seq_len(ptsNum-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = ptsNum)))
        # design points differences
        xdiff <- (design[rcoord[, 2L], ] - design[rcoord[, 1L], ])**2
        Beta <- matrix(fi, nrow = choose(ptsNum, 2), ncol = contiDim, byrow = TRUE)
        Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
        Hmatrix[rcoord] <- Htemp; Hmatrix <- Hmatrix + t(Hmatrix); Hmatrix <- exp(-Hmatrix)
        
        # to prevent computational sigular problem, if the condition number of Hmatrix is too small (<= 1e-6)
        # we add a small nugget value to the diagonal line of Hmatrix
        if(rcond(Hmatrix) <= 1e-6){
          Hmatrix <- Hmatrix + diag(rep(nugget, ptsNum))
        }
        
        #-------------------------------------------------------#
        
        # inverse of H.matrix T.matrix and R.matrix (should not have singular problem now)
        HinvMatrix <- solveCpp(Hmatrix)
        
        #-------------------------------------------------------#
        
        # estimators
        FVec <- rep(1L, ptsNum)
        betaVecHat <- solveCpp(t(FVec)%*%HinvMatrix%*%FVec)%*%t(FVec)%*%HinvMatrix%*%yVec
        sigmaSqHat <- (t(yVec-FVec%*%betaVecHat)%*%HinvMatrix%*%(yVec-FVec%*%betaVecHat)) / ptsNum
        
        #-------------------------------------------------------#
        
        # return objective
        subject <- ptsNum*log(sigmaSqHat) + log(det(Hmatrix))
        
        return(subject)
      }
      
      negLogLikVec[levelIndex] <- as.numeric(GPobjFun(par = x, data = dataTmp))
      
    }
    
    return(sum(negLogLikVec))
    
  }

  
  if(is.null(range)){
    if(is.null(dim(xdiff))){
      upp_bound <- rep(log10(-log(0.2)/max(xdiff)), contiDim)
      low_bound <- rep(log10(-log(0.8)/max(xdiff)), contiDim)
    }else{
      upp_bound <- rep(log10(-log(0.2)/max(apply(xdiff, 1,sum))), contiDim)
      low_bound <- rep(log10(-log(0.8)/max(apply(xdiff, 1,sum))), contiDim)
    }
  }else{
    upp_bound <- range[(length(range)/2+1):length(range)]
    low_bound <- range[1:(length(range)/2)]
  }
  
  initParam <- (0.8*upp_bound+0.2*low_bound)
  
  # Use getPSOInfo() to change the PSO options
  alg_setting <- getPSOInfo(nSwarm = 32, maxIter = 5, psoType = "quantum")
  res_c_large <- globpso(objFunc = GPnegLogLik, lower = low_bound, upper = upp_bound,
                         init = initParam, PSO_INFO = alg_setting, seed = 10000,
                         verbose = F)
  estimatedParameters <- res_c_large$par
  
  # combine all estimations
  estimatedList <- list(fi = estimatedParameters)
  
  return(list(estimatedList = estimatedList, bound = c(low_bound, upp_bound), nugget = nugget))
}





####################################################
########## initial GP prediction function ##########
####################################################
global.GP.prediction <- function(GP.model, data, newPoint){
  
  # transform the data to matrix to prevent troubles
  data <- as.matrix(data)
  
  # define some dimensions about data
  contiDim <- dim(data)[2]-1
  ptsNum <- dim(data)[1]
  
  # seperate resonse
  design <- data[, c(1:contiDim), drop=FALSE]
  yVec <- data[, contiDim+1]
  
  # set the result parameters above as the parameters we will use later
  fiVecHat <- GP.model$fi
  nugget <- GP.model$nugget
  
  # make estimation of correlation matrix of continuous positions (HmatrixHat)
  HmatrixHat <- matrix(0, ncol = ptsNum, nrow = ptsNum)
  # matrix index
  rcoord <- cbind(rep(seq_len(ptsNum-1L), times = rev(seq_len(ptsNum-1L))), 
                  unlist(lapply(X=rev(seq_len(ptsNum-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = ptsNum)))
  # design points differences
  xdiff <- (design[rcoord[, 2L], ] - design[rcoord[, 1L], ])**2
  Beta <- matrix(fiVecHat, nrow = choose(ptsNum, 2), ncol = contiDim, byrow = TRUE)
  Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
  HmatrixHat[rcoord] <- Htemp; HmatrixHat <- HmatrixHat + t(HmatrixHat); HmatrixHat <- exp(-HmatrixHat)
  if(rcond(HmatrixHat) <= 1e-6){
    HmatrixHat <- HmatrixHat + diag(rep(1e-10, ptsNum))
  }
  
  # make estimation of inverse correlation matrix of whole data (HmatrixHat)
  HinvMatrixHat <- solveCpp(HmatrixHat)
  
  # parameters estimations of GP model, include beta and sigma
  FVec <- rep(1L,dim(HinvMatrixHat)[1])
  betaVecHat <- (solveCpp(t(FVec)%*%HinvMatrixHat%*%FVec)) %*% t(FVec) %*% HinvMatrixHat %*% yVec
  sigmaSqHat <- (t(yVec-FVec%*%betaVecHat)%*%HinvMatrixHat%*%(yVec-FVec%*%betaVecHat)) / ptsNum
  
  # combine all estimations
  estimatedList <- list(beta = as.numeric(betaVecHat), sigmaSq = as.numeric(sigmaSqHat))
  
  if(is.vector(newPoint)){
    # design matrix of newPoint
    newPoint <- matrix(newPoint, nrow = 1)
  }else{
    newPoint <- as.matrix(newPoint)
  }
  
  # rVecHat function
  tmpFun <- function(x){
    npt <- x
    rVecHat <- rep(0, ptsNum)
    for (k in (1:ptsNum)){
      rVecHat[k] <- exp(-sum(sapply(c(1:contiDim), FUN = function(x) return(((10^fiVecHat[x])*((data[k,x]-as.numeric(npt[x]))^2))))))
    }
    return(rVecHat)
  }
  
  yPredict <- as.numeric(apply(newPoint, MARGIN = 1, FUN = function(x) return(betaVecHat + (((t(tmpFun(x))) %*% (HinvMatrixHat)) %*% (yVec-FVec%*%betaVecHat)))))
  
  # return estimated list and prediction function
  return(list(estimatedList = estimatedList, predictValue = yPredict))
}









##############################################################
########## tree fitting function (corr matrix ver.) ##########
##############################################################
treeFit.corrMat <- function(data, initRange = NULL, parallel = TRUE, BCDloopNum = 3, mcmcIterNum = 50, nugget = NULL){
  nodeData <- data
  cateNum <- sum(unlist(lapply(nodeData, FUN=is.factor)))
  contiVarNum <- ncol(nodeData)-(cateNum+1)
  contiDim <- contiVarNum
  contiVarVec <- paste0(c("x."), 1:contiVarNum)
  cateVarVec <- paste0(c("c."), 1:cateNum)
  # define some dimensions about data
  coerceData <- coerceFun.ori(data = data, design = data[,1:contiVarNum], includeSubCate = T)
  names(coerceData) <- c(contiVarVec, "c.T", "Y")
  cateList <- unique(data[,((contiVarNum+1):(contiVarNum+cateNum)),drop=FALSE])
  
  levelNum <- max(unique(coerceData[,(contiVarNum+1)]))
  ptsNum <- dim(data)[1]
  contiPosNum <- ptsNum / levelNum
  # initial parameter values
  initGP <- global.GP.model(design = data[,1:contiVarNum, drop=FALSE], data = coerceData, range=initRange, nugget=nugget)
  initParams <- initGP$estimatedList
  fiVecHat <- initParams$fi
  initGPRange <- initGP$bound
  
  # detect the union points
  unionDesign <- unique(coerceData[, c(1:contiVarNum), drop=FALSE])
  rownames(unionDesign) <- NULL; names(unionDesign) <- contiVarVec
  if(is.null(dim(unionDesign))){
    len <- length(unionDesign)
  }else{
    len <- dim(unionDesign)[1]
  }
  
  # make estimation of correlation matrix of continuous positions
  commonRmatrix <- matrix(0, ncol = len, nrow = len)
  # matrix index
  rcoord <- cbind(rep(seq_len(len-1L), times = rev(seq_len(len-1L))), 
                  unlist(lapply(X=rev(seq_len(len-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = len)))
  # design points differences
  xdiff <- (unionDesign[rcoord[, 2L], ] - unionDesign[rcoord[, 1L], ])**2
  Beta <- matrix(fiVecHat, nrow = choose(len, 2), ncol = contiVarNum, byrow = TRUE)
  Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
  commonRmatrix[rcoord] <- Htemp; commonRmatrix <- commonRmatrix + t(commonRmatrix); commonRmatrix <- exp(-commonRmatrix)
  if(is.null(nugget)){
    commonRmatrix <- commonRmatrix + diag(rep(computeNug(commonRmatrix), len))
  }else{
    commonRmatrix <- commonRmatrix + diag(rep(nugget, len))
  }
  
  
  invCommonRmatrix <- solveCpp(commonRmatrix)
  oneVec <- matrix(rep(1, len), ncol = 1)
  
  unionLogic <- (len != contiPosNum)
  
  # determine whether to use union design setting
  if(unionLogic){
    # compute union response (interpolation should work here)
    if(parallel){
      
      library(parallel)
      cpu.cores <- detectCores()
      cl = makeCluster(cpu.cores)
      registerDoParallel(cl)
      
      LLLTmp = foreach (lev = 1:levelNum, .combine='c', .packages = c("globpso", "plgp"), .export = c("GP.model", "GP.prediction", "computeNug") ) %dopar% {
        params <- GP.model(coerceData[(coerceData$c.T==lev),-(contiVarNum+1)], initParam = fiVecHat, range = initGPRange, nugget = nugget)
        predTmp <- GP.prediction(GP.model = params, newPoint = unionDesign)
        list(params, predTmp)
      }
      
      paramsList <- LLLTmp[2*(1:levelNum)-1]
      predList <- LLLTmp[2*(1:levelNum)]
      
      library(abind)
      arrayBind <- function(...){
        return(abind(..., along=3))
      }
      
      yImpArr = foreach (lev = 1:levelNum, .combine = "arrayBind") %:%
        foreach (mcmcIter = 1:mcmcIterNum, .combine = 'rbind', .packages = c("mvtnorm"), .export = c("GP.prediction", "computeNug")) %dopar% {
          rmvnorm(1, mean = predList[[lev]]$predictValue, sigma = (predList[[lev]]$predSigmaSq + t(predList[[lev]]$predSigmaSq))/2)
        }
      
      covAndCorrArr = foreach (mcmcIter = 1:mcmcIterNum, .combine = 'c', .export = c("computeNug")) %dopar% {
        muVec <- (t(oneVec) %*% invCommonRmatrix %*% yImpArr[mcmcIter,,]) / as.numeric(t(oneVec) %*% invCommonRmatrix %*% oneVec)
        muMatrix <- matrix(rep(muVec, each = len), ncol = levelNum)
        covMatrix <- t(yImpArr[mcmcIter,,] - muMatrix) %*% invCommonRmatrix %*% (yImpArr[mcmcIter,,] - muMatrix) / contiPosNum
        covMatrix <- (covMatrix + t(covMatrix))/2
        covMatrix <- covMatrix + diag(computeNug(covMatrix), dim(covMatrix)[1])
        corMatrix <- cov2cor(covMatrix)
        corMatrix <- (corMatrix + t(corMatrix))/2
        list(muVec, covMatrix, corMatrix)
      }
      
      muVecMat <- do.call(rbind, covAndCorrArr[3*(1:mcmcIterNum)-2])
      covMatArr <- do.call(arrayBind, covAndCorrArr[3*(1:mcmcIterNum)-1])
      corrMatArr <- do.call(arrayBind, covAndCorrArr[3*(1:mcmcIterNum)])
      
      meanCovMat <- apply(covMatArr, c(1,2), mean)
      meanCorrMat <- apply(corrMatArr, c(1,2), mean)
      
      meanCovMat <- (meanCovMat + t(meanCovMat))/2
      meanCorrMat <- (meanCorrMat + t(meanCorrMat))/2
      
      # continue to BCD loop from optimize theta's and mu's
      predValueMat <- sapply(1:levelNum, FUN = function(x) return(predList[[x]]$predictValue))
      
      tmpTaus <- initParams$fi
      tmpCovMat <- meanCovMat
      
      if(BCDloopNum > 1){
        
        BCDlowBound <- initGPRange[1:(length(initGPRange)/2)]
        BCDuppBound <- initGPRange[(length(initGPRange)/2+1):length(initGPRange)]
        
        for(BCDiter in 1:(BCDloopNum-1)){
          
          # optimize theta's
          optFun <- function(x){return(BCDlikelihood(taus = x, xdiff = xdiff, rcoord = rcoord, covMat = tmpCovMat, yMat = predValueMat, nugget = nugget)$negLik)}
          alg_setting <- getPSOInfo(nSwarm = 32, maxIter = 5, psoType = "quantum")
          res_c_large <- globpso(objFunc = optFun, lower = BCDlowBound, upper = BCDuppBound,
                                 init = tmpTaus, PSO_INFO = alg_setting, seed = 10000,
                                 verbose = F)
          # update theta's
          tmpTaus <- as.numeric(res_c_large$par)
          
          # optimize H mat & T mat
          LLLTmp = foreach (lev = 1:levelNum, .combine='c', .packages = c("globpso", "plgp"), .export = c("computeNug", "GP.model", "GP.prediction", "computeNug") ) %dopar% {
            params <- GP.model(coerceData[(coerceData$c.T==lev),-(contiVarNum+1)], initParam = tmpTaus, range = initGPRange, nugget = nugget)
            predTmp <- GP.prediction(GP.model = params, newPoint = unionDesign)
            list(params, predTmp)
          }
          
          paramsList <- LLLTmp[2*(1:levelNum)-1]
          predList <- LLLTmp[2*(1:levelNum)]
          
          yImpArr = foreach (lev = 1:levelNum, .combine = "arrayBind") %:%
            foreach (mcmcIter = 1:mcmcIterNum, .combine = 'rbind', .packages = c("mvtnorm"), .export = c("computeNug", "GP.prediction", "computeNug")) %dopar% {
              rmvnorm(1, mean = predList[[lev]]$predictValue, sigma = (predList[[lev]]$predSigmaSq + t(predList[[lev]]$predSigmaSq))/2)
            }
          
          
          # update commonRmatrix
          commonRmatrix <- matrix(0, ncol = len, nrow = len)
          Beta <- matrix(tmpTaus, nrow = choose(len, 2), ncol = contiVarNum, byrow = TRUE)
          Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
          commonRmatrix[rcoord] <- Htemp; commonRmatrix <- commonRmatrix + t(commonRmatrix); commonRmatrix <- exp(-commonRmatrix)
          if(is.null(nugget)){
            commonRmatrix <- commonRmatrix + diag(rep(computeNug(commonRmatrix), len))
          }else{
            commonRmatrix <- commonRmatrix + diag(rep(nugget, len))
          }
          invCommonRmatrix <- solveCpp(commonRmatrix)
          
          covAndCorrArr = foreach (mcmcIter = 1:mcmcIterNum, .combine = 'c', .export = c("computeNug")) %dopar% {
            muVec <- (t(oneVec) %*% invCommonRmatrix %*% yImpArr[mcmcIter,,]) / as.numeric(t(oneVec) %*% invCommonRmatrix %*% oneVec)
            muMatrix <- matrix(rep(muVec, each = len), ncol = levelNum)
            covMatrix <- t(yImpArr[mcmcIter,,] - muMatrix) %*% invCommonRmatrix %*% (yImpArr[mcmcIter,,] - muMatrix) / contiPosNum
            covMatrix <- (covMatrix + t(covMatrix))/2
            covMatrix <- covMatrix + diag(computeNug(covMatrix), dim(covMatrix)[1])
            corMatrix <- cov2cor(covMatrix)
            corMatrix <- (corMatrix + t(corMatrix))/2
            list(muVec, covMatrix, corMatrix)
          }
          
          muVecMat <- do.call(rbind, covAndCorrArr[3*(1:mcmcIterNum)-2])
          covMatArr <- do.call(arrayBind, covAndCorrArr[3*(1:mcmcIterNum)-1])
          corrMatArr <- do.call(arrayBind, covAndCorrArr[3*(1:mcmcIterNum)])
          
          meanCovMat <- apply(covMatArr, c(1,2), mean)
          meanCorrMat <- apply(corrMatArr, c(1,2), mean)
          
          meanCovMat <- (meanCovMat + t(meanCovMat))/2
          meanCorrMat <- (meanCorrMat + t(meanCorrMat))/2
          
          # continue to BCD loop from optimize theta's and mu's
          predValueMat <- sapply(1:levelNum, FUN = function(x) return(predList[[x]]$predictValue))
          
          tmpCovMat <- meanCovMat
        } # end of BCD loop
      }
      
      yMatrix <- predValueMat
      fiVecHat <- tmpTaus
      covMatrix <- meanCovMat
      corMatrix <- meanCorrMat
      
      stopCluster(cl)
      
    }else{
      
      
      LLLTmp = foreach (lev = 1:levelNum, .combine='c', .packages = c("globpso", "plgp"), .export = c("GP.model", "GP.prediction", "computeNug") ) %do% {
        params <- GP.model(coerceData[(coerceData$c.T==lev),-(contiVarNum+1)], initParam = fiVecHat, range = initGPRange, nugget = nugget)
        predTmp <- GP.prediction(GP.model = params, newPoint = unionDesign)
        list(params, predTmp)
      }
      
      paramsList <- LLLTmp[2*(1:levelNum)-1]
      predList <- LLLTmp[2*(1:levelNum)]
      
      library(abind)
      arrayBind <- function(...){
        return(abind(..., along=3))
      }
      
      yImpArr = foreach (lev = 1:levelNum, .combine = "arrayBind") %:%
        foreach (mcmcIter = 1:mcmcIterNum, .combine = 'rbind', .packages = c("mvtnorm"), .export = c("computeNug", "GP.prediction", "computeNug")) %do% {
          rmvnorm(1, mean = predList[[lev]]$predictValue, sigma = (predList[[lev]]$predSigmaSq + t(predList[[lev]]$predSigmaSq))/2)
        }
      
      covAndCorrArr = foreach (mcmcIter = 1:mcmcIterNum, .combine = 'c', .export = c("computeNug")) %do% {
        muVec <- (t(oneVec) %*% invCommonRmatrix %*% yImpArr[mcmcIter,,]) / as.numeric(t(oneVec) %*% invCommonRmatrix %*% oneVec)
        muMatrix <- matrix(rep(muVec, each = len), ncol = levelNum)
        covMatrix <- t(yImpArr[mcmcIter,,] - muMatrix) %*% invCommonRmatrix %*% (yImpArr[mcmcIter,,] - muMatrix) / contiPosNum
        covMatrix <- (covMatrix + t(covMatrix))/2
        covMatrix <- covMatrix + diag(computeNug(covMatrix), dim(covMatrix)[1])
        corMatrix <- cov2cor(covMatrix)
        corMatrix <- (corMatrix + t(corMatrix))/2
        list(muVec, covMatrix, corMatrix)
      }
      
      muVecMat <- do.call(rbind, covAndCorrArr[3*(1:mcmcIterNum)-2])
      covMatArr <- do.call(arrayBind, covAndCorrArr[3*(1:mcmcIterNum)-1])
      corrMatArr <- do.call(arrayBind, covAndCorrArr[3*(1:mcmcIterNum)])
      
      meanCovMat <- apply(covMatArr, c(1,2), mean)
      meanCorrMat <- apply(corrMatArr, c(1,2), mean)
      
      meanCovMat <- (meanCovMat + t(meanCovMat))/2
      meanCorrMat <- (meanCorrMat + t(meanCorrMat))/2
      
      # continue to BCD loop from optimize theta's and mu's
      predValueMat <- sapply(1:levelNum, FUN = function(x) return(predList[[x]]$predictValue))
      
      tmpTaus <- initParams$fi
      tmpCovMat <- meanCovMat
      
      if(BCDloopNum > 1){
        
        BCDlowBound <- initGPRange[1:(length(initGPRange)/2)]
        BCDuppBound <- initGPRange[(length(initGPRange)/2+1):length(initGPRange)]
        
        for(BCDiter in 1:(BCDloopNum-1)){
          
          # optimize theta's
          optFun <- function(x){return(BCDlikelihood(taus = x, xdiff = xdiff, rcoord = rcoord, covMat = tmpCovMat, yMat = predValueMat, nugget = nugget)$negLik)}
          alg_setting <- getPSOInfo(nSwarm = 32, maxIter = 5, psoType = "quantum")
          res_c_large <- globpso(objFunc = optFun, lower = BCDlowBound, upper = BCDuppBound,
                                 init = tmpTaus, PSO_INFO = alg_setting, seed = 10000,
                                 verbose = F)
          
          # update theta's
          tmpTaus <- as.numeric(res_c_large$par)
          
          # optimize H mat & T mat
          LLLTmp = foreach (lev = 1:levelNum, .combine='c', .packages = c("globpso", "plgp"), .export = c("GP.model", "GP.prediction", "computeNug") ) %do% {
            params <- GP.model(coerceData[(coerceData$c.T==lev),-(contiVarNum+1)], initParam = tmpTaus, range = initGPRange, nugget = nugget)
            predTmp <- GP.prediction(GP.model = params, newPoint = unionDesign)
            list(params, predTmp)
          }
          
          paramsList <- LLLTmp[2*(1:levelNum)-1]
          predList <- LLLTmp[2*(1:levelNum)]
          
          yImpArr = foreach (lev = 1:levelNum, .combine = "arrayBind") %:%
            foreach (mcmcIter = 1:mcmcIterNum, .combine = 'rbind', .packages = c("mvtnorm"), .export = c("GP.prediction", "computeNug")) %do% {
              rmvnorm(1, mean = predList[[lev]]$predictValue, sigma = (predList[[lev]]$predSigmaSq + t(predList[[lev]]$predSigmaSq))/2)
            }
          
          # update commonRmatrix
          commonRmatrix <- matrix(0, ncol = len, nrow = len)
          Beta <- matrix(tmpTaus, nrow = choose(len, 2), ncol = contiVarNum, byrow = TRUE)
          Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
          commonRmatrix[rcoord] <- Htemp; commonRmatrix <- commonRmatrix + t(commonRmatrix); commonRmatrix <- exp(-commonRmatrix)
          if(is.null(nugget)){
            commonRmatrix <- commonRmatrix + diag(rep(computeNug(commonRmatrix), len))
          }else{
            commonRmatrix <- commonRmatrix + diag(rep(nugget, len))
          }
          invCommonRmatrix <- solveCpp(commonRmatrix)
          
          covAndCorrArr = foreach (mcmcIter = 1:mcmcIterNum, .combine = 'c', .export = c("computeNug")) %do% {
            muVec <- (t(oneVec) %*% invCommonRmatrix %*% yImpArr[mcmcIter,,]) / as.numeric(t(oneVec) %*% invCommonRmatrix %*% oneVec)
            muMatrix <- matrix(rep(muVec, each = len), ncol = levelNum)
            covMatrix <- t(yImpArr[mcmcIter,,] - muMatrix) %*% invCommonRmatrix %*% (yImpArr[mcmcIter,,] - muMatrix) / contiPosNum
            covMatrix <- (covMatrix + t(covMatrix))/2
            covMatrix <- covMatrix + diag(computeNug(covMatrix), dim(covMatrix)[1])
            corMatrix <- cov2cor(covMatrix)
            corMatrix <- (corMatrix + t(corMatrix))/2
            list(muVec, covMatrix, corMatrix)
          }
          
          muVecMat <- do.call(rbind, covAndCorrArr[3*(1:mcmcIterNum)-2])
          covMatArr <- do.call(arrayBind, covAndCorrArr[3*(1:mcmcIterNum)-1])
          corrMatArr <- do.call(arrayBind, covAndCorrArr[3*(1:mcmcIterNum)])
          
          meanCovMat <- apply(covMatArr, c(1,2), mean)
          meanCorrMat <- apply(corrMatArr, c(1,2), mean)
          
          meanCovMat <- (meanCovMat + t(meanCovMat))/2
          meanCorrMat <- (meanCorrMat + t(meanCorrMat))/2
          
          # continue to BCD loop from optimize theta's and mu's
          predValueMat <- sapply(1:levelNum, FUN = function(x) return(predList[[x]]$predictValue))
          
          tmpCovMat <- meanCovMat
        } # end of BCD loop
      }
      
      yMatrix <- predValueMat
      fiVecHat <- tmpTaus
      corMatrix <- meanCorrMat
      covMatrix <- meanCovMat
    } # end of parallel, MCMC and BCD loop
    
    unionData <- cbind(unionDesign[rep(1:len, levelNum),], cateList[rep(1:levelNum, each=len),], as.vector(predValueMat))
    rownames(unionData) <- NULL; names(unionData) <- names(data)
    unionData <- as.data.frame(unionData)
    
    coerceUnionData <- coerceFun.ori(data = unionData, design = unionData[,1:contiVarNum], includeSubCate = T)
    names(coerceUnionData) <- c(contiVarVec, "c.T", "Y")
    
    muVec <- matrix((t(oneVec) %*% invCommonRmatrix %*% yMatrix) / as.numeric(t(oneVec) %*% solveCpp(commonRmatrix) %*% oneVec), ncol = levelNum)
    
  }else{
    
    unionData <- data
    coerceUnionData <- coerceFun.ori(data = unionData, design = unionData[,1:contiVarNum], includeSubCate = T)
    names(coerceUnionData) <- c(contiVarVec, "c.T", "Y")
    
    yMat <- matrix(unionData$Y, ncol = levelNum)
    
    muVec <- (t(oneVec) %*% invCommonRmatrix %*% yMat) / as.numeric(t(oneVec) %*% invCommonRmatrix %*% oneVec)
    muMatrix <- matrix(rep(muVec, each = len), ncol = levelNum)
    covMatrix <- t(yMat - muMatrix) %*% invCommonRmatrix %*% (yMat - muMatrix) / contiPosNum
    covMatrix <- (covMatrix + t(covMatrix))/2
    covMatrix <- covMatrix + diag(computeNug(covMatrix), dim(covMatrix)[1])
    corMatrix <- cov2cor(covMatrix)
    corMatrix <- (corMatrix + t(corMatrix))/2
    
    # continue to BCD loop from optimize theta's and mu's
    tmpTaus <- initParams$fi
    tmpCovMat <- covMatrix
    
    if(BCDloopNum > 1){
      
      BCDlowBound <- initGPRange[1:(length(initGPRange)/2)]
      BCDuppBound <- initGPRange[(length(initGPRange)/2+1):length(initGPRange)]
      
      for(BCDiter in 1:(BCDloopNum-1)){
        
        # optimize theta's
        optFun <- function(x){return(BCDlikelihood(taus = x, xdiff = xdiff, rcoord = rcoord, covMat = tmpCovMat, yMat = yMat, nugget = nugget)$negLik)}
        alg_setting <- getPSOInfo(nSwarm = 32, maxIter = 5, psoType = "quantum")
        res_c_large <- globpso(objFunc = optFun, lower = BCDlowBound, upper = BCDuppBound,
                               init = tmpTaus, PSO_INFO = alg_setting, seed = 10000,
                               verbose = F)
        # update theta's
        tmpTaus <- as.numeric(res_c_large$par)
        # optimize H mat & T mat
        # update commonRmatrix
        commonRmatrix <- matrix(0, ncol = len, nrow = len)
        Beta <- matrix(tmpTaus, nrow = choose(len, 2), ncol = contiVarNum, byrow = TRUE)
        Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
        commonRmatrix[rcoord] <- Htemp; commonRmatrix <- commonRmatrix + t(commonRmatrix); commonRmatrix <- exp(-commonRmatrix)
        if(is.null(nugget)){
          commonRmatrix <- commonRmatrix + diag(rep(computeNug(commonRmatrix), len))
        }else{
          commonRmatrix <- commonRmatrix + diag(rep(nugget, len))
        }
        invCommonRmatrix <- solveCpp(commonRmatrix)
        
        muVec <- (t(oneVec) %*% invCommonRmatrix %*% yMat) / as.numeric(t(oneVec) %*% invCommonRmatrix %*% oneVec)
        muMatrix <- matrix(rep(muVec, each = len), ncol = levelNum)
        covMatrix <- t(yMat - muMatrix) %*% invCommonRmatrix %*% (yMat - muMatrix) / contiPosNum
        covMatrix <- (covMatrix + t(covMatrix))/2
        covMatrix <- covMatrix + diag(computeNug(covMatrix), dim(covMatrix)[1])
        corMatrix <- cov2cor(covMatrix)
        corMatrix <- (corMatrix + t(corMatrix))/2
        tmpCovMat <- covMatrix
      } # end of BCD loop
    }
    
    yMatrix <- yMat
    fiVecHat <- tmpTaus
    corMatrix <- corMatrix
    covMatrix <- covMatrix
  }
  
  ########################
  
  # to do vec and node vec
  toDoPosVec <- c(1)
  nodePosVec <- c(1)
  # return information function - node text
  nodeTextFun <- function(x){
    splitFit <- get(paste0("node.", x, ".split"))
    if(splitFit$is.leaf){
      return("leaf")
    }else{
      return(splitFit$splitCateVar)
    }
  }
  # return information function - parent node
  parentNodeFun <- function(x){
    if(x == 1){
      return(0)
    }else if((x>1) && (x %% 2==0)){
      return(x/2)
    }else{
      return(x/2-0.5)
    }
  }
  # return information function - parent split levels
  parentSplitLevelFun <- function(x){
    if(x > 1){
      if((x%%2)==0){
        parentPos <- x/2
        return(paste0(paste(get(paste0("node.", parentPos, ".split"))$leftSubCombination, collapse = ' or ')))
      }else{
        parentPos <- x/2 - 0.5
        return(paste0(paste(get(paste0("node.", parentPos, ".split"))$rightSubCombination, collapse = ' or ')))
      }
    }else{
      return("root")
    }
  }
  # return information function - parent split logic
  parentSplitLogicFun <- function(x){
    if(x > 1){
      if((x%%2)==0){
        parentPos <- x/2
        return(get(paste0("node.", parentPos, ".split"))$leftSubLogic)
      }else{
        parentPos <- x/2 - 0.5
        return(get(paste0("node.", parentPos, ".split"))$rightSubLogic)
      }
    }else{
      return(rep(TRUE, dim(data)[1]))
    }
  }
  # do tree
  minMaxCorr <- numeric(0)
  leafNumberVec <- coerceFun.ori(data = data, design = data[,1:contiVarNum], includeSubCate = T)$c.T
  assign(paste0("node.0.split"), NULL)
  while(!is.null(toDoPosVec)){
    tmpTpDoPosVec <- toDoPosVec
    toDoPosVec <- NULL
    for(pos in tmpTpDoPosVec){
      #print(pos)
      if(pos > 1){
        if((pos%%2)==0){
          parentPos <- pos/2
          assign(paste0("node.", pos, ".split"), splitNode.general(data = data, logicVector = get(paste0("node.", parentPos, ".split"))$leftSubLogic, currentPos = pos,
                                                                   leafNumVec = leafNumberVec, crossCorrMat = corMatrix))
          if(get(paste0("node.", pos, ".split"))$is.leaf == TRUE){
            minMaxCorr <- c(minMaxCorr, 2)
          }else{
            tmp <- get(paste0("node.", pos, ".split"))$splitCorrMax
            minMaxCorr <- c(minMaxCorr, tmp)
          }
        }else{
          parentPos <- pos/2 - 0.5
          assign(paste0("node.", pos, ".split"), splitNode.general(data = data, logicVector = get(paste0("node.", parentPos, ".split"))$rightSubLogic, currentPos = pos,
                                                                   leafNumVec = leafNumberVec, crossCorrMat = corMatrix))
          if(get(paste0("node.", pos, ".split"))$is.leaf == TRUE){
            minMaxCorr <- c(minMaxCorr, 2)
          }else{
            tmp <- get(paste0("node.", pos, ".split"))$splitCorrMax
            minMaxCorr <- c(minMaxCorr, tmp)
          }
        }
      }else{
        assign(paste0("node.", pos, ".split"), splitNode.general(data = data, currentPos = pos, leafNumVec = leafNumberVec, crossCorrMat = corMatrix))
        if(get(paste0("node.", pos, ".split"))$is.leaf){
          minMaxCorr <- c(minMaxCorr, 2)
        }else{
          tmp <- get(paste0("node.", pos, ".split"))$splitCorrMax
          minMaxCorr <- c(minMaxCorr, tmp)
        }
      }
      if(!get(paste0("node.", pos, ".split"))$is.leaf){
        toDoPosVec <- c(toDoPosVec, get(paste0("node.", pos, ".split"))$leftNodePosition, get(paste0("node.", pos, ".split"))$rightNodePosition)
      }
    }
    nodePosVec <- c(nodePosVec, toDoPosVec)
  }
  nodeText <- sapply(nodePosVec, FUN=nodeTextFun)
  parentNode <- sapply(nodePosVec, FUN=parentNodeFun)
  parentSplitLevel <- sapply(nodePosVec, FUN=parentSplitLevelFun)
  parentSplitLogic <- lapply(nodePosVec, FUN=parentSplitLogicFun)
  classTable <- do.call(rbind, lapply(which(nodeText == "leaf"), FUN=function(x){
    return(cbind.data.frame(unique(data[parentSplitLogic[[x]],which(names(nodeData) %in% cateVarVec)]), unique(coerceData[parentSplitLogic[[x]], "c.T"]),
                            ifelse(length(cateVarVec) == 1, (length(unique(data[parentSplitLogic[[x]],which(names(nodeData) %in% cateVarVec)]))[1] == 1), (dim(unique(data[parentSplitLogic[[x]],which(names(nodeData) %in% cateVarVec)]))[1] == 1)),
                            nodePosVec[x], stringsAsFactors = FALSE))}))
  
  rownames(classTable) <- NULL
  colnames(classTable) <- c(cateVarVec, "coerceIndex", "pure", "node")
  
  
  unionParentSplitLogic <- list()
  for(i in 1:length(parentSplitLogic)){
    unionParentSplitLogic <- c(unionParentSplitLogic, list(coerceUnionData$c.T %in% unique(coerceData$c.T[parentSplitLogic[[i]]])))
  }
  
  returnList <- list(nodePosition = nodePosVec, parentNodePosition = parentNode, nodeSplitVar = nodeText,
                     nodeMinMaxCorr = minMaxCorr, parentSplitLevel = parentSplitLevel, parentSplitLogic = parentSplitLogic,
                     unionParentSplitLogic = unionParentSplitLogic,
                     coerceData = coerceData, unionData = unionData, unionDesign = unionDesign, coerceUnionData = coerceUnionData,
                     classTable = classTable, tauMatrix = corMatrix, tauCovMatrix = covMatrix,
                     muVec = muVec, fiParams = fiVecHat, initGPRange = initGPRange, nugget = nugget)
  return(returnList)
}










####################################
########## LOOCV function ##########
####################################
LOOCVFun <- function(data, param, muHat, tauMatrix = NULL){
  # transform the data to matrix to prevent troubles
  data <- as.matrix(data)
  
  # define some dimensions about data
  contiDim <- dim(data)[2]-1
  ptsNum <- dim(data)[1]
  
  # seperate resonse
  design <- data[, c(1:contiDim), drop=FALSE]
  yVec <- data[, contiDim+1]
  
  # set the result parameters above as the parameters we will use later
  fiVecHat <- param
  
  if(is.null(tauMatrix)){
    kernelDim <- ptsNum
  }else{
    kernelDim <- ptsNum/(dim(tauMatrix)[1])
  }
  
  # make estimation of correlation matrix of continuous positions (HmatrixHat)
  HmatrixHat <- matrix(0, ncol = kernelDim, nrow = kernelDim)
  # matrix index
  rcoord <- cbind(rep(seq_len(kernelDim-1L), times = rev(seq_len(kernelDim-1L))),
                  unlist(lapply(X=rev(seq_len(kernelDim-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = kernelDim)))
  # design points differences
  xdiff <- (design[rcoord[, 2L], ] - design[rcoord[, 1L], ])**2
  Beta <- matrix(fiVecHat, nrow = choose(kernelDim, 2), ncol = contiDim, byrow = TRUE)
  Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
  HmatrixHat[rcoord] <- Htemp; HmatrixHat <- HmatrixHat + t(HmatrixHat); HmatrixHat <- exp(-HmatrixHat)
  HmatrixHat <- HmatrixHat + diag(computeNug(HmatrixHat), kernelDim)
  
  # make estimation of inverse correlation matrix of whole data (HmatrixHat)
  HinvMatrixHat <- solveCpp(HmatrixHat)
  
  if(!is.null(tauMatrix)){
    tauInvMatrix <- solveCpp(tauMatrix+diag(computeNug(tauMatrix), nrow(tauMatrix)))
    invKernelMatrix <- kronecker(tauInvMatrix, HinvMatrixHat)
  }else{
    invKernelMatrix <- HinvMatrixHat
  }
  
  LamMatrix <- diag(diag(1/invKernelMatrix))
  
  
  # parameters estimations of GP model, include beta and sigma
  FVec <- rep(1L, dim(invKernelMatrix)[1])
  muVecHat <- rep(muHat, each = kernelDim)
  
  return(sum((LamMatrix%*%invKernelMatrix%*%(yVec-muVecHat))^2))
}



########################################
########## combining function ##########
########################################
combineFun <- function(tree, data)
{
  # pre-define data informations
  nodeData <- data
  cateNum <- sum(unlist(lapply(nodeData, FUN=is.factor)))
  contiVarNum <- ncol(nodeData)-(cateNum+1)
  contiDim <- contiVarNum
  
  # pre-define original tree informations
  nodePositionVec <- tree$nodePosition
  parentNodePositionVec <- tree$parentNodePosition
  nodeSplitVar <- tree$nodeSplitVar
  nodeMinMaxCorrVec <- tree$nodeMinMaxCorr
  parentSplitLevelVec <- tree$parentSplitLevel
  parentSplitLogicList <- tree$parentSplitLogic
  unionParentSplitLogicList <- tree$unionParentSplitLogic
  classTable <- tree$classTable
  crossCorrMatrix <- tree$tauMatrix
  crossCovMatrix <- tree$tauCovMatrix
  leafMuVec <- tree$muVec
  commonFi <- tree$fiParams
  unionData <- tree$unionData
  coerceData <- tree$coerceData
  unionDesign <- tree$unionDesign
  coerceUnionData <- tree$coerceUnionData
  initGPRange <- tree$initGPRange
  nugget <- tree$nugget
  
  rootNodeContainLeaves <- unique(parentNodePositionVec[nodeSplitVar == "leaf"])[order(unique(parentNodePositionVec[nodeSplitVar == "leaf"]))]
  candidateRoot <- rootNodeContainLeaves[table(parentNodePositionVec[nodeSplitVar == "leaf"]) == 2]
  combineRoot <- numeric(0)
  combineLeaves <- numeric(0)
  for(root in candidateRoot){
    
    combineLeafNodePosition <- nodePositionVec[parentNodePositionVec == root]
    combineLeftLeafNodePosition <- nodePositionVec[parentNodePositionVec == root][1]
    combineRightLeafNodePosition <- nodePositionVec[parentNodePositionVec == root][2]
    combineCoerceIndex <- classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root]][order(classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root]])]
    combineCoerceIndexLeft <- classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root][1]][order(classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root][1]])]
    combineCoerceIndexRight <- classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root][2]][order(classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root][2]])]
    combineCrossCorrMatrix <- crossCorrMatrix[combineCoerceIndex, combineCoerceIndex]
    combineCrossCorrMatrixLeft <- as.matrix(crossCorrMatrix[combineCoerceIndexLeft, combineCoerceIndexLeft])
    combineCrossCorrMatrixRight <- as.matrix(crossCorrMatrix[combineCoerceIndexRight, combineCoerceIndexRight])
    
    combineLOOCV <- LOOCVFun(coerceUnionData[unionParentSplitLogicList[which(nodePositionVec %in% root)][[1]],c(1:contiVarNum, dim(coerceUnionData)[2])],
                             tree$fiParams, leafMuVec[combineCoerceIndex], as.matrix(combineCrossCorrMatrix))
    separateLOOCV <- LOOCVFun(coerceUnionData[unionParentSplitLogicList[nodePositionVec == combineLeftLeafNodePosition][[1]],c(1:contiVarNum, dim(coerceUnionData)[2])],
                              tree$fiParams, leafMuVec[combineCoerceIndexLeft], combineCrossCorrMatrixLeft) +
      LOOCVFun(coerceUnionData[unionParentSplitLogicList[nodePositionVec == combineRightLeafNodePosition][[1]],c(1:contiVarNum, dim(coerceUnionData)[2])],
               tree$fiParams, leafMuVec[combineCoerceIndexRight], combineCrossCorrMatrixRight)
    
    if(combineLOOCV < separateLOOCV){
      combineRoot <- c(combineRoot, root)
      combineLeaves <- c(combineLeaves, combineLeafNodePosition)
    }
  }
  
  
  if(length(combineRoot) != 0)
  {
    removePosition <- (parentNodePositionVec %in% combineRoot)
    # prune the tree (specifically, combine these two nodes)
    nodePositionVecTmp <- nodePositionVec[!removePosition]
    parentNodePositionVec <- parentNodePositionVec[!removePosition]
    nodeSplitVar[nodePositionVec %in% combineRoot] <- "leaf"
    nodeSplitVar <- nodeSplitVar[!removePosition]
    nodeMinMaxCorrVec[nodePositionVec %in% combineRoot] <- 2
    nodeMinMaxCorrVec <- nodeMinMaxCorrVec[!removePosition]
    parentSplitLevelVec <- parentSplitLevelVec[!removePosition]
    parentSplitLogicList <- parentSplitLogicList[!removePosition]
    unionParentSplitLogicList <- unionParentSplitLogicList[!removePosition]
    classTable[(classTable$node %in% combineLeaves), "pure"] <- FALSE
    nodePositionVec <- nodePositionVecTmp
    
    for(r in 1:length(combineRoot)){
      classTable[(classTable$node %in% combineLeaves[(2*(r-1)+1):(2*r)]), "node"] <- combineRoot[r]
    }
    
    returnList <- list(nodePosition = nodePositionVec, parentNodePosition = parentNodePositionVec,
                       nodeSplitVar = nodeSplitVar, nodeMinMaxCorr = nodeMinMaxCorrVec,
                       parentSplitLevel = parentSplitLevelVec, parentSplitLogic = parentSplitLogicList,
                       unionParentSplitLogic = unionParentSplitLogicList,
                       coerceData = coerceData, unionData = unionData, unionDesign = unionDesign, coerceUnionData = coerceUnionData,
                       classTable = classTable, tauMatrix = crossCorrMatrix, tauCovMatrix = crossCovMatrix, #invCommonRmatrix = invCommonRmatrix,
                       muVec = leafMuVec, fiParams = commonFi, initGPRange = initGPRange, nugget = nugget)
    
    return(list(prune = TRUE, tree = returnList))
  }
  else
  {
    # if combining LOOCV is not smaller than separating LOOCV, return the original tree.
    return(list(prune = FALSE, tree = tree))
  }
}





######################################
########## pruning function ##########
######################################
pruneFun <- function(tree, data)
{
  tmpTree <- combineFun(tree, data)
  while(tmpTree$prune)
  {
    if(length(unique(tmpTree$tree$classTable$node)) == 1)
      break
    tmpTree <- combineFun(tmpTree$tree, data)
  }
  returnTree <- tmpTree$tree
  return(returnTree)
}








###################################################
########## leaf GP/QQGP fitting function ##########
###################################################
leafGPFun <- function(data, tree, NRF = FALSE, GPRange = NULL, QQfiRange = NULL, QQthetaRange = NULL){
  
  # pre-set data
  nodeData <- data
  Y <- nodeData$Y
  cateNum <- sum(unlist(lapply(nodeData, FUN=is.factor)))
  contiVarNum <- ncol(nodeData)-(cateNum+1)
  initParams <- tree$fiParams
  nugget <- tree$nugget
  
  if(is.null(GPRange)){
    GPRange <- tree$initGPRange
  }
  
  # leaf node number
  leafNodeVec <- tree$classTable$node
  nodePosition <- tree$nodePosition
  
  # create null lists
  GPList <- c()
  GPdataList <- c()
  QQGPList <- c()
  QQGPdesignList <- c()
  QQGPdataList <- c()
  
  if(!NRF){
    
    # fit GP or QQGP model depending whether the leaf is pure or not
    for(num in unique(leafNodeVec)){
      if(tree$classTable$pure[tree$classTable$node == num][1]){
        # if it is a pure leaf, subset only the continuous part of data and the response
        assign(paste0("node.", num, ".data"), data[tree$parentSplitLogic[nodePosition == num][[1]], c(1:contiVarNum, (contiVarNum+cateNum+1))])
        GPdataList <- c(GPdataList, list(get(paste0("node.", num, ".data"))))
        # since it is a pure leaf, use GP to fit
        assign(paste0("node.", num, ".model"), GP.model(data = get(paste0("node.", num, ".data")), initParam = as.numeric(initParams), range = GPRange, nugget = nugget))
        GPList <- c(GPList, list(get(paste0("node.", num, ".model"))))
      }else{
        # if it is not a pure leaf, subset all parts of data and the response
        assign(paste0("node.", num, ".data"), data[tree$parentSplitLogic[nodePosition == num][[1]],])
        QQGPdataList <- c(QQGPdataList, list(get(paste0("node.", num, ".data"))))
        assign(paste0("node.", num, ".design"), unique(data[tree$parentSplitLogic[nodePosition == num][[1]], c(1:contiVarNum), drop=FALSE]))
        QQGPdesignList <- c(QQGPdesignList, list(get(paste0("node.", num, ".design"))))
        # since it is not a pure leaf, use QQGP to fit
        assign(paste0("node.", num, ".model"), QQGP.ori.model(data = coerceFun.ori(data = data[tree$parentSplitLogic[nodePosition == num][[1]],], design = unique(data[tree$parentSplitLogic[nodePosition == num][[1]],c(1:contiVarNum),drop=FALSE]), includeSubCate = T),
                                                              design = get(paste0("node.", num, ".design")), initParam = as.numeric(initParams), fiRange = QQfiRange, thetaRange = QQthetaRange, nugget = nugget)) 
        QQGPList <- c(QQGPList, list(get(paste0("node.", num, ".model"))))
      }
    }
    
    
  }else{
    
    # fit GP or QQGP model depending whether the leaf is pure or not
    for(num in unique(leafNodeVec)){
      if(tree$classTable$pure[tree$classTable$node == num][1]){
        # if it is a pure leaf, subset only the continuous part of data and the response
        assign(paste0("node.", num, ".data"), data[tree$parentSplitLogic[nodePosition == num][[1]], c(1:contiVarNum, (contiVarNum+cateNum+1))])
        GPdataList <- c(GPdataList, list(get(paste0("node.", num, ".data"))))
        # since it is a pure leaf, use GP to fit
        assign(paste0("node.", num, ".model"), GP.model(data = get(paste0("node.", num, ".data")), initParam = as.numeric(initParams)), nugget = nugget)
        GPList <- c(GPList, list(get(paste0("node.", num, ".model"))))
      }else{
        # if it is not a pure leaf, subset all parts of data and the response
        assign(paste0("node.", num, ".data"), data[tree$parentSplitLogic[nodePosition == num][[1]],])
        QQGPdataList <- c(QQGPdataList, list(get(paste0("node.", num, ".data"))))
        assign(paste0("node.", num, ".design"), unique(data[tree$parentSplitLogic[nodePosition == num][[1]], c(1:contiVarNum), drop=FALSE]))
        QQGPdesignList <- c(QQGPdesignList, list(get(paste0("node.", num, ".design"))))
        # since it is not a pure leaf, use QQGP to fit
        assign(paste0("node.", num, ".model"), QQGP.ori.model(data = coerceFun.ori(data = data[tree$parentSplitLogic[nodePosition == num][[1]],], design = unique(data[tree$parentSplitLogic[nodePosition == num][[1]],c(1:contiVarNum),drop=FALSE]), includeSubCate = T),
                                                              design = get(paste0("node.", num, ".design")), initParam = as.numeric(initParams), Tmatrix = tree$tauMatrix[(tree$classTable$node %in% num), (tree$classTable$node %in% num)], nugget = nugget))
        QQGPList <- c(QQGPList, list(get(paste0("node.", num, ".model"))))
      }
    }
    
  }

  # rename the list names
  if(!is.null(GPList)){
    names(GPList) <- paste0("node.", tree$classTable$node[tree$classTable$pure], ".model")
  }
  if(!is.null(GPdataList)){
    names(GPdataList) <- paste0("node.", tree$classTable$node[tree$classTable$pure], ".data")
  }
  if(!is.null(QQGPList)){
    names(QQGPList) <- paste0("node.", unique(tree$classTable$node[!tree$classTable$pure]), ".model")
  }
  if(!is.null(QQGPdesignList)){
    names(QQGPdesignList) <- paste0("node.", unique(tree$classTable$node[!tree$classTable$pure]), ".design")
  }
  if(!is.null(QQGPdataList)){
    names(QQGPdataList) <- paste0("node.", unique(tree$classTable$node[!tree$classTable$pure]), ".data")
  }
  
  return(list(GP.model = GPList, GP.data = GPdataList, QQGP.model = QQGPList, QQGP.design = QQGPdesignList, QQGP.data = QQGPdataList))
}






###################################################
########## leaf GP/QQGP predict function ##########
###################################################
predictFun <- function(newPoint, tree, leaf.GP, parallel=F){
  
  # pre-set data
  nodeData <- newPoint
  cateNum <- sum(unlist(lapply(nodeData, FUN=is.factor)))
  contiVarNum <- ncol(nodeData)-(cateNum+1)
  # identify which node should the new point belong
  nodeTable <- tree$classTable[,c(1:cateNum, dim(tree$classTable)[2])]
  suppressMessages(
    node <- left_join(newPoint[,c((contiVarNum+1):(contiVarNum+cateNum)),drop=F], nodeTable)$node
  )
  pure <- sapply(node, FUN = function(x) return(tree$classTable$pure[tree$classTable$node == x][1]))
  predictValueVec <- numeric(dim(newPoint)[1])
  predictVarianceVec <- numeric(dim(newPoint)[1])
  
  if(parallel){
    library(parallel)
    cpu.cores <- detectCores()
    cl = makeCluster(cpu.cores)
    registerDoParallel(cl)
    tmplist = foreach (nnn = unique(node), .combine = 'c', .packages = c("plgp", "dplyr"),
                       .export = c("tree", "node", "pure", "contiVarNum", "nodeData", "GP.prediction", "QQGP.ori.prediction", "coerceFun.ori")) %dopar% {
                         # predict value depending whether the leaf is pure
                         if(pure[node == nnn][1]){
                           GP.predict.model <- leaf.GP$GP.model[[which(tree$classTable$node[tree$classTable$pure] == nnn)]]
                           GP.predict.data <- leaf.GP$GP.data[[which(tree$classTable$node[tree$classTable$pure] == nnn)]]
                           predictPoint <- newPoint[,c(1:contiVarNum, ncol(nodeData)), drop=FALSE]
                           predTmp <- GP.prediction(GP.model = GP.predict.model, newPoint = predictPoint)
                           predictValue <- predTmp$predictValue
                           predictVariance <- predTmp$predSigmaSq
                         }else{
                           QQGP.predict.model <- leaf.GP$QQGP.model[[which(unique(tree$classTable$node[!tree$classTable$pure]) == nnn)]]
                           QQGP.predict.data <- leaf.GP$QQGP.data[[which(unique(tree$classTable$node[!tree$classTable$pure]) == nnn)]]
                           QQGP.predict.design <- leaf.GP$QQGP.design[[which(unique(tree$classTable$node[!tree$classTable$pure]) == nnn)]]
                           compareTable <- cbind(QQGP.predict.data[,((contiVarNum+1):(contiVarNum+cateNum)), drop = FALSE], coerceFun.ori(data = QQGP.predict.data, design = QQGP.predict.design, includeSubCate = T)$c.T)
                           names(compareTable) <- c(names(compareTable)[1:cateNum], "c.T")
                           compareTable <- unique(compareTable)
                           suppressMessages(QQGPcomb <- left_join(newPoint[(node == nnn), c((contiVarNum+1):(contiVarNum+cateNum)), drop=F], compareTable)$c.T)
                           predictPoint <- cbind(newPoint[which(node == nnn), c(1:contiVarNum)], QQGPcomb, newPoint[which(node == nnn), (contiVarNum+cateNum+1)])
                           predTmp <- QQGP.ori.prediction(QQGP.model = QQGP.predict.model, newPoint = predictPoint)
                           predictValue <- predTmp$predictValue
                           predictVariance <- predTmp$predSigmaSq
                         }
                         list(predictValue, predictVariance)
                       }
    for( no in 1:length(unique(node)) ){
      predictValueVec[node == unique(node)[no]] <- tmplist[2*no-1][[1]]
      predictVarianceVec[node == unique(node)[no]] <- tmplist[2*no][[1]]
    }
    stopCluster(cl)
  }else{
    for(nnn in unique(node)){
      # predict value depending whether the leaf is pure
      if(pure[node == nnn][1]){
        GP.predict.model <- leaf.GP$GP.model[[which(tree$classTable$node[tree$classTable$pure] == nnn)]]
        GP.predict.data <- leaf.GP$GP.data[[which(tree$classTable$node[tree$classTable$pure] == nnn)]]
        predictPoint <- newPoint[which(node == nnn),c(1:contiVarNum, ncol(nodeData)), drop=FALSE]
        predTmp <- GP.prediction(GP.model = GP.predict.model, newPoint = predictPoint)
        predictValueVec[node == nnn] <- predTmp$predictValue
        predictVarianceVec[node == nnn] <- predTmp$predSigmaSq
      }else{
        QQGP.predict.model <- leaf.GP$QQGP.model[[which(unique(tree$classTable$node[!tree$classTable$pure]) == nnn)]]
        QQGP.predict.data <- leaf.GP$QQGP.data[[which(unique(tree$classTable$node[!tree$classTable$pure]) == nnn)]]
        QQGP.predict.design <- leaf.GP$QQGP.design[[which(unique(tree$classTable$node[!tree$classTable$pure]) == nnn)]]
        compareTable <- cbind(QQGP.predict.data[,((contiVarNum+1):(contiVarNum+cateNum)), drop = FALSE], coerceFun.ori(data = QQGP.predict.data, design = QQGP.predict.design, includeSubCate = T)$c.T)
        names(compareTable) <- c(names(compareTable)[1:cateNum], "c.T")
        compareTable <- unique(compareTable)
        suppressMessages(QQGPcomb <- left_join(newPoint[(node == nnn), c((contiVarNum+1):(contiVarNum+cateNum)), drop=F], compareTable)$c.T)
        predictPoint <- cbind(newPoint[which(node == nnn), c(1:contiVarNum)], QQGPcomb, newPoint[which(node == nnn), (contiVarNum+cateNum+1)])
        predTmp <- QQGP.ori.prediction(QQGP.model = QQGP.predict.model, newPoint = predictPoint)
        predictValueVec[node == nnn] <- predTmp$predictValue
        predictVarianceVec[node == nnn] <- predTmp$predSigmaSq
      }
    }
  }
  return(list(predictValueVec = predictValueVec, predictVarianceVec = pmax(0, predictVarianceVec)))
}






################################################################
########## LOOCV summary function (for visualization) ##########
################################################################
LOOCVsummary <- function(tree, data)
{
  # pre-define data informations
  nodeData <- data
  cateNum <- sum(unlist(lapply(nodeData, FUN=is.factor)))
  contiVarNum <- ncol(nodeData)-(cateNum+1)
  contiDim <- contiVarNum
  
  # pre-define original tree informations
  nodePositionVec <- tree$nodePosition
  parentNodePositionVec <- tree$parentNodePosition
  nodeSplitVar <- tree$nodeSplitVar
  nodeMinMaxCorrVec <- tree$nodeMinMaxCorr
  parentSplitLevelVec <- tree$parentSplitLevel
  parentSplitLogicList <- tree$parentSplitLogic
  unionParentSplitLogicList <- tree$unionParentSplitLogic
  classTable <- tree$classTable
  crossCorrMatrix <- tree$tauMatrix
  leafMuVec <- tree$muVec
  commonFi <- tree$fiParams
  unionData <- tree$unionData
  coerceData <- tree$coerceData
  unionDesign <- tree$unionDesign
  coerceUnionData <- tree$coerceUnionData
  initGPRange <- tree$initGPRange
  
  combineRoot <- numeric(0)
  combineLeaves <- numeric(0)
  
  LOOCV.save <- rep(0, max(nodePositionVec))
  for(root in rev(unique(parentNodePositionVec)[-1])){
    
    
    combineLeafNodePosition <- nodePositionVec[parentNodePositionVec == root]
    combineLeftLeafNodePosition <- nodePositionVec[parentNodePositionVec == root][1]
    combineRightLeafNodePosition <- nodePositionVec[parentNodePositionVec == root][2]
    combineCoerceIndex <- classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root]][order(classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root]])]
    combineCoerceIndexLeft <- classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root][1]][order(classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root][1]])]
    combineCoerceIndexRight <- classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root][2]][order(classTable$coerceIndex[classTable$node %in% nodePositionVec[parentNodePositionVec == root][2]])]
    combineCrossCorrMatrix <- crossCorrMatrix[combineCoerceIndex, combineCoerceIndex]
    combineCrossCorrMatrixLeft <- as.matrix(crossCorrMatrix[combineCoerceIndexLeft, combineCoerceIndexLeft])
    combineCrossCorrMatrixRight <- as.matrix(crossCorrMatrix[combineCoerceIndexRight, combineCoerceIndexRight])
    
    if(root == 1){
      LOOCV.save[1] <- LOOCVFun(coerceUnionData[unionParentSplitLogicList[which(nodePositionVec %in% root)][[1]],c(1:contiVarNum, dim(coerceUnionData)[2])],
                                tree$fiParams, leafMuVec[combineCoerceIndex], as.matrix(combineCrossCorrMatrix))
    }
    LOOCV.save[parentNodePositionVec==root][1] <- LOOCVFun(coerceUnionData[unionParentSplitLogicList[nodePositionVec == combineLeftLeafNodePosition][[1]],c(1:contiVarNum, dim(coerceUnionData)[2])],
                                                           tree$fiParams, leafMuVec[combineCoerceIndexLeft], combineCrossCorrMatrixLeft)
    LOOCV.save[parentNodePositionVec==root][2] <- LOOCVFun(coerceUnionData[unionParentSplitLogicList[nodePositionVec == combineRightLeafNodePosition][[1]],c(1:contiVarNum, dim(coerceUnionData)[2])],
                                                           tree$fiParams, leafMuVec[combineCoerceIndexRight], combineCrossCorrMatrixRight)
    
    combineRoot <- c(combineRoot, root)
    combineLeaves <- c(combineLeaves, combineLeafNodePosition)
    
    removePosition <- (parentNodePositionVec %in% combineRoot)
    # prune the tree (specifically, combine these two nodes)
    nodePositionVecTmp <- nodePositionVec[!removePosition]
    parentNodePositionVec <- parentNodePositionVec[!removePosition]
    nodeSplitVar[nodePositionVec %in% combineRoot] <- "leaf"
    nodeSplitVar <- nodeSplitVar[!removePosition]
    nodeMinMaxCorrVec[nodePositionVec %in% combineRoot] <- 2
    nodeMinMaxCorrVec <- nodeMinMaxCorrVec[!removePosition]
    parentSplitLevelVec <- parentSplitLevelVec[!removePosition]
    parentSplitLogicList <- parentSplitLogicList[!removePosition]
    unionParentSplitLogicList <- unionParentSplitLogicList[!removePosition]
    classTable[(classTable$node %in% combineLeaves), "pure"] <- FALSE
    nodePositionVec <- nodePositionVecTmp
    
    for(r in 1:length(combineRoot)){
      classTable[(classTable$node %in% combineLeaves[(2*(r-1)+1):(2*r)]), "node"] <- combineRoot[r]
    }
  }
  return(LOOCV.save)
}



############################################################
########## plot Tree function (with LOOCV values) ##########
############################################################
plotTree <- function(tree, cateVar = NULL){
  posVec <- tree$nodePosition
  depth <- max(trunc(log(posVec)/log(2)))
  yCenter <- 3*max(trunc(log(posVec)/log(2)))-4*trunc(log(posVec)/log(2))+2
  yLow <- yCenter-0.5
  yHigh <- yCenter+0.5
  xCenter <- ((4*(2^max(trunc(log(posVec)/log(2)))) / (2^trunc(log(posVec)/log(2)))) * (0.5+(posVec-2^(trunc(log(posVec)/log(2))))))
  xLeft <- xCenter-1.5
  xRight <- xCenter+1.5
  parVec <- tree$parentNodePosition
  parVec[1] <- 1
  yParent <- 3*max(trunc(log(posVec)/log(2)))-4*trunc(log(parVec)/log(2))+2
  yParent[1] <- yParent[1] + 3
  xParent <- ((4*(2^max(trunc(log(posVec)/log(2)))) / (2^trunc(log(parVec)/log(2)))) * (0.5+(parVec-2^(trunc(log(parVec)/log(2))))))
  slevel <- gsub("\\.", "=", tree$parentSplitLevel)
  splitVar <- tree$nodeSplitVar
  if(!is.null(cateVar)){
    slevel <- gsub("c", cateVar, slevel)
    splitVar <- gsub("c.", cateVar, tree$nodeSplitVar)
  }
  
  LOOCV.out <- LOOCVsummary(tree,data)
  rVec <- numeric(length(tree$nodeMinMaxCorr))
  for(i in 1:length(tree$nodeMinMaxCorr)){
    if(tree$nodeMinMaxCorr[i] != 2){
      rVec[i] <- paste0(format(round(tree$nodeMinMaxCorr, 2), nsmall = 2)[i], "\n", "(",format(LOOCV.out[i],digits=4),")")
    }else{
      rVec[i] <- paste0(splitVar[i], "\n", "(", format(LOOCV.out[i],digits=4), ")")
    }
  }
  d <- data.frame(x1=xLeft, x2=xRight, y1=yLow, y2=yHigh, r=rVec)
  df <- data.frame(xCenter = xCenter, xParent = xParent, yCenter = yCenter, yParent = yParent, slevel = slevel)
  returnPlot <- ggplot() + 
    geom_segment(aes(x = xCenter, y = yCenter, xend = xCenter, yend = yParent), data = df, color="black") +
    geom_segment(aes(x = xCenter, y = yParent, xend = xParent, yend = yParent), data = df, color="black") +
    geom_rect(data=d, mapping=aes(xmin=xCenter-1.5, xmax=xCenter+1.5, ymin=(yCenter+yParent)/2-0.3, ymax=(yCenter+yParent)/2+0.3),
              color="white", fill="white", alpha=1) + 
    geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", fill="white", alpha=1) +
    geom_text(data=d, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=r), size=4) +
    geom_text(data=df, aes(x=xCenter, y=(yCenter+yParent)/2, label=slevel), size=4) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  return(returnPlot)
}



