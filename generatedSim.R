library(DiceKriging)
library(SLHD)
library(lhs)
library(caret)
library(GPareto)
library(MASS)
if(.Platform$OS.type == "windows"){
  pathlink = "D:/Dropbox/RB/TGP/"
} else{
  pathlink = "/Users/SagaLin/Dropbox/TBGP"
}

setwd(pathlink)
source("TBGP.R")
library(plgp)
library(mvtnorm)
summaryTable <- list()

##### design function #####
iterNum = 100
ptsNum <- 20
num_trainingSet <- ptsNum/2
contiDim <- 1
X <- matrix(seq(0, 1, length=ptsNum), ncol=1)
all_data <- data.frame(rbind(X,X,X,X))
all_data[,dim(all_data)[2]+1] <- c(rep(1,ptsNum),rep(2,ptsNum),rep(3,ptsNum),rep(4,ptsNum))
levelVec <- unique(all_data[, contiDim+1])
levelNum <- length(levelVec)
levelVecTotal <- all_data[, contiDim+1]
design <- all_data[, 1:contiDim, drop=FALSE]
nugget <- sqrt(.Machine$double.eps) 

GPR = c(rep(0,1),rep(4,1))
QQGPR = c(0,2*pi)
fiSet <- 1.5#seq(1,2,0.1)
tauSet <- c(0)
#Tmatrix <- matrix(c(1.0,0.8,0.6,0.0,
#                    0.8,1.0,0.7,0.0,
#                    0.6,0.7,1.0,0.0,
#                    0.0,0.0,0.0,1.0), ncol=4)
# Tmatrix <- matrix(c(1.0,0.9,0.7,0.3,
#                     0.9,1.0,0.8,0.6,
#                     0.7,0.8,1.0,0.7,
#                     0.3,0.6,0.7,1.0), ncol=4)
# Tmatrix <- matrix(c(1.0,0.8,0.0,0.0,
#                     0.8,1.0,0.0,0.0,
#                     0.0,0.0,1.0,0.9,
#                     0.0,0.0,0.9,1.0), ncol=4)
# Tmatrix <- matrix(c(1.0,0.0,0.0,0.0,
#                     0.0,1.0,0.0,0.0,
#                     0.0,0.0,1.0,0.0,
#                     0.0,0.0,0.0,1.0), ncol=4)


# Tmatrix <- matrix(c(1.0,0.8,0.6,0.4,
#                     0.8,1.0,0.6,0.4,
#                     0.6,0.6,1.0,0.2,
#                     0.4,0.4,0.2,1.0), ncol=4)

Tmatrix <- matrix(c(1.0,0.4,0.0,0.0,
                    0.4,1.0,0.0,0.0,
                    0.0,0.0,1.0,0.6,
                    0.0,0.0,0.6,1.0), ncol=4)

contiPosNum <- ifelse(is.null(dim(all_data)), length(design), dim(all_data)[1])
totalPtsNum <- contiPosNum
store_data <- list()

datY <- matrix(0, ncol = totalPtsNum, nrow = nrow(Tmatrix)*length(fiSet)*iterNum)

  fiIt = 1
  for(fi in fiSet){
    tauIt = 1
    for(tau in tauSet){
      Hmatrix <- matrix(0, ncol = totalPtsNum, nrow = totalPtsNum)
      rcoord <- cbind(rep(seq_len(totalPtsNum-1L), times = rev(seq_len(totalPtsNum-1L))), 
                      unlist(lapply(X=rev(seq_len(totalPtsNum-1L)), FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = totalPtsNum)))
      xdiff <- (design[rcoord[, 2L], ] - design[rcoord[, 1L], ])**2
      Beta <- matrix(fi, nrow = choose(totalPtsNum, 2), ncol = contiDim, byrow = TRUE)
      Htemp <- 10^Beta * xdiff; Htemp <- rowSums(Htemp)
      Hmatrix[rcoord] <- Htemp; Hmatrix <- Hmatrix + t(Hmatrix); Hmatrix <- exp(-Hmatrix)
      Rmatrix <- matrix(0, ncol = totalPtsNum, nrow = totalPtsNum)
      for(row in 1:totalPtsNum){
        for(col in 1:totalPtsNum){
          Rmatrix[row, col] <- Hmatrix[row, col] * Tmatrix[levelVecTotal[row], levelVecTotal[col]]
        }
      }
      if(nugget == 0){
        Rmatrix <- Rmatrix + diag(rep(computeNug(Rmatrix), totalPtsNum))
      }else{
        Rmatrix <- Rmatrix + diag(rep(nugget, totalPtsNum))
      }
      set.seed(999999)
      datY[(iterNum*(tauIt-1)+iterNum*length(tauSet)*(fiIt-1)+1):(iterNum*(tauIt)+iterNum*length(tauSet)*(fiIt-1)), ] <- mvrnorm(iterNum, mu = rep(0, times = totalPtsNum), Sigma = 10*Rmatrix)
  
      tauIt = tauIt + 1
    }
    fiIt = fiIt + 1
  }





for(it in 1:iterNum){
  cat(paste0('Iter ',it, "\n"))
  fiIt = 1
for(fi in fiSet){
  tauIt = 1
  for(tau in tauSet){
    cat(paste0('  fi = ', fi, ", tau = NULL", "..."))
  
  Y <- datY[iterNum*(tauIt-1)+iterNum*length(tauSet)*(fiIt-1)+it,]

  tauIt = tauIt + 1
  
  plot(range(X), range(Y), xlab = "x", ylab = "y", type = "n",
       main = paste0("Data generated via mixed-input GP \ncross corr. = ", tau, " and fi = ", fi))
  lines(all_data[1:ptsNum,contiDim], Y[1:ptsNum], col = 'blue', lwd = 1.5)
  lines(all_data[(ptsNum+1):(2*ptsNum),contiDim], Y[(ptsNum+1):(2*ptsNum)], col = 'red', lwd = 1.5)
  lines(all_data[(2*ptsNum+1):(3*ptsNum),contiDim], Y[(2*ptsNum+1):(3*ptsNum)], col = 'green', lwd = 1.5)
  lines(all_data[(3*ptsNum+1):(4*ptsNum),contiDim], Y[(3*ptsNum+1):(4*ptsNum)], col = 'orange', lwd = 1.5)

  cor(Y[1:ptsNum],Y[(ptsNum+1):(2*ptsNum)])
  all_data[,3] <- Y
  store_data <- append(store_data,list(all_data))
  true_df <- all_data
  true_df[,contiDim+1] <- as.factor(paste0("c1.",true_df[,contiDim+1]))
  colnames(true_df) <- c('x.1','c.1','Y')
  for(num_training in num_trainingSet){
    index <- sample(nrow(true_df)/levelNum, size=num_training)
    data <- true_df[c(index,index+ptsNum,index+2*ptsNum,index+3*ptsNum),]
    data.test <- true_df[-c(index,index+ptsNum,index+2*ptsNum,index+3*ptsNum),]

    #############################
    ######## Fit model ##########
    ##### ctGP #####
    tree <- treeFit.corrMat(data = data, initRange = GPR,
                            BCDloopNum = 50, mcmcIterNum = 100, nugget = sqrt(.Machine$double.eps))
    treeFinal <- pruneFun(tree, data)
    leaf.GP <- leafGPFun(data = data, GPRange = GPR, QQfiRange = GPR, QQthetaRange = QQGPR,
                         tree = treeFinal, NRF = F)
    diff <- predictFun(newPoint = data.test, treeFinal, leaf.GP, parallel = F)$predictValueVec - data.test[,dim(data.test)[2]]
    ctGPrmse <- sqrt(sum(diff^2)/length(diff))
    treeFinal$fiParams
    treeFinal$tauMatrix
    ##### QQGP #####
    QQGPdata <- coerceFun.ori(data = data, design = data[,1,drop=FALSE])
    QQGPfit <- QQGP.ori.model(design = QQGPdata[,c(1)], data = QQGPdata,
                              initParam = c(4), fiRange = GPR, thetaRange = QQGPR, Tmatrix = NULL,
                              nugget = sqrt(.Machine$double.eps))
    pred <- QQGP.ori.prediction(QQGP.model = QQGPfit, newPoint = coerceFun.ori(data = data.test, design = data.test[,1]))
    diff <- pred$predictValue - data.test[,dim(data.test)[2]]
    QQGPrmse <- sqrt(sum(diff^2)/dim(data.test)[1])
    QQGPfit$fi
    QQGPfit$Tmat
    ## indep GP
    GPdat1 <- QQGPdata[1:num_training, c(1,3)]
    GPdat2 <- QQGPdata[(num_training+1):(2*num_training), c(1,3)]
    GPdat3 <- QQGPdata[(2*num_training+1):(3*num_training), c(1,3)]
    GPdat4 <- QQGPdata[(3*num_training+1):(4*num_training), c(1,3)]
    
    GPmodel1 <- GP.model(GPdat1, initParam = c(4), range = GPR, nugget = sqrt(.Machine$double.eps))
    GPmodel2 <- GP.model(GPdat2, initParam = c(4), range = GPR, nugget = sqrt(.Machine$double.eps))
    GPmodel3 <- GP.model(GPdat3, initParam = c(4), range = GPR, nugget = sqrt(.Machine$double.eps))
    GPmodel4 <- GP.model(GPdat4, initParam = c(4), range = GPR, nugget = sqrt(.Machine$double.eps))
    
    predDat1 <- data.test[1:(ptsNum-num_training), c(1,3)]
    predDat2 <- data.test[((ptsNum-num_training)+1):(2*(ptsNum-num_training)), c(1,3)]
    predDat3 <- data.test[(2*(ptsNum-num_training)+1):(3*(ptsNum-num_training)), c(1,3)]
    predDat4 <- data.test[(3*(ptsNum-num_training)+1):(4*(ptsNum-num_training)), c(1,3)]
    
    predVal1 <- GP.prediction(GPmodel1, predDat1)
    predVal2 <- GP.prediction(GPmodel2, predDat2)
    predVal3 <- GP.prediction(GPmodel3, predDat3)
    predVal4 <- GP.prediction(GPmodel4, predDat4)
    
    diff1 <- predDat1$Y - predVal1$predictValue
    diff2 <- predDat2$Y - predVal2$predictValue
    diff3 <- predDat3$Y - predVal3$predictValue
    diff4 <- predDat4$Y - predVal4$predictValue
    
    GPrmse <- sqrt(sum(c(diff1, diff2, diff3, diff4)^2)/dim(data.test)[1])

    summaryTable <- append(summaryTable,list(c(num_training,fi,treeFinal$fiParams,QQGPfit$fi,
                                            0,QQGPfit$Tmat[1,2],treeFinal$tauMatrix[1,2],
                                            cor(data[,3][1:num_training],data[,3][num_training+(1:num_training)]),
                                            QQGPrmse,ctGPrmse,GPrmse)))
    cat(paste0(" Done ! \n"))
    }
  }
  fiIt = fiIt + 1
}
}
result <- data.frame(matrix(unlist(summaryTable),ncol=11,byrow=TRUE))
colnames(result) <- c('train_point','fi','fi_ctGP','fi_QQGP','tau','tmat_QQGP','BCD','Pearson','RMSE_QQ','RMSE_ctGP','RMSE_iGP')
result <- result %>% mutate(across(c(9,10,11), round, 6))


write.csv(result,'result_rep100.csv', row.names = FALSE)

