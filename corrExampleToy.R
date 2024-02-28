library(tgp)
library(DiceKriging)
library(SLHD)
library(LVGP)
library(lhs)

if(.Platform$OS.type == "windows"){
  pathlink = "D:/Dropbox/RB/TGP/"
} else{
  pathlink = "/Users/SagaLin/Dropbox/RB/TGP"
}

setwd(pathlink)



# source("GPmodel5L.R")
# source("globalGPinit.R")
# source("QQGPmodelOrigin5LNRF.R")
# source("coerceFun.R")
# source("TreeFitModelpartCorrVer.R")
# source("plotTree.R")
# source("classFun.R")
# source("leafGP5L.R")
# source("update_functions_0605.R")

source("QQGPmodelOrigin5LevelQQ.R")
source("TBGP.R")
GPR = c(-1, 1)
QQGPR = c(0, pi/2)






response <- function(x, z1, z2){
  if((z1 == 1)&(z2 == 1)){
    # function 1
    tmp <- cos(3.5*pi*x)*exp(-2*x)#sin(3*x)*exp(-x)
    return(tmp)
  }else if((z1 == 1)&(z2 == 2)){
    # function 2
    tmp <- cos(3.5*pi*x)*exp(-0.5*x)#sin(3*x)
    return(tmp)
  }else if((z1 == 2)&(z2 == 1)){
    # function 3
    tmp <- x-0.5
    return(tmp)
  }else if((z1 == 2)&(z2 == 2)){
    # function 4
    tmp <- 2*x-1
    return(tmp)
  }
}


responseLVGPmath1 <- function(x1, x2, z){
  if(z == 1){
    # function 1
    tmp <- 7 * sin(2 * pi * x1 - pi) + sin(2 * pi * x2 - pi)
    return(tmp)
  }else if(z == 2){
    # function 2
    tmp <- 7 * sin(2 * pi * x1 - pi) + 13 * sin(2 * pi * x2 - pi)
    return(tmp)
  }else if(z == 3){
    # function 3
    tmp <- 7 * sin(2 * pi * x1 - pi) + 1.5 * sin(2 * pi * x2 - pi)
    return(tmp)
  }else if(z == 4){
    # function 4
    tmp <- 7 * sin(2 * pi * x1 - pi) + 9 * sin(2 * pi * x2 - pi)
    return(tmp)
  }else if(z == 5){
    # function 5
    tmp <- 7 * sin(2 * pi * x1 - pi) + 4.5 * sin(2 * pi * x2 - pi)
    return(tmp)
  }
}






# testing data
x.test.vec <- seq(0, 1, by = 0.01)
len.test <- length(x.test.vec)


yvec.test <- numeric(4*len.test)
yvec.test[1:len.test] <- sapply(x.test.vec, function(x) return(response(x, 1, 1)))
yvec.test[(len.test+1):(2*len.test)] <- sapply(x.test.vec, function(x) return(response(x, 1, 2)))
yvec.test[(2*len.test+1):(3*len.test)] <- sapply(x.test.vec, function(x) return(response(x, 2, 1)))
yvec.test[(3*len.test+1):(4*len.test)] <- sapply(x.test.vec, function(x) return(response(x, 2, 2)))

data.test <- data.frame(x.1 = rep(x.test.vec, 4),
                        c.1 = as.factor(c(rep(paste0("c1.",1:2), each = 2*len.test))),
                        c.2 = as.factor(rep(c(rep(c("c2.1", "c2.2"), each = len.test)), 2)),
                        Y = yvec.test)


numList <- 6
resultList <- list()


for(num in numList){
  
  pts = num
  iternum = 100
  
  TBGPmseRecord <- numeric(iternum)
  TBGPtimeRecord <- numeric(iternum)
  TBGPpredictTimeRecord <- numeric(iternum)
  
  tgpmseRecord <- numeric(iternum)
  tgptimeRecord <- numeric(iternum)
  
  LVGPmseRecord <- numeric(iternum)
  LVGPtimeRecord <- numeric(iternum)
  
  QQGPmseRecord <- numeric(iternum)
  QQGPtimeRecord <- numeric(iternum)
  QQGPpredictTimeRecord <- numeric(iternum)
  
  GPmseRecord <- numeric(iternum)
  GPtimeRecord <- numeric(iternum)
  GPpredictTimeRecord <- numeric(iternum)
  
  
  set.seed(10000)
  designList <- lapply(1:iternum, FUN = function(x) return(maximinLHS(pts, 1, maxIter = 1000000)))
  
  for(iter in 1:iternum)
  {
    x.tr <- designList[[iter]]
    len <- dim(x.tr)[1]
    
    yvec <- numeric(4*len)
    yvec[1:len] <- sapply(x.tr, function(x) return(response(x, 1, 1)))
    yvec[(len+1):(2*len)] <- sapply(x.tr, function(x) return(response(x, 1, 2)))
    yvec[(2*len+1):(3*len)] <- sapply(x.tr, function(x) return(response(x, 2, 1)))
    yvec[(3*len+1):(4*len)] <- sapply(x.tr, function(x) return(response(x, 2, 2)))
    data <- data.frame(x.1 = rep(x.tr[,1], 4),
                       c.1 = as.factor(c(rep(paste0("c1.",1:2), each = 2*len))),
                       c.2 = as.factor(rep(c(rep(c("c2.1", "c2.2"), each = len)), 2)),
                       Y = yvec)
    
    #################################################
    # TGP tree / GP fit
    TBGPTime1 <- proc.time()
    tree1 <- treeFit.corrMat(data = data,  nugget = sqrt(.Machine$double.eps))
    treeFinal <- pruneFun(tree1, data)
    plotTree(treeFinal)
    leaf.GP1 <- leafGPFun(data = data, tree = treeFinal, GPRange = GPR, QQfiRange = GPR, QQthetaRange = QQGPR)
    TBGPTime2 <- proc.time()
    TBGPTime2 - TBGPTime1
    
    assign(paste0("plotNum", num, "Iter", iter, "tree"), treeFinal)
    assign(paste0("plotNum", num, "Iter", iter, "GP"), leaf.GP1)
    ###################
    plotTree(treeFinal)
    
    print(paste0("TBGP has been built, starting predicting."))
    
    TBGPpredictTime1 <- proc.time()
    #diff <- sapply(1:dim(data.test)[1], FUN = function(x) return(predictFun(newPoint = data.test[x, ], tree = treeFinal, leaf.GP = leaf.GP1) - data.test[x,dim(data.test)[2]]))
    diff <- predictFun(newPoint = data.test, treeFinal, leaf.GP1)$predictValueVec - data.test[,dim(data.test)[2]]
    mse <- sum(diff^2)/dim(data.test)[1]
    TBGPpredictTime2 <- proc.time()
    TBGPperform <- c(mse = sqrt(mse), CPUtime = (TBGPTime2 - TBGPTime1)[3], PredictTime = (TBGPpredictTime2 - TBGPpredictTime1)[3])
    
    
    TBGPmseRecord[iter] <- TBGPperform[1]
    TBGPtimeRecord[iter] <- TBGPperform[2]
    TBGPpredictTimeRecord[iter] <- TBGPperform[3]
    
    print(paste0("TBGP has been completed."))
    
    #################################################
    # original QQGP fit
    QQGPdata <- coerceFun.ori(data = data, design = data[,1, drop=F])
    #print(paste0("QQGP has begun."))
    QQGPTime1 <- proc.time()
    set.seed(99999)
    QQGPfit <- QQGP.ori.model.QQ5level(design = QQGPdata[,1, drop=F], data = QQGPdata, initParam = c(0))
    QQGPTime2 <- proc.time()
    QQGPTime2 - TBGPTime1
    ###################
    QQGPpredictTime1 <- proc.time()
    diff <- QQGP.ori.prediction(QQGP.model = QQGPfit, newPoint = coerceFun.ori(data = data.test, design = data.test[,1, drop=F]))$predictValue - data.test[,dim(data.test)[2]]
    mse <- sum(diff^2)/dim(data.test)[1]
    QQGPpredictTime2 <- proc.time()
    QQGPperform <- c(mse = sqrt(mse), CPUtime = (QQGPTime2 - QQGPTime1)[3], PredictTime = (QQGPpredictTime2 - QQGPpredictTime1)[3])


    QQGPmseRecord[iter] <- QQGPperform[1]
    QQGPtimeRecord[iter] <- QQGPperform[2]
    QQGPpredictTimeRecord[iter] <- QQGPperform[3]

    #print(paste0("RMSE: ", round(sqrt(mse), 2), ", time: ", round((QQGPTime2 - QQGPTime1)[3]), ", predtict Time: ", round((QQGPpredictTime2 - QQGPpredictTime1)[3])))
    print(paste0("QQGP has been completed."))
    
    #################################################
    # original GP fit
    #print(paste0("GP has begun."))
    GPmseRecordSingle <- numeric(101*4)
    
    GPTime1 <- proc.time()
    for(t in 1:4){
      assign(paste0("GPdata", t), data[(t-1)*6+(1:6), c(1,4)])
      assign(paste0("GPmodel",t), GP.model(get(paste0("GPdata", t))))
    }
    GPTime2 <- proc.time()
    
    GPpredictTime1 <- proc.time()
    for(t in 1:4){
      assign(paste0("GPpredict", t), GP.prediction(get(paste0("GPmodel",t)), newPoint = data.test[(t-1)*101+(1:101), c(1,4)]))
      GPmseRecordSingle[(t-1)*101+(1:101)] <- get(paste0("GPpredict", t))$predictValue - data.test[(t-1)*101+(1:101), c(4)]
    }
    GPpredictTime2 <- proc.time()
    
    mse <- sum(GPmseRecordSingle^2)/dim(data.test)[1]
    
    GPperform <- c(mse = sqrt(mse), CPUtime = (GPTime2 - GPTime1)[3], PredictTime = (GPpredictTime2 - GPpredictTime1)[3])
    
    GPmseRecord[iter] <- GPperform[1]
    GPtimeRecord[iter] <- GPperform[2]
    GPpredictTimeRecord[iter] <- GPperform[3]
    
    #print(paste0("RMSE: ", round(sqrt(mse), 2), ", time: ", round((GPTime2 - GPTime1)[3]), ", predtict Time: ", round((GPpredictTime2 - GPpredictTime1)[3])))
    print(paste0("GP has been completed."))
    
    ##############################################################################
    # tgp training data setting
    c.1.level <- paste0("c1.", 1:2)
    c.2.level <- paste0("c2.", 1:2)
    Imat <- matrix(0, ncol = 3, nrow = dim(data)[1])
    for(i in 1:dim(data)[1]){
      Imat[i,2*(which(data[i,]$c.1 == c.1.level)-1)+which(data[i,]$c.2 == c.2.level)-1] <- 1
    }
    Xdat <- cbind(data$x.1, Imat)
    Y <- data$Y
    ##############################################################################
    # tgp testing data setting
    pred.Imat <- matrix(0, ncol = 3, nrow = dim(data.test)[1])
    for(i in 1:dim(data.test)[1]){
      pred.Imat[i,2*(which(data.test[i,]$c.1 == c.1.level)-1)+which(data.test[i,]$c.2 == c.2.level)-1] <- 1
    }
    pred.Xdat <- cbind(data.test$x.1, pred.Imat)
    Y.test <- data.test$Y

    tgpTime1 <- proc.time()
    fit1 <- btgp(X = Xdat, XX= pred.Xdat, Z = Y, verb = 0, meanfn = "constant", basemax = 4)
    tgpTime2 <- proc.time()
    #tgp.trees(fit1)
    tgpmse <- sum((fit1$ZZ.mean - data.test[, dim(data.test)[2]])^2)/dim(pred.Imat)[1]
    # tgp mse
    tgpperform <- c(mse=sqrt(tgpmse), CPUtime = (tgpTime2 - tgpTime1)[3])

    tgpmseRecord[iter] <- tgpperform[1]
    tgptimeRecord[iter] <- tgpperform[2]

    print(paste0("tgp has been completed."))
    ##############################################################################
    # LVGP training data setting
    X_tr <- matrix(0, ncol = 2, nrow = dim(Xdat)[1])
    X_tr[,1:dim(x.tr)[2]] <- Xdat[,1:dim(x.tr)[2]]
    for(row in (len+1):dim(Xdat)[1]){
      X_tr[row, (dim(x.tr)[2]+1)] <- which(Xdat[row, ((dim(x.tr)[2]+1):(dim(x.tr)[2]+dim(Imat)[2]))] == 1) + 1
    }
    X_tr[1:len, (dim(x.tr)[2]+1)] <- 1
    Y_tr <- Y
    ##############################################################################
    # LVGP testing data setting
    X_te <- matrix(0, ncol = 2, nrow = dim(pred.Xdat)[1])
    X_te[,1] <- pred.Xdat[,1]
    for(row in (len.test+1):dim(pred.Xdat)[1]){
      X_te[row, 2] <- which(pred.Xdat[row, ((1+1):(1+dim(pred.Imat)[2]))] == 1) + 1
    }
    X_te[1:len.test, (1+1)] <- 1
    Y_te <- Y.test
    n_te <- nrow(X_te)

    LVGPTime1 <- proc.time()
    model <- LVGP_fit(X_tr, Y_tr, ind_qual = c(2))
    LVGPTime2 <- proc.time()
    LVGPTime2 - LVGPTime1

    output <- LVGP_predict(X_te, model)
    Y_hat <- output$Y_hat

    LVGPmse <- sum((Y_hat-Y_te)^2)/n_te
    LVGPperform <- c(mse=sqrt(LVGPmse), CPUtime = (LVGPTime2 - LVGPTime1)[3])

    LVGPmseRecord[iter] <- LVGPperform[1]
    LVGPtimeRecord[iter] <- LVGPperform[2]

    print(paste0("LVGP has been completed."))
    print(paste0("Iter ", iter, " has been done."))
    cat("\n")
  }
  
  
  resultList <- c(resultList, list(list( ctGP = c(MSE.mean = mean(TBGPmseRecord), MSE.SD = sd(TBGPmseRecord), TIME.mean = mean(TBGPtimeRecord), TIME.SD = sd(TBGPtimeRecord), predictTIME.mean = mean(TBGPpredictTimeRecord), predictTIME.SD = sd(TBGPpredictTimeRecord)),
                                         QQGP = c(MSE.mean = mean(QQGPmseRecord), MSE.SD = sd(QQGPmseRecord), TIME.mean = mean(QQGPtimeRecord), TIME.SD = sd(QQGPtimeRecord), predictTIME.mean = mean(QQGPpredictTimeRecord), predictTIME.SD = sd(QQGPpredictTimeRecord)),
                                         iGP = c(MSE.mean = mean(GPmseRecord), MSE.SD = sd(GPmseRecord), TIME.mean = mean(GPtimeRecord), TIME.SD = sd(GPtimeRecord), predictTIME.mean = mean(GPpredictTimeRecord), predictTIME.SD = sd(GPpredictTimeRecord)),
                                         tgp = c(MSE.mean = mean(tgpmseRecord), MSE.SD = sd(tgpmseRecord), TIME.mean = mean(tgptimeRecord), TIME.SD = sd(tgptimeRecord)),
                                         LVGP = c(MSE.mean = mean(LVGPmseRecord), MSE.SD = sd(LVGPmseRecord), TIME.mean = mean(LVGPtimeRecord), TIME.SD = sd(LVGPtimeRecord)),
                                         is.ctGP.best = ((mean(TBGPmseRecord) < mean(tgpmseRecord)) & (mean(TBGPmseRecord) < mean(LVGPmseRecord))
                                                       & (mean(TBGPmseRecord) < mean(GPmseRecord)) & (mean(TBGPmseRecord) < mean(QQGPmseRecord))))))
  print(paste0(num, " pts design has been done."))
  cat("\n")
}


names(resultList) <- paste0("designPts", numList)

resultList





rmse = data.frame(rmse = c(LVGPmseRecord, tgpmseRecord, QQGPmseRecord, TBGPmseRecord, GPmseRecord),
                  model = factor(rep(rep(c("LVGP", "tgp", "QQGP", "ctGP", "iGP"), each = 100),5), level=c("LVGP", "tgp", "QQGP", "ctGP", "iGP")))


library("ggplot2")

ggplot(aes(y=model,x=rmse),data=rmse) + geom_boxplot()




