library(tgp)
library(DiceKriging)
library(SLHD)
library(LVGP)
library(lhs)


if(.Platform$OS.type == "windows"){
  pathlink = "C:/Users/SagaLin/Desktop/TBGP"
} else{
  pathlink = "/Users/SagaLin/Dropbox/RB/TGP"
}

setwd(pathlink)


source("TBGP.R")


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
x.test <- maximinLHS(n=1000, k=2, maxIter = 100000)
len.test <- dim(x.test)[1]
catePart.test <- data.frame(z = c(rep(1:5, each = len.test)))


yvec.test <- numeric(5*len.test)
yvec.test[1:len.test] <- apply(x.test, 1, function(x) return(responseLVGPmath1(x[1], x[2], 1)))
yvec.test[(len.test+1):(2*len.test)] <- apply(x.test, 1, function(x) return(responseLVGPmath1(x[1], x[2], 2)))
yvec.test[(2*len.test+1):(3*len.test)] <- apply(x.test, 1, function(x) return(responseLVGPmath1(x[1], x[2], 3)))
yvec.test[(3*len.test+1):(4*len.test)] <- apply(x.test, 1, function(x) return(responseLVGPmath1(x[1], x[2], 4)))
yvec.test[(4*len.test+1):(5*len.test)] <- apply(x.test, 1, function(x) return(responseLVGPmath1(x[1], x[2], 5)))

data.test <- data.frame(x.1 = rep(x.test[,1], 5), x.2 = rep(x.test[,2], 5),
                        c.1 = c(rep(paste0("z.1",1:5), each = len.test)),
                        Y = yvec.test)

data.test$c.1 <- as.factor(data.test$c.1)

iternum = 100
testFunc <- "5L"

levelNum <- 5
numList <- 13
resultList <- list()


GPR = c(-2, -2, 0, 0)
GPRi = c(-1.4, -1.4, -0.7, -0.7)
GGPR = c(-2, -2, 3, 3)
QQGPR = c(0.15, 0.5)

for(num in numList){
  
  pts = num
  
  
  for(comm in num){
    
    if(comm == num){type <- "Common"}else{type <- "Union"}
    
    commonPts = comm
    indepPts = num - commonPts
    
    TBGPmseRecord <- numeric(iternum)
    TBGPtimeRecord <- numeric(iternum)
    TBGPpredictTimeRecord <- numeric(iternum)
    
    QQGPmseRecord <- numeric(iternum)
    QQGPtimeRecord <- numeric(iternum)
    QQGPpredictTimeRecord <- numeric(iternum)
    
    tgpmseRecord <- numeric(iternum)
    tgptimeRecord <- numeric(iternum)
    
    LVGPmseRecord <- numeric(iternum)
    LVGPtimeRecord <- numeric(iternum)
    
    GPmseRecord <- numeric(iternum)
    GPtimeRecord <- numeric(iternum)
    GPpredictTimeRecord <- numeric(iternum)
    
    set.seed(999999)
    if(commonPts == 0){
      commonDesignList <- lapply(1:iternum, FUN = function(x) return(NULL))
    }else{

      commonDesignList <- lapply(1:iternum, FUN = function(x) return(maximinLHS(commonPts, 2, maxIter = 10000)))
    }
    
    if(indepPts != 0){
      indepDesignList <- lapply(1:(levelNum*iternum), FUN = function(x) return(maximinLHS(indepPts, 2, maxIter = 1000000)))
    }else{
      indepDesignList <- lapply(1:(levelNum*iternum), FUN = function(x) return(NULL))
    }
    
    
    
    for(iter in 1:iternum)
    {
      #x.tr <- designList[[iter]]
      #len <- dim(x.tr)[1]
      
      x.tr <- NULL
      for(case in 1:levelNum){
        x.tmp <- rbind(commonDesignList[[iter]], indepDesignList[[levelNum*(iter-1)+case]])
        x.tr <- rbind(x.tr, x.tmp)
      }
      len <- dim(x.tr)[1]/levelNum
      
      catePart <- data.frame(z = c(rep(1:5, each = len)))
      
      yvec <- numeric(5*len)
      for(tlev in 1:5){
        yvec[((tlev-1)*len+1):(tlev*len)] <- apply(x.tr[((tlev-1)*len+1):(tlev*len),], 1, function(x) return(responseLVGPmath1(x[1], x[2], tlev)))
      }
      
      data <- data.frame(x.1 = x.tr[,1],
                         x.2 = x.tr[,2],
                         c.1 = c(rep(paste0("z.1",1:5), each = len)),
                         Y = yvec)
      
      data$c.1 <- as.factor(data$c.1)
      
      #################################################
      # TGP tree / GP fit
      print(paste0("TBGP has begun."))
      TBGPTime1 <- proc.time()
      #set.seed(99999)
      set.seed(1234)
      tree1 <- treeFit.corrMat(data = data, initRange = GPRi,
                               BCDloopNum = 3, mcmcIterNum = 50, nugget = sqrt(.Machine$double.eps))
      treeFinal <- pruneFun(tree1, data)
      leaf.GP1 <- leafGPFun(data = data, GPRange = GPR, QQfiRange = GPR, QQthetaRange = QQGPR,
                            tree = treeFinal, NRF = F)
      TBGPTime2 <- proc.time()
      TBGPTime2 - TBGPTime1
      
      plotTree(tree1)
      # save 7X6
      plotTree(treeFinal)
      # save 5X4
      
      ###################
      assign(paste0("plotNum", num, "Iter", iter, "treeFirst"), tree1)
      assign(paste0("plotNum", num, "Iter", iter, "tree"), treeFinal)
      assign(paste0("plotNum", num, "Iter", iter, "GP"), leaf.GP1)
      ###################
      plotTree(treeFinal)
      print(paste0("combined any leaf: ", any(!get(paste0("plotNum", num, "Iter", iter, "tree"))$classTable$pure)))
      print(paste0("root QQ: ", length(unique(get(paste0("plotNum", num, "Iter", iter, "tree"))$classTable$node)) == 1))
      
      print(paste0("Building completed."))
      TBGPpredictTime1 <- proc.time()
      diff <- predictFun(newPoint = data.test, treeFinal, leaf.GP1, parallel = F)$predictValueVec - data.test[,dim(data.test)[2]]
      mse <- sum(diff^2)/dim(data.test)[1]
      TBGPpredictTime2 <- proc.time()
      TBGPperform <- c(mse = sqrt(mse), CPUtime = (TBGPTime2 - TBGPTime1)[3], PredictTime = (TBGPpredictTime2 - TBGPpredictTime1)[3])
      
      print(paste0("RMSE: ", round(sqrt(mse), 2), ", time: ", round((TBGPTime2 - TBGPTime1)[3]), ", predtict Time: ", round((TBGPpredictTime2 - TBGPpredictTime1)[3])))
      
      TBGPmseRecord[iter] <- TBGPperform[1]
      TBGPtimeRecord[iter] <- TBGPperform[2]
      TBGPpredictTimeRecord[iter] <- TBGPperform[3]
      
      
      print(paste0("TBGP has been completed."))
      ################################################
      # original QQGP fit
      QQGPdata <- coerceFun.ori(data = data, design = data[,1:2])
      print(paste0("QQGP has begun."))
      QQGPTime1 <- proc.time()
      set.seed(99999)
      QQGPfit <- #QQGP.ori.model.QQ5level(design = QQGPdata[,1:2], data = QQGPdata, initParam = c(0, 0))
                  QQGP.ori.model(design = QQGPdata[,c(1:2)], data = QQGPdata,
                                 initParam = rep(-.5, 2), fiRange = GPR, thetaRange = QQGPR, Tmatrix = NULL,
                                 nugget = sqrt(.Machine$double.eps))
      QQGPTime2 <- proc.time()
      QQGPTime2 - TBGPTime1
      ###################
      QQGPpredictTime1 <- proc.time()
      diff <- QQGP.ori.prediction(QQGP.model = QQGPfit, newPoint = coerceFun.ori(data = data.test, design = data.test[,1:2]))$predictValue - data.test[,dim(data.test)[2]]
      mse <- sum(diff^2)/dim(data.test)[1]
      QQGPpredictTime2 <- proc.time()
      QQGPperform <- c(mse = sqrt(mse), CPUtime = (QQGPTime2 - QQGPTime1)[3], PredictTime = (QQGPpredictTime2 - QQGPpredictTime1)[3])


      QQGPmseRecord[iter] <- QQGPperform[1]
      QQGPtimeRecord[iter] <- QQGPperform[2]
      QQGPpredictTimeRecord[iter] <- QQGPperform[3]
      print(paste0("RMSE: ", round(sqrt(mse), 2), ", time: ", round((QQGPTime2 - QQGPTime1)[3]), ", predtict Time: ", round((QQGPpredictTime2 - QQGPpredictTime1)[3])))
      print(paste0("QQGP has been completed."))
      ##########################################################################
      # tgp training data setting
      c.1.level <- paste0("z.1", 1:5)
      Imat <- matrix(0, ncol = 4, nrow = dim(data)[1])
      for(i in 1:dim(data)[1]){
        Imat[i,(which(data[i,]$c.1 == c.1.level)-1)] <- 1
      }
      Xdat <- cbind(data$x.1, data$x.2, Imat)
      Y <- data$Y
      ##############################################################################
      # tgp testing data setting
      pred.Imat <- matrix(0, ncol = 4, nrow = dim(data.test)[1])
      for(i in 1:dim(data.test)[1]){
        pred.Imat[i,(which(data.test[i,]$c.1 == c.1.level)-1)] <- 1
      }
      pred.Xdat <- cbind(data.test$x.1, data.test$x.2, pred.Imat)
      Y.test <- data.test$Y

      print(paste0("tGP has begun."))
      tgpTime1 <- proc.time()
      set.seed(99999)
      fit1 <- btgp(X = Xdat, XX= pred.Xdat, Z = Y, verb = 0, meanfn = "constant", basemax = 6,
                   nug.p=0, gd=c(1e-8, 1))
      tgpTime2 <- proc.time()
      #tgp.trees(fit1)
      tgpmse <- sum((fit1$ZZ.mean - data.test[, dim(data.test)[2]])^2)/dim(pred.Imat)[1]
      # tgp mse
      tgpperform <- c(mse=sqrt(tgpmse), CPUtime = (tgpTime2 - tgpTime1)[3])
      print(paste0("RMSE: ", round(sqrt(tgpmse), 2), ", time: ", round((tgpTime2 - tgpTime1)[3])))
      tgpmseRecord[iter] <- tgpperform[1]
      tgptimeRecord[iter] <- tgpperform[2]

      print(paste0("tgp has been completed."))
      ##############################################################################
      # LVGP training data setting
      X_tr <- matrix(0, ncol = (dim(x.tr)[2]+1), nrow = dim(Xdat)[1])
      X_tr[,1:dim(x.tr)[2]] <- Xdat[,1:dim(x.tr)[2]]
      for(row in (len+1):dim(Xdat)[1]){
        X_tr[row, (dim(x.tr)[2]+1)] <- which(Xdat[row, ((dim(x.tr)[2]+1):(dim(x.tr)[2]+dim(Imat)[2]))] == 1) + 1
      }
      X_tr[1:len, (dim(x.tr)[2]+1)] <- 1
      Y_tr <- Y
      ##############################################################################
      # LVGP testing data setting
      X_te <- matrix(0, ncol = (dim(x.test)[2]+1), nrow = dim(pred.Xdat)[1])
      X_te[,1:dim(x.test)[2]] <- pred.Xdat[,1:dim(x.test)[2]]
      for(row in (len.test+1):dim(pred.Xdat)[1]){
        X_te[row, (dim(x.test)[2]+1)] <- which(pred.Xdat[row, ((dim(x.test)[2]+1):(dim(x.test)[2]+dim(pred.Imat)[2]))] == 1) + 1
      }
      X_te[1:len.test, (dim(x.test)[2]+1)] <- 1
      Y_te <- Y.test
      n_te <- nrow(X_te)

      print(paste0("LVGP has begun."))
      LVGPTime1 <- proc.time()
      model <- LVGP_fit(X_tr, Y_tr, ind_qual = c(dim(x.test)[2]+1), parallel = T)
      LVGPTime2 <- proc.time()
      LVGPTime2 - LVGPTime1

      output <- LVGP_predict(X_te, model)
      Y_hat <- output$Y_hat

      LVGPmse <- sum((Y_hat-Y_te)^2)/n_te
      LVGPperform <- c(mse=sqrt(LVGPmse), CPUtime = (LVGPTime2 - LVGPTime1)[3])
      print(paste0("RMSE: ", round(sqrt(LVGPmse), 2), ", time: ", round((LVGPTime2 - LVGPTime1)[3])))
      LVGPmseRecord[iter] <- LVGPperform[1]
      LVGPtimeRecord[iter] <- LVGPperform[2]

      print(paste0("LVGP has been completed."))

      #################################################
      # original GP fit
      print(paste0("GP has begun."))
      GPmseRecordSingle <- numeric(1000*5)
      
      GPTime1 <- proc.time()
      for(t in 1:5){
        assign(paste0("GPdata", t), data[(t-1)*numList+(1:numList), c(1:2,4)])
        assign(paste0("GPmodel",t), GP.model(get(paste0("GPdata", t)), range = GGPR,
                                             initParam = c(3,3), nugget = sqrt(.Machine$double.eps)))
      }
      GPTime2 <- proc.time()
      
      GPpredictTime1 <- proc.time()
      for(t in 1:5){
        assign(paste0("GPpredict", t), GP.prediction(get(paste0("GPmodel",t)), newPoint = data.test[(t-1)*1000+(1:1000), c(1:2,4)]))
        GPmseRecordSingle[(t-1)*1000+(1:1000)] <- get(paste0("GPpredict", t))$predictValue - data.test[(t-1)*1000+(1:1000), c(4)]
      }
      GPpredictTime2 <- proc.time()
      
      mse <- sum(GPmseRecordSingle^2)/dim(data.test)[1]
      
      GPperform <- c(mse = sqrt(mse), CPUtime = (GPTime2 - GPTime1)[3], PredictTime = (GPpredictTime2 - GPpredictTime1)[3])
      
      GPmseRecord[iter] <- GPperform[1]
      GPtimeRecord[iter] <- GPperform[2]
      GPpredictTimeRecord[iter] <- GPperform[3]
      
      print(paste0("RMSE: ", round(sqrt(mse), 2), ", time: ", round((GPTime2 - GPTime1)[3]), ", predtict Time: ", round((GPpredictTime2 - GPpredictTime1)[3])))
      print(paste0("GP has been completed."))
      
      print(paste0("Iter ", iter, " has been done."))
      cat("\n")
      
    }
    
    
    resultList <- c(resultList, list(list( TBGP = c(MSE.mean = mean(TBGPmseRecord), MSE.SD = sd(TBGPmseRecord), TIME.mean = mean(TBGPtimeRecord), TIME.SD = sd(TBGPtimeRecord), predictTIME.mean = mean(TBGPpredictTimeRecord), predictTIME.SD = sd(TBGPpredictTimeRecord)),
                                           QQGP = c(MSE.mean = mean(QQGPmseRecord), MSE.SD = sd(QQGPmseRecord), TIME.mean = mean(QQGPtimeRecord), TIME.SD = sd(QQGPtimeRecord), predictTIME.mean = mean(QQGPpredictTimeRecord), predictTIME.SD = sd(QQGPpredictTimeRecord)),
                                           tgp = c(MSE.mean = mean(tgpmseRecord), MSE.SD = sd(tgpmseRecord), TIME.mean = mean(tgptimeRecord), TIME.SD = sd(tgptimeRecord), predictTIME.mean = NA, predictTIME.SD = NA),
                                           LVGP = c(MSE.mean = mean(LVGPmseRecord), MSE.SD = sd(LVGPmseRecord), TIME.mean = mean(LVGPtimeRecord), TIME.SD = sd(LVGPtimeRecord), predictTIME.mean = NA, predictTIME.SD = NA),
                                           GP = c(MSE.mean = mean(GPmseRecord), MSE.SD = sd(GPmseRecord), TIME.mean = mean(GPtimeRecord), TIME.SD = sd(GPtimeRecord), predictTIME.mean = mean(GPpredictTimeRecord), predictTIME.SD = sd(GPpredictTimeRecord)))))
    print(paste0(num, " pts design has been done."))
    cat("\n")
  }
}

MSE.result <- data.frame(TBGP = TBGPmseRecord, QQGP = QQGPmseRecord, tGP = tgpmseRecord, LVGP = LVGPmseRecord, GP = GPmseRecord)
fitTime.result <- data.frame(TBGP = TBGPtimeRecord, QQGP = QQGPtimeRecord, tGP = tgptimeRecord, LVGP = LVGPtimeRecord, GP = GPtimeRecord)
predTime.result <- data.frame(TBGP = TBGPpredictTimeRecord, QQGP = QQGPpredictTimeRecord, GP = GPpredictTimeRecord)

write.csv(MSE.result, file = paste0(pathlink, "/results/", testFunc, numList[1], "MSE.csv"))
write.csv(fitTime.result, file = paste0(pathlink, "/results/", testFunc, numList[1], "FitTime.csv"))
write.csv(predTime.result, file = paste0(pathlink, "/results/", testFunc, numList[1], "PredTime.csv"))

resultList
write.csv(resultList, file = paste0(pathlink, "/results/", testFunc, numList[1], "summary.csv"))

library("reshape2")
require(gridExtra)
m <- ggplot( melt(data.frame(GPmseRecord, TBGPmseRecord, QQGPmseRecord, tgpmseRecord, LVGPmseRecord)), aes(x = variable, y = value)) +
  geom_boxplot(width=0.5)+ labs(x = "model",y="RMSE")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(labels=c("iGP","ctGP", "QQGP","tgp","LVGP"))

t <- ggplot( melt(data.frame(GPtimeRecord+GPpredictTimeRecord, TBGPtimeRecord+TBGPpredictTimeRecord, QQGPtimeRecord+QQGPpredictTimeRecord,tgptimeRecord, LVGPtimeRecord)), aes(x = variable, y = value)) +
  geom_boxplot(width=0.5)+ labs(x = "model",y="Time")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(labels=c("iGP","ctGP", "QQGP","tgp","LVGP"))

grid.arrange(m,t, ncol=2)
