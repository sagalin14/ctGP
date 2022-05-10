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


testFunc = "Borehole"
unionOrNot = "common"

num <- 20
type <- unionOrNot



source("TBGP.R")
source("QQGPmodelOrigin5LevelQQ")



t.list <- list(c(2,2,2), c(2,2,3), c(2,3,3), c(3,3,3), c(3,3,4), c(3,4,4), c(4,4,4), c(4,4,5), c(4,5,5), c(5,5,5))

for(tl in 1:length(t.list)){



t1LevelNum <- t.list[[tl]][1]
t2LevelNum <- t.list[[tl]][2]
t3LevelNum <- t.list[[tl]][3]


print(paste0(num, " pts design", ", Level ", t1LevelNum, " ", t2LevelNum, " ", t3LevelNum, " type ", type, " has begun."))

iternum = 5

GPR = c(rep(-4.5, 5), rep(-1, 5)) #c(rep(-4, 6), rep(-0.5, 6))
QQGPR = c(2*pi/5, 3*pi/5)


responseBorehole <- function(xx, t)
{
  rwB <- seq(0.05, 0.15, length.out = t1LevelNum)
  HlB <- seq(700, 820, length.out = t2LevelNum)
  LB <- seq(1120, 1680, length.out = t3LevelNum)
  levelList <- expand.grid(LB,HlB,rwB)[,c(3,2,1)]
  
  rw <- levelList[t,1]
  r  <- xx[1]
  Tu <- xx[2]
  Hu <- xx[3]
  Tl <- xx[4]
  Hl <- levelList[t,2]
  L  <- levelList[t,3]#xx[5]
  Kw <- xx[5]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}




inputScaleFun <- function(xx)
{
  r  <- (50000-100)*xx[1]+100
  Tu <- (115600-63070)*xx[2]+63070
  Hu <- (1110-990)*xx[3]+990
  Tl <- (116-63.1)*xx[4]+63.1
  #L  <- (1680-1120)*xx[5]+1120
  Kw <- (12045-9855)*xx[5]+9855
  
  return(c(r, Tu, Hu, Tl, Kw))
}



levelFun <- function(t)
{
  rwB <- paste0("c1.", c(1:t1LevelNum))
  HlB <- paste0("c2.", c(1:t2LevelNum))
  LB <- paste0("c3.", c(1:t3LevelNum))
  levelList <- expand.grid(LB,HlB,rwB)[,c(3,2,1)]
  names(levelList) <- c("c.1", "c.2", "c.3")
  
  return(levelList[t,])
}








# testing data
set.seed(999999)
x.test <- maximinLHS(n=1000, k=5, maxIter = 100000)
len.test <- dim(x.test)[1]
catePart.test <- data.frame(z = c(rep(1:(t1LevelNum*t2LevelNum*t3LevelNum), each = len.test)))


yvec.test <- numeric(t1LevelNum*t2LevelNum*t3LevelNum*len.test)
for(tlev in 1:(t1LevelNum*t2LevelNum*t3LevelNum)){
  yvec.test[((tlev-1)*len.test+1):(tlev*len.test)] <- apply(x.test, 1, function(x) return(responseBorehole(xx=inputScaleFun(x), t=tlev)))
}


data.test <- data.frame(x.1 = x.test[,1],
                        x.2 = x.test[,2],
                        x.3 = x.test[,3],
                        x.4 = x.test[,4],
                        x.5 = x.test[,5],
                        #x.6 = x.test[,6],
                        c.1 = do.call(rbind, lapply(rep(1:(t1LevelNum*t2LevelNum*t3LevelNum), each=len.test), levelFun))[,1],
                        c.2 = do.call(rbind, lapply(rep(1:(t1LevelNum*t2LevelNum*t3LevelNum), each=len.test), levelFun))[,2],
                        c.3 = do.call(rbind, lapply(rep(1:(t1LevelNum*t2LevelNum*t3LevelNum), each=len.test), levelFun))[,3],
                        Y = yvec.test)

levelNum <- (t1LevelNum*t2LevelNum*t3LevelNum)
numList <- 20
resultList <- list()


for(num in numList){
  
  pts = num
  
  
  for(comm in c(num)){
    
    if(comm == num){type <- "Comm"}else{type <- "Union"}
    
    commonPts = comm
    indepPts = num - commonPts
    
    TBGPmseRecord <- numeric(iternum)
    TBGPtimeRecord <- numeric(iternum)
    TBGPpredictTimeRecord <- numeric(iternum)
    
    GPmseRecord <- numeric(iternum)
    GPtimeRecord <- numeric(iternum)
    GPpredictTimeRecord <- numeric(iternum)
    
    QQGPmseRecord <- numeric(iternum)
    QQGPtimeRecord <- numeric(iternum)
    QQGPpredictTimeRecord <- numeric(iternum)
    
    tgpmseRecord <- numeric(iternum)
    tgptimeRecord <- numeric(iternum)
    
    LVGPmseRecord <- numeric(iternum)
    LVGPtimeRecord <- numeric(iternum)
    
    
    
    set.seed(99999)
    if(commonPts == 0){
      commonDesignList <- lapply(1:iternum, FUN = function(x) return(NULL))
    }else{
      commonDesignList <- lapply(1:iternum, FUN = function(x) return(maximinLHS(commonPts, 5, maxIter = 10000000)))
    }
    
    if(indepPts != 0){
      indepDesignList <- lapply(1:(levelNum*iternum), FUN = function(x) return(maximinLHS(indepPts, 5, maxIter = 1000000)))
    }else{
      indepDesignList <- lapply(1:(levelNum*iternum), FUN = function(x) return(NULL))
    }
    
    
    
    for(iter in 1:iternum)
    {
      x.tr <- NULL
      for(case in 1:levelNum){
        x.tmp <- rbind(commonDesignList[[iter]], indepDesignList[[levelNum*(iter-1)+case]])
        x.tr <- rbind(x.tr, x.tmp)
      }
      len <- dim(x.tr)[1]/levelNum
      
      catePart <- data.frame(z = c(rep(1:(t1LevelNum*t2LevelNum*t3LevelNum), each = len)))
      
      yvec <- numeric(t1LevelNum*t2LevelNum*t3LevelNum*len)
      for(tlev in 1:(t1LevelNum*t2LevelNum*t3LevelNum)){
        yvec[((tlev-1)*len+1):(tlev*len)] <- apply(x.tr[((tlev-1)*len+1):(tlev*len),], 1, function(x) return(responseBorehole(xx=inputScaleFun(x), t=tlev)))
      }
      
      data <- data.frame(x.1 = x.tr[,1],
                         x.2 = x.tr[,2],
                         x.3 = x.tr[,3],
                         x.4 = x.tr[,4],
                         x.5 = x.tr[,5],
                         #x.6 = x.tr[,6],
                         c.1 = do.call(rbind, lapply(rep(1:(t1LevelNum*t2LevelNum*t3LevelNum), each=len), levelFun))[,1],
                         c.2 = do.call(rbind, lapply(rep(1:(t1LevelNum*t2LevelNum*t3LevelNum), each=len), levelFun))[,2],
                         c.3 = do.call(rbind, lapply(rep(1:(t1LevelNum*t2LevelNum*t3LevelNum), each=len), levelFun))[,3],
                         Y = yvec)
      
      #################################################
      # TGP tree / GP fit
      print(paste0("TBGP has begun."))
      TBGPTime1 <- proc.time()
      set.seed(99999)
      tree1 <- treeFit.corrMat(data = data, initRange = GPR,
                               BCDloopNum = 3, mcmcIterNum = 50, nugget = sqrt(.Machine$double.eps))
      treeFinal <- pruneFun(tree1, data)
      leaf.GP1 <- leafGPFun(data = data, GPRange = GPR, QQfiRange = GPR, QQthetaRange = QQGPR,
                            tree = treeFinal, NRF = F)
      TBGPTime2 <- proc.time()
      TBGPTime2 - TBGPTime1

      ###################
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
      #################################################
      # original QQGP fit
      QQGPdata <- coerceFun.ori(data = data, design = data[,1:5])
      print(paste0("QQGP has begun."))
      QQGPTime1 <- proc.time()
      set.seed(99999)
      QQGPfit <- QQGP.ori.model.QQ5level(design = QQGPdata[,1:5], data = QQGPdata, initParam = rep(-1,5))
      QQGPTime2 <- proc.time()
      QQGPTime2 - TBGPTime1
      ###################
      QQGPpredictTime1 <- proc.time()
      diff <- QQGP.ori.prediction(QQGP.model = QQGPfit, newPoint = coerceFun.ori(data = data.test, design = data.test[,1:5]))$predictValue - data.test[,dim(data.test)[2]]
      mse <- sum(diff^2)/dim(data.test)[1]
      QQGPpredictTime2 <- proc.time()
      QQGPperform <- c(mse = sqrt(mse), CPUtime = (QQGPTime2 - QQGPTime1)[3], PredictTime = (QQGPpredictTime2 - QQGPpredictTime1)[3])


      QQGPmseRecord[iter] <- QQGPperform[1]
      QQGPtimeRecord[iter] <- QQGPperform[2]
      QQGPpredictTimeRecord[iter] <- QQGPperform[3]

      print(paste0("RMSE: ", round(sqrt(mse), 2), ", time: ", round((QQGPTime2 - QQGPTime1)[3]), ", predtict Time: ", round((QQGPpredictTime2 - QQGPpredictTime1)[3])))
      print(paste0("QQGP has been completed."))

      #################################################
      # original GP fit
      print(paste0("GP has begun."))
      GPmseRecordSingle <- numeric(1000*(t1LevelNum*t2LevelNum))

      GPTime1 <- proc.time()
      for(t in 1:(t1LevelNum*t2LevelNum*t3LevelNum)){
        assign(paste0("GPdata", t), data[(t-1)*20+(1:20), c(1:5,9)])
        assign(paste0("GPmodel",t), GP.model(get(paste0("GPdata", t))))
      }
      GPTime2 <- proc.time()

      GPpredictTime1 <- proc.time()
      for(t in 1:(t1LevelNum*t2LevelNum*t3LevelNum)){
        assign(paste0("GPpredict", t), GP.prediction(get(paste0("GPmodel",t)), newPoint = data.test[(t-1)*1000+(1:1000), c(1:5,9)]))
        GPmseRecordSingle[(t-1)*1000+(1:1000)] <- get(paste0("GPpredict", t))$predictValue - data.test[(t-1)*1000+(1:1000), c(9)]
      }
      GPpredictTime2 <- proc.time()

      mse <- sum(GPmseRecordSingle^2)/dim(data.test)[1]

      GPperform <- c(mse = sqrt(mse), CPUtime = (GPTime2 - GPTime1)[3], PredictTime = (GPpredictTime2 - GPpredictTime1)[3])

      GPmseRecord[iter] <- GPperform[1]
      GPtimeRecord[iter] <- GPperform[2]
      GPpredictTimeRecord[iter] <- GPperform[3]

      print(paste0("RMSE: ", round(sqrt(mse), 2), ", time: ", round((GPTime2 - GPTime1)[3]), ", predtict Time: ", round((GPpredictTime2 - GPpredictTime1)[3])))
      print(paste0("GP has been completed."))
      
      
      
      ############################################################################
      # tgp training data setting
      c.1.level <- paste0("c1.", 1:t1LevelNum)
      c.2.level <- paste0("c2.", 1:t2LevelNum)
      c.3.level <- paste0("c3.", 1:t3LevelNum)

      Imat <- matrix(0, ncol = t1LevelNum*t2LevelNum*t3LevelNum-1, nrow = dim(data)[1])
      for(i in 1:dim(data)[1]){
        #print(which(data[i,]$c.1 == c.1.level) + t1LevelNum * (which(data[i,]$c.2 == c.2.level)-1) + t1LevelNum * t2LevelNum * (which(data[i,]$c.3 == c.3.level)-1) -1)
        Imat[i,(which(data[i,]$c.1 == c.1.level) + t1LevelNum * (which(data[i,]$c.2 == c.2.level)-1) + t1LevelNum * t2LevelNum * (which(data[i,]$c.3 == c.3.level)-1) -1)] <- 1
      }
      Xdat <- cbind(data$x.1, data$x.2, data$x.3, data$x.4, data$x.5, Imat)
      Y <- data$Y
      ##############################################################################
      # tgp testing data setting
      pred.Imat <- matrix(0, ncol = t1LevelNum*t2LevelNum*t3LevelNum-1, nrow = dim(data.test)[1])
      for(i in 1:dim(data.test)[1]){
        pred.Imat[i,(which(data.test[i,]$c.1 == c.1.level) + t1LevelNum * (which(data.test[i,]$c.2 == c.2.level)-1) + t1LevelNum * t2LevelNum * (which(data.test[i,]$c.3 == c.3.level)-1)-1)] <- 1
      }
      pred.Xdat <- cbind(data.test$x.1, data.test$x.2, data.test$x.3, data.test$x.4, data.test$x.5, pred.Imat)
      Y.test <- data.test$Y

      print(paste0("tGP has begun."))
      tgpTime1 <- proc.time()
      set.seed(99999)
      fit1 <- btgp(X = Xdat, XX= pred.Xdat, Z = Y, verb = 0, meanfn = "constant", basemax = 4+t1LevelNum*t2LevelNum*t3LevelNum,
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
      #############################################################################
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
      print(paste0("Iter ", iter, " has been done."))
      cat("\n")
      
      MSE.result <- data.frame(GP = GPmseRecord, TBGP = TBGPmseRecord, QQGP = QQGPmseRecord, tGP = tgpmseRecord, LVGP = LVGPmseRecord)
      fitTime.result <- data.frame(GP = GPtimeRecord, TBGP = TBGPtimeRecord, QQGP = QQGPtimeRecord, tGP = tgptimeRecord, LVGP = LVGPtimeRecord)
      predTime.result <- data.frame(GP = GPpredictTimeRecord, TBGP = TBGPpredictTimeRecord, QQGP = QQGPpredictTimeRecord)
      
      write.csv(MSE.result, file = paste0(pathlink, "/results/3q/", testFunc, numList[1], "Level", t1LevelNum, t2LevelNum, t3LevelNum, "RMSE", type, ".csv"))
      write.csv(fitTime.result, file = paste0(pathlink, "/results/3q/", testFunc, numList[1], "Level", t1LevelNum, t2LevelNum, t3LevelNum, "FitTime", type, ".csv"))
      write.csv(predTime.result, file = paste0(pathlink, "/results/3q/", testFunc, numList[1], "Level", t1LevelNum, t2LevelNum, t3LevelNum, "PredTime", type, ".csv"))
      
      
    }
    
    print(paste0("GP: ", round(mean(GPmseRecord), 2), ",  TBGP: ", round(mean(TBGPmseRecord), 2), ",  QQGP: ", round(mean(QQGPmseRecord), 2), ",  tgp: ", round(mean(tgpmseRecord), 2), ",  LVGP: ", round(mean(LVGPmseRecord), 2)))
    
    
    print(paste0(num, " pts design", ", Level ", t1LevelNum, " ", t2LevelNum, " ", t3LevelNum, " type ", type, " has been done."))
    cat("\n")
  }
}



}


