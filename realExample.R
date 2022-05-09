library(tgp)
library(DiceKriging)
library(SLHD)
library(LVGP)
library(lhs)

library(RcppArmadillo)


# Create Armadillo function
Rcpp::cppFunction(depends = "RcppArmadillo", code = '
                  Rcpp::NumericMatrix solveCpp(SEXP A) {
                  arma::mat A_c = (arma::mat)Rcpp::as<arma::mat>(A);
                  arma::mat B_c = inv_sympd(A_c);
                  return wrap(B_c);
                  }')

if(.Platform$OS.type == "windows"){
  pathlink = "C:/Users/SagaLin/Desktop/TBGP"
} else{
  pathlink = "/Users/SagaLin/Dropbox/RB/TGP"
}

setwd(pathlink)



source("TBGP.R")




c.1.level <- paste0("c1.", 1:4)
c.2.level <- paste0("c2.", 1:4)
c.3.level <- paste0("c3.", 1:2)
c.vector <- expand.grid(c.3.level, c.2.level, c.1.level)[,c(3,2,1)]
names(c.vector) <- c("c.1", "c.2", "c.3")



testData <- read.csv("testData.csv")
testData$BaseMaterial <- as.character(testData$BaseMaterial)
testData$FinMaterial <- as.character(testData$FinMaterial)
testData$FanType <- as.character(testData$FanType)
for(i in 1:4){
  testData[(testData$BaseMaterial ==  unique(testData$BaseMaterial)[i]), 4] <- paste0("c1.", i)
  testData[(testData$FinMaterial ==  unique(testData$FinMaterial)[i]), 5] <- paste0("c2.", i)
}
for(i in 1:2){
  testData[(testData$FanType ==  unique(testData$FanType)[i]), 6] <- paste0("c3.", i)
}
names(testData) <- c("x.1", "x.2", "x.3", "c.1", "c.2", "c.3", "Y")
testData$c.1 <- factor(testData$c.1)
testData$c.2 <- factor(testData$c.2)
testData$c.3 <- factor(testData$c.3)

data.test <- testData






Pup = 0
Plow = -2
GPR = c(rep(Plow, 3), rep(Pup, 3))
GGPR = c(rep(-2, 3), rep(2, 3))
QQGPR = c(pi/4, 3*pi/4)


# Data settings
levelNum <- 32L
contiDim <- 3L

numList <- c(16)
resultList <- list()
nodeList <- list()

TBGPmseRecord <- numeric(1)
TBGPtimeRecord <- numeric(1)
TBGPpredictTimeRecord <- numeric(1)

QQGPmseRecord <- numeric(1)
QQGPtimeRecord <- numeric(1)
QQGPpredictTimeRecord <- numeric(1)

tgpmseRecord <- numeric(1)
tgptimeRecord <- numeric(1)

LVGPmseRecord <- numeric(1)
LVGPtimeRecord <- numeric(1)

GPmseRecord <- numeric(1)
GPtimeRecord <- numeric(1)
GPpredictTimeRecord <- numeric(1)

tData <- read.csv("trainData.csv")
tData$BaseMaterial <- as.character(tData$BaseMaterial)
tData$FinMaterial <- as.character(tData$FinMaterial)
tData$FanType <- as.character(tData$FanType)
for(i in 1:4){
  tData[(tData$BaseMaterial ==  unique(tData$BaseMaterial)[i]), 4] <- paste0("c1.", i)
  tData[(tData$FinMaterial ==  unique(tData$FinMaterial)[i]), 5] <- paste0("c2.", i)
}
for(i in 1:2){
  tData[(tData$FanType ==  unique(tData$FanType)[i]), 6] <- paste0("c3.", i)
}
names(tData) <- c("x.1", "x.2", "x.3", "c.1", "c.2", "c.3", "Y")
tData$c.1 <- factor(tData$c.1)
tData$c.2 <- factor(tData$c.2)
tData$c.3 <- factor(tData$c.3)
data <- tData



#################################################
# TGP tree / GP fit
print(paste0("TBGP has begun."))
TBGPTime1 <- proc.time()
set.seed(99999)
tree1 <- treeFit.corrMat(data = data, initRange = GPR,
                        BCDloopNum = 3, mcmcIterNum = 50, nugget = sqrt(.Machine$double.eps))
print("initial complete.")
treeFinal <- pruneFun(tree1, data)
print("pruning complete.")
leaf.GP1 <- leafGPFun(data = data, GPRange = GPR, QQfiRange = GPR, QQthetaRange = QQGPR,
                     tree = treeFinal, NRF = F)
TBGPTime2 <- proc.time()
TBGPTime2 - TBGPTime1

plotTree(treeFinal)

print(paste0("Building completed."))
TBGPpredictTime1 <- proc.time()
diff <- predictFun(newPoint = data.test, treeFinal, leaf.GP1, parallel = F)$predictValueVec - data.test[,dim(data.test)[2]]
mse <- sum(diff^2)/dim(data.test)[1]
TBGPpredictTime2 <- proc.time()
TBGPperform <- c(mse = mse, CPUtime = (TBGPTime2 - TBGPTime1)[3], PredictTime = (TBGPpredictTime2 - TBGPpredictTime1)[3])

TBGPmseRecord <- TBGPperform[1]
TBGPtimeRecord <- TBGPperform[2]
TBGPpredictTimeRecord <- TBGPperform[3]

print(paste0("TBGP has been completed."))




###############################################
# QQGP fit
QQGPdata <- coerceFun.ori(data = data, design = data[,1:3])
print(paste0("QQGP has begun."))
QQGPTime1 <- proc.time()
set.seed(99999)
QQGPfit <- QQGP.ori.model(design = QQGPdata[,c(1:3)], data = QQGPdata,
                 initParam = rep(-2, 3), fiRange = GPR, thetaRange = QQGPR, Tmatrix = NULL,
                 nugget = sqrt(.Machine$double.eps))
QQGPTime2 <- proc.time()
QQGPTime2 - TBGPTime1

QQGPpredictTime1 <- proc.time()
diff <- QQGP.ori.prediction(QQGP.model = QQGPfit, newPoint = coerceFun.ori(data = data.test, design = data.test[,1:3]))$predictValue - data.test[,dim(data.test)[2]]
mse <- sum(diff^2)/dim(data.test)[1]
QQGPpredictTime2 <- proc.time()
QQGPperform <- c(mse = mse, CPUtime = (QQGPTime2 - QQGPTime1)[3], PredictTime = (QQGPpredictTime2 - QQGPpredictTime1)[3])


QQGPmseRecord <- QQGPperform[1]
QQGPtimeRecord <- QQGPperform[2]
QQGPpredictTimeRecord <- QQGPperform[3]

print(paste0("QQGP has been completed."))


##############################################################################
# tgp training data setting
Imat <- matrix(0, ncol = 31, nrow = dim(data)[1])
for(i in 1:dim(data)[1]){
  Imat[i,(which(data[i,]$c.1 == c.1.level) + 4 * (which(data[i,]$c.2 == c.2.level)-1) + 16 * (which(data[i,]$c.3 == c.3.level)-1)-1)] <- 1
}

Xdat <- cbind(data$x.1, data$x.2, data$x.3, Imat)
Y <- data$Y

# tgp testing data setting
pred.Imat <- matrix(0, ncol = 31, nrow = dim(data.test)[1])

for(i in 1:dim(data.test)[1]){
  pred.Imat[i,(which(data.test[i,]$c.1 == c.1.level) + 4 * (which(data.test[i,]$c.2 == c.2.level)-1) + 16 * (which(data.test[i,]$c.3 == c.3.level)-1)-1)] <- 1
}
pred.Xdat <- cbind(data.test$x.1, data.test$x.2, data.test$x.3, pred.Imat)
Y.test <- data.test$Y


print(paste0("tGP has begun."))

tgpTime1 <- proc.time()
fit1 <- btgp(X = Xdat, XX= pred.Xdat, Z = Y, verb = 0, meanfn = "constant", basemax = 34, nug.p=0, gd=c(1e-8, 1))
tgpTime2 <- proc.time()
tgpmse <- sum((fit1$ZZ.mean - data.test[, dim(data.test)[2]])^2)/dim(pred.Imat)[1]
tgpperform <- c(mse=tgpmse, CPUtime = (tgpTime2 - tgpTime1)[3])

tgpmseRecord <- tgpperform[1]
tgptimeRecord <- tgpperform[2]

print(paste0("tGP has been completed."))


##############################################################################
# LVGP training data setting
len <- 16
X_tr <- matrix(0, ncol = (3+1), nrow = dim(Xdat)[1])
X_tr[,1:3] <- Xdat[,1:3]
for(row in (len+1):dim(Xdat)[1]){
 X_tr[row, (3+1)] <- which(Xdat[row, ((3+1):(3+dim(Imat)[2]))] == 1) + 1
}
X_tr[1:len, (3+1)] <- 1
Y_tr <- Y

# LVGP testing data setting
X_te <- matrix(0, ncol = (3+1), nrow = dim(pred.Xdat)[1])
X_te[,1:3] <- pred.Xdat[,1:3]
for(row in 17:dim(pred.Xdat)[1]){
   X_te[row, (3+1)] <- which(pred.Xdat[row, ((3+1):(3+dim(pred.Imat)[2]))] == 1) + 1
}
X_te[1:16, (3+1)] <- 1
Y_te <- Y.test
n_te <- nrow(X_te)

print(paste0("LVGP has begun."))

LVGPTime1 <- proc.time()
model <- LVGP_fit(X_tr, Y_tr, ind_qual = c(dim(X_tr)[2]), parallel = TRUE)
LVGPTime2 <- proc.time()
LVGPTime2 - LVGPTime1

output <- LVGP_predict(X_te, model)
Y_hat <- output$Y_hat

LVGPmse <- sum((Y_hat-Y_te)^2)/n_te
LVGPperform <- c(mse=LVGPmse, CPUtime = (LVGPTime2 - LVGPTime1)[3])

LVGPmseRecord <- LVGPperform[1]
LVGPtimeRecord <- LVGPperform[2]

print(paste0("LVGP has been completed."))


print(paste0("GP has begun."))
GPmseRecordSingle <- numeric(16*32)

GPTime1 <- proc.time()
for(t in 1:32){
  assign(paste0("GPdata", t), data[(t-1)*numList+(1:numList), c(1:3,7)])
  assign(paste0("GPmodel",t), GP.model(get(paste0("GPdata", t)), range = GGPR,
                                       initParam = rep(0, 3), nugget = sqrt(.Machine$double.eps)))
}
GPTime2 <- proc.time()

GPpredictTime1 <- proc.time()
for(t in 1:32){
  assign(paste0("GPpredict", t), GP.prediction(get(paste0("GPmodel",t)), newPoint = data.test[(t-1)*16+(1:16), c(1:3,7)]))
  GPmseRecordSingle[(t-1)*16+(1:16)] <- get(paste0("GPpredict", t))$predictValue - data.test[(t-1)*16+(1:16), c(7)]
}
GPpredictTime2 <- proc.time()

mse <- sum(GPmseRecordSingle^2)/dim(data.test)[1]

GPperform <- c(mse = mse, CPUtime = (GPTime2 - GPTime1)[3], PredictTime = (GPpredictTime2 - GPpredictTime1)[3])

GPmseRecord <- GPperform[1]
GPtimeRecord <- GPperform[2]
GPpredictTimeRecord <- GPperform[3]

print(paste0("GP has been completed."))


RMSE.result <- data.frame(GP = sqrt(GPmseRecord), TBGP = sqrt(TBGPmseRecord), QQGP = sqrt(QQGPmseRecord),
                          tGP = sqrt(tgpmseRecord), LVGP = sqrt(LVGPmseRecord),
                          row.names = 'RMSE')
fitTime.result <- data.frame(GP = GPtimeRecord, TBGP = TBGPtimeRecord, QQGP = QQGPtimeRecord, tGP = tgptimeRecord, LVGP = LVGPtimeRecord)
predTime.result <- data.frame(GP = GPpredictTimeRecord, TBGP = TBGPpredictTimeRecord, QQGP = QQGPpredictTimeRecord)


RMSE.result
