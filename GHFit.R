# import the 3-third package
library(ghyp)
library(R.matlab)
# clear the environment 
rm(list = ls())
# read the mat file
tmp=readMat("C:\\Users\\yubin\\Desktop\\GH Distribution\\matlab_x0=1point3.mat")
Data <- matrix(unlist(tmp), ncol = 2)
QList=as.numeric(Data[,1])
# fit the optimal q
QListTrain=QList[sample(length(QList), 8000)]
QListTest=QList[-sample(length(QList), 8000)]
QListTest=QListTest[!is.na(QListTest)]
QFit=fit.ghypuv(data=QListTrain,opt.pars = c(mu = FALSE, sigma = FALSE))
hist(QFit,col="green")
QDensity=dghyp(x=QListTest,object = QFit)
# output the test set log likelihood 
print(sum(log(QDensity))/length(QDensity))

# fit the thetaq
ThetaList=as.numeric(Data[,2])
ThetaListTrain=ThetaList[sample(length(ThetaList), 8000)]
ThetaListTest=ThetaList[-sample(length(ThetaList), 8000)]
ThetaListTest=ThetaListTest[!is.na(ThetaListTest)]
ThetaFit=fit.ghypuv(data=ThetaListTrain,opt.pars = c(mu = FALSE, sigma = FALSE))
hist(ThetaFit,col="pink")
ThetaDensity=dghyp(x=ThetaListTest,object = ThetaFit)
# output the test set log likelihood 
print(sum(log(ThetaDensity))/length(ThetaDensity))
