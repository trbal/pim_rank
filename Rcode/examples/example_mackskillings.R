#first clear the environment
rm(list = ls())
#load all needed packages and functions
library(pim)
library(MASS)
library(asbio)
##your working directory should be the Rcode folder
source("functions/mackskillings_functions.R")
source("functions/rankpim.R")
#load data
load("data/MS.RData")

#quick check of the data
summary(data)
table(data$lexicality, data$language)


#the tests
##classic MS
dataNL <- data$rt[which(data$language=="NL")]
dataEN <- data$rt[which(data$language=="EN")]
dataFR <- data$rt[which(data$language=="FR")]
dataTOT <- cbind(dataNL, dataEN, dataFR)
trt <- c(rep(1,50),rep(2,50))

MS.test(dataTOT, trt, reps=50)

##pim version
fit <- mackskillings.pim(rt~lexicality|language, data=data)
summary(fit)
confint(fit)