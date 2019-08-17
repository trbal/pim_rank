#first clear the environment
rm(list = ls())
#load all needed packages and functions
library(pim)
library(MASS)
library(clinfun)
##your working directory should be the Rcode folder
source("functions/jt_functions.R")
#load data
load("data/jt.RData")


summary(data)
data$nsyl


#the tests
##classic KW
jonckheere.test(data$rt, data$nsyl)
##pim version
fit <- jt.pim(rt~nsyl, data=data)
summary(fit)
