#first clear the environment
rm(list = ls())
#load all needed packages and functions
library(pim)
library(MASS)
##your working directory should be the Rcode folder
source("functions/kruskall_functions.R")
source("functions/rankpim.R")
#load data
load("data/kruskall.RData")


summary(data)

#the tests
##classic KW
kruskal.test(rt~nsyl, data = data)
##pim version
#pseudo - 1 syl - 2 syl
fit <- kw.pim(rt~nsyl, data = data)
summary(fit)
confint(fit)

#lexicality
fit <- kw.pim(rt~lexicality, data = data)
summary(fit)
confint(fit)
