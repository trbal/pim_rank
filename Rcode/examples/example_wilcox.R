#first clear the environment
rm(list = ls())
#load all needed packages and functions
library(pim)
library(MASS)
##your working directory should be the Rcode folder
source("functions/wilcox_functions.R")
#load data
load("data/kruskall.RData")


summary(data)

#the tests
##classic KW
wilcox.test(rt~lexicality, data = data)
##pim version
fit <- wilcox.pim(rt~lexicality, data = data)
summary(fit)
confint(fit)

