#first clear the environment
rm(list = ls())
#load all needed packages and functions
library(pim)
library(MASS)
library(Skillings.Mack)
##your working directory should be the Rcode folder
source("functions/skillingsmack_functions.R")
source("functions/rankpim.R")

#load data
load("data/SM.RData")

#the tests
##classic SM
Ski.Mack(data$rt, groups = data$lexicality, blocks = data$language)

##pim version
fit <- skillingsmack.pim(rt~lexicality|language, data=data)
summary(fit)
confint(fit)
Ski.Mack
