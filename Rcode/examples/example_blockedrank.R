#first clear the environment
rm(list = ls())
#load all needed packages and functions
library(pim)
library(MASS)
##your working directory should be the Rcode folder
source("functions/blockedrank_functions.R")
source("functions/rankpim.R")
#load data
load("data/blockedrank.RData")

#quick check of the data
summary(data)
table(data$lexicality, data$language)


#the test
##pim version
fit <- blockedrank.pim(rt~lexicality|language, data=data)
summary(fit)
confint(fit)

##MS-version
rm(data)
load("data/MS.RData")
fit <- blockedrank.pim(rt~lexicality|language, data=data)
summary(fit)
confint(fit)


##Friedman-version
rm(data)
load("data/friedman.RData")
fit <- blockedrank.pim(rt~lexicality|language, data=data)
summary(fit)
confint(fit)
