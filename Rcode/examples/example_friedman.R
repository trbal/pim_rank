#first clear the environment
rm(list = ls())
#load all needed packages and functions
library(pim)
library(MASS)
##your working directory should be the Rcode folder
source("functions/friedman_functions.R")
source("functions/rankpim.R")
#load data
load("data/friedman.RData")

#quick check of the data
summary(data)
table(data$lexicality, data$language)


#the tests
##classic MS
friedman.test(rt ~ lexicality|language, data = data)

##pim version
fit <- friedman.pim(rt~lexicality|language, data=data)
summary(fit)
confint(fit)
