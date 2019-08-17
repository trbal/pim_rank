library(pim)
library(MASS)

source("functions/rankcovar_functions.R")
source("functions/newmm.R")
source("functions/pimrankf.R")

load("data/sampleblockdata.RData")

fit <- rank_covar.pim(rt ~ lexicality | block, data = data.sm)
summary(fit)


fit2 <- rank_covar.pim(rt ~ lexicality, data = data.sm)
summary(fit2)


fit3 <- rank_covar.pim(rt ~ nsyl | block, data = data.sm)
summary(fit3)


fit4 <- rank_covar.pim(rt ~ nsyl, data = data.sm)
summary(fit4)
