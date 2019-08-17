rank.pim <- function(data, x, y, compare, vcov.estim, link = "identity"){
  formula <- as.formula(paste(y, "~", x))
  
  pim(formula, data = data, compare = compare,
      link = link, vcov.estim = vcov.estim, model = "marginal")
}