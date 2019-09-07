jt.pim <- function(formula, data, alternative = "two.sided",...){
  if(missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  if(!((length(formula[[2]]) == 1) && (length(formula[[3]]) == 1)))
    stop("'formula' should be of the form response ~ group")
  if(missing(data))
    stop("'data' is missing")
  
  options(warn=-1)
  
  y <- formula[[2]]
  x <- formula[[3]]
  
  if(!(is.ordered(eval(parse(text = paste("data$",x, sep="")))) && is.factor(eval(parse(text = paste("data$",x, sep=""))))))
    stop(paste("Following variables should be an ordered factor:",
               pim:::.lpaste(as.character(x))))
  
  text_form <- as.character(formula)
  pnames_old <- c(text_form[3], text_form[2])
  text_form[3] <- paste("I((L(", text_form[3], ") < R(", text_form[3],
                        ")) - (L(", text_form[3], ") > R(", text_form[3], "))) + 1", sep="")
  new_formula <- as.formula(paste(text_form[2], text_form[1], text_form[3]))
  
  x <- eval(parse(text = paste("data$",x, sep="")))
  y <- eval(parse(text = paste("data$",y, sep="")))
  
  groups <- lapply(levels(x), function(object) y[which(x==object)])
  k <- length(groups)
  ns <- sapply(groups, length)
  JT <- 0
  for(i in c(1:(k-1))){
    for(j in c(1:ns[i])){
      for(cnt in c((i+1):k)){
        JT <- JT + sum(groups[[i]][j] < groups[[cnt]]) + sum(groups[[i]][j] == groups[[cnt]])/2
      }
    }
  }
  
  pim.score <- pim(formula = new_formula, data = data,
                   link = "identity", vcov.estim = score.vcov,
                   compare = "all")
  
  pim.wald <- pim(formula = new_formula, data = data,
                  link = "identity", vcov.estim = sandwich.vcov,
                  compare = "all")  
  
  overall_test <- (coef(pim.score)[2])/sqrt(vcov(pim.score)[2,2])
  if(alternative=="two.sided"){
    pval_overall <- (1 - pnorm(abs(overall_test)))*2
  } else {
    if(alternative=="greater"){
      pval_overall <- (1-pnorm(overall_test))
    } else {
      pval_overall <- pnorm(overall_test)
    }
  }
  
  cf <- coef(pim.score)
  se.wald <- summary(pim.wald)@se
  
  names(cf)[2] <- names(se.wald)[2] <-
    paste("P(",pnames_old[2],"i <= ",pnames_old[2],"j|",pnames_old[1],"i < ", pnames_old[1], "j)", sep="")
  cf[2] <- cf[2] + 0.5
  
  new('jt.pim',
      formula = formula,
      coef = cf[2],
      se.standard = summary(pim.score)@se[2],
      se.wald = se.wald[2],
      zval = summary(pim.score)@zval[2],
      pr = summary(pim.score)@pr[2],
      JT = JT,
      zval_overall = overall_test,
      pval_overall = pval_overall,
      alternative = alternative)
}





setClass(
  'jt.pim',
  slots = c(formula = 'formula',
            coef = 'numeric',
            se.standard = 'numeric',
            se.wald = 'numeric',
            zval = 'numeric',
            pr = 'numeric',
            JT = 'numeric',
            zval_overall = 'numeric',
            pval_overall = 'numeric',
            alternative = 'character'
  )
)

summary.jt.pim <- function(object,method,
                           ...){
  if(missing(method) || (method == "default")){
    method = 'default'
  } else {
    if(method == "Wald"){
      method == "Wald"
    } else {
      stop("'method' is not correctly defined")
    }
  }
  
  new("jt.pim.summary",
      formula=object@formula,
      coef = object@coef,
      se.standard = object@se.standard,
      zval = object@zval,
      pr = object@pr,
      JT = object@JT,
      zval_overall = object@zval_overall,
      pval_overall = object@pval_overall,
      alternative = object@alternative,
      method = method
  )
}

setClass(
  'jt.pim.summary',
  slots=c(formula='formula',
          coef = 'numeric',
          se.standard = 'numeric',
          zval = 'numeric',
          pr = 'numeric',
          JT = 'numeric',
          zval_overall = 'numeric',
          pval_overall = 'numeric',
          alternative = 'character',
          method = 'character'
  ),
  validity = function(object){
    
    if(!pim:::.equal.lengths(
      object@coef,
      object@se.standard,
      object@zval,
      object@pr
    )){
      stop("coef, se, zval and pr should be of equal length")
    } else {
      TRUE
    }
    
  }
)

print.jt.pim.summary <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  orig <- paste(deparse(formula(x@formula)))
  cat("Summary of following PIM : \n Jonckheere-Terpstra \n\nFormula: ", orig, "\n\n")
  
  Tab <- cbind(
    Estimate = x@coef, 
    "Std. Error" = x@se.standard,
    "z value" = x@zval, 
    "Pr(>|z|)" = x@pr
  )
  
  if(x@method == 'Wald'){
    Tab[,2] = x@se.wald
    Tab[,3] = x@zwald
    Tab[,4] = x@prwald
  }
  
  cat("\n")
  printCoefmat(Tab, digits = digits)
  
  cat("\nJonckheere-Terpstra Statistic = ", x@JT,"\n")
  cat("z-value = ", x@zval_overall, "p-value = ", x@pval_overall, "\n\n")
  
  cat("Null hypothesis: P(Yi <= Yj|Xi < Xj) = 0.5\n")
  cat("Alternative = ", x@alternative, "\n")
  
  if(x@method == 'Wald')
    cat("\nKeep in mind:\n The standard errors, z-values and corresponding p-values are based on the Wald-type covariance matrix. Difference in significance may occur.")
  
}


setMethod("show",
         "jt.pim.summary",
         function(object){
           print.jt.pim.summary(object)
         })

print.jt.pim <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 ...){
  orig <- paste(deparse(x@formula))
  coefs <- x@coef
  
  cat('\nProbabilistic Index Model:\n Jonckheere-Terpstra\n\nFormula: ',orig,
      "\n\n")
  
  if (length(coefs)) {
    cat("Coefficients:\n")
    print.default(format(coefs, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  
  cat("\nJonckheere-Terpstra Statistic =", x@JT, ", z-value =", x@zval[2], ", p-value =", x@pval_overall, "\n")
  cat("Null hypothesis: P(Yi <= Yj|Xi < Xj) = 0.5\n")
  cat("Alternative = ", x@alternative, "\n")
  invisible(NULL)
}

setMethod('show',
          'jt.pim',
          function(object){print(object)})

setMethod('print',
          'jt.pim',
          print.jt.pim)

coef.jt.pim <- coef.jt.pim.summary <- function(object,...){
  object@coef
}

setMethod('coef',
          'jt.pim',
          coef.jt.pim)


setMethod('coef',
          'jt.pim.summary',
          coef.jt.pim.summary)

confint.jt.pim <- function(object, parm, level = 0.95, method, ...){
  # This code is almost literally copied from confint.default
  # because that one doesn't work with S4
  # kinda stupid, but well... 
  cf <- coef(object)
  pnames <- names(cf)
  if(missing(parm))
    parm <- pnames
  else if(is.numeric(parm))
    parm <- pnames[parm]
  
  if(missing(method) || (method == "Wald")){
    method = "Wald" 
  } else {
    if(method == "default")
      method = method
    else
      stop("'method' is not correctly defined")
  }
  
  
  a <- (1-level)/2
  a <- c(a, 1-a)
  
  pct <- stats:::format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                             pct))
  if(method=="Wald")
    ses <- object@se.wald[parm]
  else{
    ses <- object@se.standard[parm]
    warning("The confidence intervals are not obtained with the Wald statistic and may be inaccurate.")
  }
  
  
  ci[] <- cf[parm] + ses %o% fac
  ci
}

confint.jt.pim.summary <- function(object, parm, level = 0.95, method,...){
  # This code is almost literally copied from confint.default
  # because that one doesn't work with S4
  # kinda stupid, but well... 
  cf <- coef(object)
  pnames <- names(cf)
  if(missing(parm))
    parm <- pnames
  else if(is.numeric(parm))
    parm <- pnames[parm]
  
  if(missing(method) || (method == "Wald")){
    method = "Wald" 
  } else {
    if(method == "default")
      method = method
    else
      stop("'method' is not correctly defined")
  }
  
  a <- (1-level)/2
  a <- c(a, 1-a)
  
  pct <- stats:::format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                             pct))
  if(method=="Wald")
    ses <- object@se.wald[parm]
  else{
    ses <- object@se.standard[parm]
    warning("The confidence intervals are not obtained with the Wald statistic and may be inaccurate.")
  }
  
  ci[] <- cf[parm] + ses %o% fac
  ci
}


setMethod('confint',
          'jt.pim',
          confint.jt.pim)

setMethod('confint',
          'jt.pim.summary',
          confint.jt.pim.summary)
