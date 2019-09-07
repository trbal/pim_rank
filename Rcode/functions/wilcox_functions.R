wilcox.pim <- function(formula, data, h0 = 0.5,alternative="two.sided",...){
  if(!(alternative %in% c("two.sided","less","greater")))
    stop("'alternative' should be 'two.sided','less' or 'greater'")
  if(missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  if(!((length(formula[[2]]) == 1) && (length(formula[[3]]) == 1)))
    stop("'formula' should be of the form response ~ groups")
  if(missing(data))
    stop("'data' is missing")
  
  text_form <- as.character(formula)

  y <- text_form[2]
  x <- text_form[3]
  
  if(nlevels(factor(eval(parse(text = paste("data$",x, sep=""))))) != 2L)
    stop("'groups' must have exact 2 levels")
  
  if (any(diff(c(length(eval(parse(text = paste("data$",y, sep="")))), 
                 length(eval(parse(text = paste("data$",x, sep="")))))) != 0L))
    stop("'response' and 'groups' must have the same length")
  
  comp <- expand.grid(which(factor(eval(parse(text = paste("data$",x, sep="")))) == levels(factor(eval(parse(text = paste("data$",x, sep="")))))[1]), 
              which(factor(eval(parse(text = paste("data$",x, sep="")))) == levels(factor(eval(parse(text = paste("data$",x, sep="")))))[2]))
  
  
  # restrict comparisons within block
  #compare <- comp[eval(parse(text = paste("data$",blocks, sep="")))[comp[,1]] == 
  #                  eval(parse(text = paste("data$",blocks, sep="")))[comp[,2]],]
  
  #new_formula <- as.formula(paste(y, "~ 1"))
  pim.score <- pim(formula, data=data, compare = comp, link="identity",vcov.estim = score.vcov)
  
  lv1 <- paste("t(combn(levels(as.factor(data$",x,")),2))[,1]",sep="")
  lv2 <- paste("t(combn(levels(as.factor(data$",x,")),2))[,2]",sep="")
  newnames <- paste(eval(x),
                    eval(parse(text =lv1)),
                    eval(parse(text =lv2)),sep="")  
  names(pim.score@coef) <- colnames(pim.score@vcov) <- rownames(pim.score@vcov) <- newnames
  
  #pim.score@coef <- pim.score@coef - 0.5/ncol(model.matrix(pim.score))
  overall_test <- (coef(pim.score)-0.5)%*%ginv(sqrt(vcov(pim.score)))
  #pval_overall <- 1 - pchisq(overall_test,df_overall)
  
  pim.wald <- pim(formula, data=data, compare = comp, link="identity",vcov.estim = sandwich.vcov)
  wald <- (coef(pim.wald)-0.5)%*%ginv(sqrt(vcov(pim.wald)))
  names(pim.wald@coef) <- colnames(pim.wald@vcov) <- rownames(pim.wald@vcov) <- newnames
  if(alternative=="two.sided"){
    pval_overall <- (1 - pnorm(abs(overall_test)))*2
    p.wald <- (1 - pnorm(abs(wald)))*2
  } else {
    if(alternative=="greater"){
      pval_overall <- (1-pnorm(overall_test))
      p.wald <- (1-pnorm(wald))
    } else {
      pval_overall <- pnorm(overall_test)
      p.wald <- pnorm(wald)
    }
  }
  
  new('wilcox.pim',
      formula = formula,
      coef = coef(pim.score),
      se.standard = summary(pim.score, h0=h0)@se,
      se.wald = summary(pim.wald, h0=h0)@se,
      zval = summary(pim.score, h0=h0)@zval,
      zwald = summary(pim.wald, h0=h0)@zval,
      pr = summary(pim.score, h0=h0)@pr,
      prwald = summary(pim.wald, h0=h0)@pr,
      z_overall = as.numeric(overall_test),
      pr_z = as.numeric(pval_overall),
      wald = as.numeric(wald),
      pr_w = as.numeric(p.wald),
      h0 = h0,
      alternative = alternative
  )
}

setClass(
  'wilcox.pim',
  slots = c(formula = 'formula',
            coef = 'numeric',
            se.standard = 'numeric',
            se.wald = 'numeric',
            zval = 'numeric',
            zwald = 'numeric',
            pr = 'numeric',
            prwald = 'numeric',
            z_overall = 'numeric',
            pr_z = 'numeric',
            wald = 'numeric',
            pr_w = 'numeric',
            h0 = 'numeric',
            alternative = 'character'
  )
)

summary.wilcox.pim <- function(object,method,
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
  
  if(method == "Wald"){
    z_overall = object@wald
    pr_z = object@pr_w
  } else {
    z_overall = object@z_overall
    pr_z = object@pr_z
  }
  
  
  new("wilcox.pim.summary",
      formula=object@formula,
      coef = object@coef,
      se.standard = object@se.standard,
      se.wald = object@se.wald,
      zval = object@zval,
      zwald = object@zwald,
      pr = object@pr,
      prwald = object@prwald,
      z_overall = z_overall,
      pr_z = pr_z,
      h0 = object@h0,
      method = method,
      alternative = object@alternative
  )
}


setClass(
  'wilcox.pim.summary',
  slots=c(formula='formula',
          coef = 'numeric',
          se.standard = 'numeric',
          se.wald = 'numeric',
          zval = 'numeric',
          zwald = 'numeric',
          pr = 'numeric',
          prwald = 'numeric',
          z_overall = 'numeric',
          pr_z = 'numeric',
          h0 = 'numeric',
          method = 'character',
          alternative = 'character'
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

print.wilcox.pim.summary <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  orig <- paste(deparse(formula(x@formula)))
  cat("Summary of following PIM : \n Wilcoxon-Mann-Whitney \n\nFormula: ", orig, "\n\n")
  
  
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
  
  cat("\nz-value =", x@z_overall, ", p-value =", x@pr_z,"\n\n")
  
  cat("Null hypothesis: b = ", x@h0, "\nAlternative = ", x@alternative,"\n")
  
  if(x@method == 'Wald')
    cat("\nKeep in mind:\n The standard errors, z-values and corresponding p-values are based on the Wald-type covariance matrix. Difference in significance may occur.")
  
}


setMethod("show",
          "wilcox.pim.summary",
          function(object){
            print.wilcox.pim.summary(object)
          })


print.wilcox.pim <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  ...){
  orig <- paste(deparse(x@formula))
  coefs <- x@coef
  
  cat('\nProbabilistic Index Model:\n Wilcoxon-Mann-Whitney\n\nFormula: ',orig,
      "\n\n")
  
  if (length(coefs)) {
    cat("Coefficients:\n")
    print.default(format(coefs, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  
  cat("\nz-value =", x@z_overall, ", p-value =", x@pr_z)
  cat("\nNull hypothesis: b = ", x@h0, "\nAlternative = ", x@alternative,"\n")
  
  invisible(NULL)
}

setMethod('show',
          'wilcox.pim',
          function(object){print(object)})

setMethod('print',
          'wilcox.pim',
          print.wilcox.pim)

coef.wilcox.pim <- coef.wilcox.pim.summary <- function(object,...){
  object@coef
}

setMethod('coef',
          'wilcox.pim',
          coef.wilcox.pim)


setMethod('coef',
          'wilcox.pim.summary',
          coef.wilcox.pim.summary)


confint.wilcox.pim <- function(object, parm, level = 0.95, method, ...){
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

confint.wilcox.pim.summary <- function(object, parm, level = 0.95, method,...){
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
          'wilcox.pim',
          confint.wilcox.pim)

setMethod('confint',
          'wilcox.pim.summary',
          confint.wilcox.pim.summary)

