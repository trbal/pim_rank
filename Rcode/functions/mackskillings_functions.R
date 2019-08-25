mackskillings.pim <- function(formula, data, h0 = 0.5,...){
  if(missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  if(!((length(formula[[2]]) == 1) && (length(formula[[3]]) == 3)))
    stop("'formula' should be of the form response ~ groups|blocks")
  if(missing(data))
    stop("'data' is missing")
  
  text_form <- as.character(formula)
  text_form[3] <- gsub("\\s", "", text_form[3])
  text_form[3:4] <- unlist(strsplit(text_form[3],split='\\|'))
  
  y <- text_form[2]
  groups <- text_form[3]
  blocks <- text_form[4]
  
  if (anyNA(eval(parse(text = paste("data$",groups, sep="")))) || 
      anyNA(eval(parse(text = paste("data$",blocks, sep="")))))
    stop("NA's are not allowed in 'groups' or 'blocks'")
  if (any(diff(c(length(eval(parse(text = paste("data$",y, sep="")))), 
                 length(eval(parse(text = paste("data$",groups, sep="")))), 
                 length(eval(parse(text = paste("data$",blocks, sep="")))))) != 0L))
    stop("'response', 'groups' and 'blocks' must have the same length")
  if (any(table(eval(parse(text = paste("data$",groups, sep=""))), 
                eval(parse(text = paste("data$",blocks, sep="")))) - 
          mean(table(eval(parse(text = paste("data$",groups, sep=""))),
                     eval(parse(text = paste("data$",blocks, sep="")))))!=0))
    stop("not an equally replicated complete block design")
  
  comp <- expand.grid(1:nrow(data), 1:nrow(data))
  # restrict comparisons within block
  compare <- comp[eval(parse(text = paste("data$",blocks, sep="")))[comp[,1]] == 
                    eval(parse(text = paste("data$",blocks, sep="")))[comp[,2]],]
  
  new_formula <- as.formula(paste(y, "~", groups))
  
  pim.score <- rank.pim(data, groups, y, compare, score.vcov)
  
  #pim.score@coef <- pim.score@coef - 0.5/ncol(model.matrix(pim.score))
  overall_test <- t(coef(pim.score))%*%ginv(vcov(pim.score))%*%c(coef(pim.score))
  df_overall <- length(coef(pim.score)) - 1
  pval_overall <- 1 - pchisq(overall_test,df_overall)
  
  pim.wald <- rank.pim(data, groups, y, compare, sandwich.vcov)
  wald <- t(coef(pim.wald))%*%ginv(vcov(pim.wald))%*%c(coef(pim.wald))
  p.wald <- 1 - pchisq(wald, df_overall)
  
  new('mackskillings.pim',
      formula = formula,
      coef = coef(pim.score),
      se.standard = summary(pim.score, h0=h0)@se,
      se.wald = summary(pim.wald, h0=h0)@se,
      zval = summary(pim.score, h0=h0)@zval,
      zwald = summary(pim.wald, h0=h0)@zval,
      pr = summary(pim.score, h0=h0)@pr,
      prwald = summary(pim.wald, h0=h0)@pr,
      chi_sq = as.numeric(overall_test),
      df_cs = as.numeric(df_overall),
      pr_cs = as.numeric(pval_overall),
      wald = as.numeric(wald),
      pr_w = as.numeric(p.wald),
      h0 = h0
  )
}

setClass(
  'mackskillings.pim',
  slots = c(formula = 'formula',
            coef = 'numeric',
            se.standard = 'numeric',
            se.wald = 'numeric',
            zval = 'numeric',
            zwald = 'numeric',
            pr = 'numeric',
            prwald = 'numeric',
            chi_sq = 'numeric',
            df_cs = 'numeric',
            pr_cs = 'numeric',
            wald = 'numeric',
            pr_w = 'numeric',
            h0 = 'numeric'
  )
)

summary.mackskillings.pim <- function(object,method,
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
    chi_sq = object@wald
    pr_cs = object@pr_w
  } else {
    chi_sq = object@chi_sq
    pr_cs = object@pr_cs
  }
  
  new("mackskillings.pim.summary",
      formula=object@formula,
      coef = object@coef,
      se.standard = object@se.standard,
      se.wald = object@se.wald,
      zval = object@zval,
      zwald = object@zwald,
      pr = object@pr,
      prwald = object@prwald,
      chi_sq = chi_sq,
      df_cs = object@df_cs,
      pr_cs = pr_cs,
      h0 = object@h0,
      method = method
  )
}


setClass(
  'mackskillings.pim.summary',
  slots=c(formula='formula',
          coef = 'numeric',
          se.standard = 'numeric',
          se.wald = 'numeric',
          zval = 'numeric',
          zwald = 'numeric',
          pr = 'numeric',
          prwald = 'numeric',
          chi_sq = 'numeric',
          df_cs = 'numeric',
          pr_cs = 'numeric',
          h0 = 'numeric',
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

print.mackskillings.pim.summary <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  orig <- paste(deparse(formula(x@formula)))
  cat("Summary of following PIM : \n Mack-Skillings \n\nFormula: ", orig, "\n\n")
  
  
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
  
  cat("\nchi-squared =", x@chi_sq, ", df =", x@df_cs, ", p-value =", x@pr_cs,"\n\n")
  
  cat("Null hypothesis: b = ", x@h0, "\n")
  
  if(x@method == 'Wald')
    cat("\nKeep in mind:\n The standard errors, z-values and corresponding p-values are based on the Wald-type covariance matrix. Difference in significance may occur.")
  
}


setMethod("show",
          "mackskillings.pim.summary",
          function(object){
            print.mackskillings.pim.summary(object)
          })


print.mackskillings.pim <- function(x, digits = max(3L, getOption("digits") - 3L),
                               ...){
  orig <- paste(deparse(x@formula))
  coefs <- x@coef
  
  cat('\nProbabilistic Index Model:\n Mack-Skillings\n\nFormula: ',orig,
      "\n\n")
  
  if (length(coefs)) {
    cat("Coefficients:\n")
    print.default(format(coefs, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  
  cat("\nchi-squared =", x@chi_sq, ", df =", x@df_cs, ", p-value =", x@pr_cs)
  cat("\nNull hypothesis: b = ", x@h0, "\n")
  invisible(NULL)
}

setMethod('show',
          'mackskillings.pim',
          function(object){print(object)})

setMethod('print',
          'mackskillings.pim',
          print.mackskillings.pim)

coef.mackskillings.pim <- coef.mackskillings.pim.summary <- function(object,...){
  object@coef
}

setMethod('coef',
          'mackskillings.pim',
          coef.mackskillings.pim)


setMethod('coef',
          'mackskillings.pim.summary',
          coef.mackskillings.pim.summary)


confint.mackskillings.pim <- function(object, parm, level = 0.95, method, ...){
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

confint.mackskillings.pim.summary <- function(object, parm, level = 0.95, method,...){
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
          'mackskillings.pim',
          confint.mackskillings.pim)

setMethod('confint',
          'mackskillings.pim.summary',
          confint.mackskillings.pim.summary)
