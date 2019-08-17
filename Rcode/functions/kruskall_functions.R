kw.pim <- function(formula, data, h0 = 0.5,...){
  nodata <- missing(data)
  if(missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  if(!((length(formula[[2]]) == 1) && (length(formula[[3]]) == 1)))
    stop("'formula' should be of the form response ~ group")
  if(missing(data))
    stop("'data' is missing")
    
  
  left_hand <- formula[[2]]
  right_hand <- formula[[3]]
  
  
  f.terms <- terms(formula, simplify = TRUE)
  vars <- all.vars(formula)
  
  text_form <- as.character(formula)
  
  x <- eval(parse(text = paste("data$",text_form[3], sep="")))
  data[,text_form[3]] <- factor(x)
  
  
  pnames_old <- text_form[3]

  pim.score <- rank.pim(data, text_form[3], text_form[2], "all", score.vcov)


  overall_test <- t(coef(pim.score) - h0)%*%ginv(vcov(pim.score))%*%c(coef(pim.score)-h0)
  df <- length(coef(pim.score))-1
  pscore <- 1 - pchisq(overall_test, df = df)
  
  pim.wald <- rank.pim(data, text_form[3], text_form[2], "all", sandwich.vcov)
  
  wald <- t(coef(pim.wald) - h0)%*%ginv(vcov(pim.wald))%*%c(coef(pim.wald) - h0)
  p.wald <- 1 - pchisq(wald, ncol(vcov(pim.wald)))
  
  cf <- coef(pim.score)
  names(cf) <- paste(pnames_old, c(1:length(cf)), sep="")
  
  vcov.wald <- vcov(pim.wald)
  colnames(vcov.wald) <- rownames(vcov.wald) <- names(cf)
  
  se.wald <- summary(pim.wald, h0 = h0)@se
  names(se.wald) <- paste(pnames_old, c(1:length(cf)), sep="")
  
  new('kw.pim',
      formula = formula,
      coef = cf,
      se.standard = summary(pim.score, h0=h0)@se,
      se.wald = summary(pim.wald, h0=h0)@se,
      zval = summary(pim.score, h0=h0)@zval,
      zwald = summary(pim.wald, h0=h0)@zval,
      pr = summary(pim.score, h0=h0)@pr,
      prwald = summary(pim.wald, h0=h0)@pr,
      chi_sq = as.numeric(overall_test),
      df_cs = as.numeric(df),
      pr_cs = as.numeric(pscore),
      wald = as.numeric(wald),
      df_w = as.numeric(ncol(vcov(pim.wald))),
      pr_w = as.numeric(p.wald),
      vcov.wald = vcov.wald,
      h0 = h0
      )
}

setClass(
  'kw.pim',
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
            df_w = 'numeric',
            pr_w = 'numeric',
            vcov.wald = 'matrix',
            h0 = 'numeric'
  )
)

summary.kw.pim <- function(object,method,
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
  
  new("kw.pim.summary",
      formula=object@formula,
      coef = object@coef,
      se.standard = object@se.standard,
      se.wald = object@se.wald,
      zval = object@zval,
      zwald = object@zwald,
      pr = object@pr,
      prwald = object@prwald,
      chi_sq = object@chi_sq,
      df_cs = object@df_cs,
      pr_cs = object@pr_cs,
      h0 = object@h0,
      method = method
  )
}


setClass(
  'kw.pim.summary',
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

print.kw.pim.summary <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  orig <- paste(deparse(formula(x@formula)))
  cat("Summary of following PIM : \n Kruskal-Wallis \n\nFormula: ", orig, "\n\n")
  
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
          "kw.pim.summary",
          function(object){
            print.kw.pim.summary(object)
          })

print.kw.pim <- function(x, digits = max(3L, getOption("digits") - 3L),
                       ...){
  orig <- paste(deparse(x@formula))
  coefs <- x@coef

  cat('\nProbabilistic Index Model:\n Kruskal-Wallis\n\nFormula: ',orig,
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
          'kw.pim',
          function(object){print(object)})

setMethod('print',
          'kw.pim',
          print.kw.pim)


coef.kw.pim <- coef.kw.pim.summary <- function(object,...){
  object@coef
}

setMethod('coef',
          'kw.pim',
          coef.kw.pim)


setMethod('coef',
          'kw.pim.summary',
          coef.kw.pim.summary)


confint.kw.pim <- function(object, parm, level = 0.95, method, ...){
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

confint.kw.pim.summary <- function(object, parm, level = 0.95, method,...){
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
          'kw.pim',
          confint.kw.pim)

setMethod('confint',
          'kw.pim.summary',
          confint.kw.pim.summary)

