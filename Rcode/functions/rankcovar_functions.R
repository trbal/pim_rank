rank_covar.pim <- function(formula, data, h0 = 0.5, ...){
  if(missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  if(!((length(formula[[2]]) == 1) && ((length(formula[[3]]) == 3)||(length(formula[[3]]) == 1))))
    stop("'formula' should be of the form response ~ groups or response~groups|covariate")
  if(missing(data))
    stop("'data' is missing")
  options(warn=-1)
  text_form <- as.character(formula)
  if(grepl("\\|",text_form[3])){
    text_form[3] <- gsub("\\s", "", text_form[3])
    text_form[3:4] <- unlist(strsplit(text_form[3],split='\\|'))
    group <- text_form[3]
    covar <- text_form[4]
    covname <- covar
  } else {
    group <- text_form[3]
    covar <- NA
  }
  
  
  y <- text_form[2]
  
  
  
  new_formula <- as.formula(paste(y, "~", group))
  if(!is.na(covar)){
    formula_covar <- as.formula(paste(y,"~ I(R(", covar,") - L(", covar, "))", sep=""))
    covar <- pim(formula_covar, data=data,link="identity", compare = "unique",keep.data=TRUE)@model.matrix
    colnames(covar) <- covname
  }
  temp.mm <- new.mm.pim(data, group, covar=covar)
  temp.mm <- temp.mm[-which(temp.mm[,"dup"]==1),]
  comp <- temp.mm[,c("Var1","Var2"), drop=FALSE]
  new.mm <- temp.mm[,!colnames(temp.mm) %in% c("dup","Var1","Var2"), drop=FALSE]
  pim.prim <- pim(new_formula, data, link = "identity", compare = comp,keep.data=TRUE)
  pim.score <- pim.rankf(formula, data, pim.prim, new.mm,
                         group, covar = covar, vcov.estim = score.vcov)
  
  df_overall <- length(coef(pim.score))
  if(suppressWarnings(!is.na(covar))){
    pim.score@coef <- pim.score@coef[-length(pim.score@coef)]
    pim.score@vcov <- pim.score@vcov[, -ncol(pim.score@vcov),drop=FALSE]
    pim.score@vcov <- pim.score@vcov[ -nrow(pim.score@vcov), ,drop=FALSE]
  }
  
  #pim.score@coef <- pim.score@coef - 0.5/ncol(model.matrix(pim.score))
  overall_test <- t(coef(pim.score)-0.5)%*%ginv(vcov(pim.score))%*%c(coef(pim.score)-0.5)
  
  pval_overall <- 1 - pchisq(overall_test,df_overall)
  
  
  pim.wald <- pim.rankf(formula, data, pim.prim, new.mm,
                        group, covar = covar, vcov.estim = sandwich.vcov)
  
  if(suppressWarnings(!is.na(covar))){
    pim.wald@coef <- pim.wald@coef[-length(pim.wald@coef)]
    pim.wald@vcov <- pim.wald@vcov[, -ncol(pim.wald@vcov),drop=FALSE]
    pim.wald@vcov <- pim.wald@vcov[ -nrow(pim.wald@vcov), ,drop=FALSE]
  }
  
  new('rank_covar.pim',
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
      h0 = h0
  )
}

setClass(
  'rank_covar.pim',
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
            h0 = 'numeric'
  )
)

summary.rank_covar.pim <- function(object,method,
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
  
  
  new("rank_covar.pim.summary",
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
  'rank_covar.pim.summary',
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

print.rank_covar.pim.summary <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  orig <- paste(deparse(formula(x@formula)))
  cat("Summary of following PIM : \n Rank test with covariate \n\nFormula: ", orig, "\n\n")
  
  
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
          "rank_covar.pim.summary",
          function(object){
            print.rank_covar.pim.summary(object)
          })


print.rank_covar.pim <- function(x, digits = max(3L, getOption("digits") - 3L),
                               ...){
  orig <- paste(deparse(x@formula))
  coefs <- x@coef
  
  cat('\nProbabilistic Index Model:\n Rank test with covariate\n\nFormula: ',orig,
      "\n\n")
  
  if (length(coefs)) {
    cat("Coefficients:\n")
    print.default(format(coefs, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  
  cat("\nchi-squared =", x@chi_sq, ", df =", x@df_cs, ", p-value =", x@pr_cs)
  
  invisible(NULL)
}

setMethod('show',
          'rank_covar.pim',
          function(object){print(object)})

setMethod('print',
          'rank_covar.pim',
          print.rank_covar.pim)

coef.rank_covar.pim <- coef.rank_covar.pim.summary <- function(object,...){
  object@coef
}

setMethod('coef',
          'rank_covar.pim',
          coef.rank_covar.pim)


setMethod('coef',
          'rank_covar.pim.summary',
          coef.rank_covar.pim.summary)


confint.rank_covar.pim <- function(object, parm, level = 0.95, method, ...){
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

confint.rank_covar.pim.summary <- function(object, parm, level = 0.95, method,...){
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
          'rank_covar.pim',
          confint.rank_covar.pim)

setMethod('confint',
          'rank_covar.pim.summary',
          confint.rank_covar.pim.summary)

