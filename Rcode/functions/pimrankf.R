pim.rankf <- function(formula,
                data,
                object,
                new.modmat,
                group = group,
                covar = covar,
                na.action = getOption("na.action"),
                weights=NULL,
                keep.data = FALSE,
                ...)
  {
  
  # Check the arguments
  model <- "difference"
  
  nodata <- missing(data)
  #vcov <- match.fun(vcov)
  link <- "identity"
  
  if(is.null(na.action)) na.action <- "na.fail"
  if(!is.character(na.action)) 
    na.action <- deparse(substitute(na.action))
  
  # Check formula and extract info
  
  
  # Create the pim environment (similar to model frame)
  
  penv <- object@penv
  
  x <- new.modmat  
  
  y <- object@response
  
  res <- pim.fit(x, y, link, weights = weights, 
                 penv = as.environment(penv@poset),...)
  

  # as.environment will only pass the environment of the penv to avoid
  # copying the whole thing. makes it easier to get the poset out
  
  names(res$coef) <- colnames(x)
  
  if(!keep.data){
    x <- matrix(nrow=0,ncol=0)
    y <- numeric(0)
  } 
  
  pim:::new.pim(
    formula = object@formula,
    coef = res$coef,
    vcov = res$vcov,
    fitted = res$fitted,
    penv = penv,
    link = link,
    estimators=res$estim,
    model.matrix = x,
    na.action = na.action,
    response = y,
    keep.data = keep.data,
    model = model)
}