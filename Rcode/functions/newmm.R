new.mm.pim <- function(data, group, covar=NA,...){
  data <- droplevels(data)
  lv1 <- paste("t(combn(levels(as.factor(data$",group,")),2))[,1]",sep="")
  lv2 <- paste("t(combn(levels(as.factor(data$",group,")),2))[,2]",sep="")
  newnames <- paste(eval(group),
                    eval(parse(text =lv1)),
                    eval(parse(text =lv2)),sep="")

  new.mm <- matrix(0, ncol=length(newnames)+3,nrow=ncol(combn(nrow(data),2)))
  colnames(new.mm) <- c(newnames,"dup","Var1","Var2")  
  cnt <- 0
  for(i in 1:(nrow(data)-1)){
    lvl1 <- data[i,eval(group)]
    for(j in (i+1):nrow(data)){
      cnt <- cnt+1
      lvl2 <- data[j,eval(group)]
      nmlv <- c(paste(eval(group),lvl1,lvl2,sep=""),paste(eval(group),lvl2,lvl1,sep=""))
      defnm <- which(nmlv %in% newnames)
      new.mm[cnt,"Var1"] <- i
      new.mm[cnt,"Var2"] <- j
      if(!identical(defnm,integer(0))){
        defnm <- nmlv[defnm]
        new.mm[cnt,defnm] <- 1
        new.mm[cnt,"Var1"] <- ifelse(nmlv[1]==defnm, i, j)
        new.mm[cnt,"Var2"] <- ifelse(nmlv[1]==defnm, j, i)
      } else {
        new.mm[cnt,"dup"] <- 1
      }
    }
    
  }
  if(is.matrix(covar)){
      new.mm <- cbind(new.mm, covar)
  }
  
  return(new.mm)
}
