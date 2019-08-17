#load data per trial
load("/root/Desktop/CNSDAN_Thesis_Template_tex/data/dlp-trials.Rdata")
#load data per stimulus
load("/root/Desktop/CNSDAN_Thesis_Template_tex/data/dataNL.RData")


findnsyl <- function(x){
  dataNL[which(dataNL$spelling == data.block$spelling[j]),"nsyl"]
}

for(i in 51:58){
  cat(i,"\n")
  #extracte data from block 1 and remove NA
  data.block <- na.omit(dlp.trials[dlp.trials$block==i,c("spelling","lexicality","rt","participant")])
  data.block$participant <- droplevels(data.block$participant)
  #add nsyl in data.block
  nsyl <- c()
  dat.sp <- data.block$spelling
  for(j in 1:nrow(data.block)){
    cat(i,"  ",j,"\n")
    nsyl[j] <- dataNL[which(dataNL$spelling == dat.sp[j]),"nsyl"]
  }
  filen <- paste("fullblock",i,".RData", sep="")
  save(data.block,nsyl,file = filen)
  #remove outliers per participant per nsyl
  for(j in 1:length(levels(data.block$participant))){
    pt <- which(data.block$participant == levels(data.block$participant)[j])
    med <- median(data.block$rt[pt])
    mda <- mad(data.block$rt[pt])
    rmpt <- intersect(pt, which(abs((data.block$rt - med)/mda) > 3))
    if(!identical(rmpt,integer(0))){
      data.block <- data.block[-rmpt,]
      nsyl <- nsyl[-rmpt]
    }
      
    rownames(data.block) <- c()
  }
  
  #
  data.block <- aggregate(x = data.block$rt,by =  list(spelling = data.block$spelling, lexicality = data.block$lexicality, nsyl = nsyl),FUN = mean, simplify = TRUE, drop = TRUE)
  names(data.block)[4] <- "rt"
  
  filen <- paste("aggblock",i,".RData",sep="")
  data.fin <- data.frame(lexicality = data.block$lexicality, nsyl = data.block$nsyl, rt = data.block$rt)
  save(data.fin, file = filen)
  
}
