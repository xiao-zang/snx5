.libPaths()
start <- proc.time()

require(lme4)
require(foreach)
require(doMC)
registerDoMC(cores=16)
load("./GeneID.RData")
GeneID$Symbol <- as.character(GeneID$Symbol)
GeneID$BARCODE <- as.character(GeneID$BARCODE)

dat.sum <- read.csv("./preliminary_summary.csv", stringsAsFactors = FALSE)
dat <- list()
for(i in 1:9){
  dat.tem <- dat.sum[, c(51, 2:5, (5*i+1):(5*i+5))]
  names(dat.tem)[6:10] <- c("estimate", "std.error", 
                            "z.value", "p.value", "cell.number")
  dat.tem$replicate <- i
  dat[[i]] <- dat.tem
}
dat.col <- rbind(dat[[1]], dat[[2]], dat[[3]], 
                 dat[[4]], dat[[5]], dat[[6]], 
                 dat[[7]], dat[[8]], dat[[9]])
dat.col$virus <- rep(c("mock", "rna", "dna"), each=nrow(dat.col)/3)
alwaysMockwells <- dat.col$WellID[1:96][dat.col$Gene.Symbol[1:96] %in% c("Atg7_Mock", "NC_Mock")]
dat.col$virus[dat.col$Gene.Symbol %in% c("Atg7_Mock", "NC_Mock")] <- "mock"
dat.col$Gene.Symbol <- sub("(.*?)_Mock", "\\1", dat.col$Gene.Symbol)
rm(dat.tem)
rm(dat)
rm(dat.sum)
expwells <- dat.col$WellID[1:96][!dat.col$Gene.Symbol[1:96] %in% c("NC", "Atg7")]
Atg7wells <- dat.col$WellID[1:96][dat.col$Gene.Symbol[1:96] %in% "Atg7"]
NCwells <- dat.col$WellID[1:96][dat.col$Gene.Symbol[1:96] %in% "NC"]

libplates <- unique(dat.col$Barcode)
n.libplate <- length(libplates)
platefolders <- substring(libplates, first = 8)

effects <- 
  foreach(libplate = 1:n.libplate, .combine = rbind) %:%  
  foreach(i = 1:81, .combine = rbind) %dopar% {
    print(proc.time() - start)
    require(lme4)
    require(MIPHENO)
    platedir <- paste("./csv", platefolders[libplate], sep="/")
    repfiles <- list.files(platedir)
    plateid <- substring(platefolders[libplate], first = 2)
    csvs <- list()
    id <- GeneID[GeneID$UTSWPlateNo == plateid, c(2:4, 8:9), drop=FALSE]
    for (j in 1:9){
      csv <- read.csv(paste(platedir, repfiles[j], sep="/"), stringsAsFactors = FALSE)
      csv <- merge(csv[,1:10], id, by.x="Sci_WellID", by.y="UTSWWellNo", all.x=TRUE)
      csv$Symbol[csv$Sci_WellID %in% NCwells] <- "NC"
      csv$Symbol[csv$Sci_WellID %in% Atg7wells] <- "Atg7"
      csv$replicate <- j
      csvs[[j]] <- csv
    }
    dat.currlib <- do.call(rbind, csvs)
    dat.currlib$virus <- "mock"
    dat.currlib$virus[dat.currlib$replicate %in% 4:6] <- "rna"
    dat.currlib$virus[dat.currlib$replicate %in% 7:9] <- "dna"
    dat.currlib$virus[dat.currlib$Sci_WellID %in% alwaysMockwells] <- "mock"
    cellnum.currlib <- dat.col[substring(dat.col$Barcode, first=9)==plateid,c("Barcode", "WellID", "cell.number", "replicate")]
    numcheck <- TRUE
    if(i == 81){
      cellnum.currgene <- cellnum.currlib[cellnum.currlib$WellID %in% c(NCwells, Atg7wells),,drop=FALSE]
      cellnum.currgene$treatment <- "control"
      cellnum.currgene$treatment[cellnum.currgene$WellID %in% Atg7wells] <- "transfect"
      cellnum.currgene <- cellnum.currgene[cellnum.currgene$cell.number > 100,,drop=FALSE]
      if(nrow(cellnum.currgene)==0){
        numcheck <- FALSE
      }else{
        cellnum.currgene$virus <- "mock"
        cellnum.currgene$virus[cellnum.currgene$replicate %in% c(4,5,6)] <- "vir"
        cellnum.currgene$virus[cellnum.currgene$replicate %in% c(7,8,9)] <- "vir"
        cellnum.currgene$virus[cellnum.currgene$WellID %in% alwaysMockwells] <- "mock"
        if(!any(cellnum.currgene$virus=="mock") | !any(cellnum.currgene$treatment == "control")){
          numcheck <- FALSE
        }
      }
      if(numcheck){
        dat.currgene <- dat.currlib[dat.currlib$Sci_WellID %in% c(NCwells, Atg7wells),,drop=FALSE]
        geneinfo <- dat.currgene[dat.currgene$Sci_WellID %in% Atg7wells,,drop=FALSE][1,c(1:2, 11:14),drop=FALSE]
        names(cellnum.currgene)[2] <- "Sci_WellID"
        cellnum.currgene <- cellnum.currgene[,c(2,4)]
        dat.currgene <- merge(cellnum.currgene, dat.currgene)
      }
      
    }else{
      cellnum.currgene <- cellnum.currlib[cellnum.currlib$WellID %in% c(NCwells, expwells[i]),,drop=FALSE]
      cellnum.currgene$treatment <- "control"
      cellnum.currgene$treatment[cellnum.currgene$WellID %in% expwells[i]] <- "transfect"
      cellnum.currgene <- cellnum.currgene[cellnum.currgene$treatment == "transfect" | cellnum.currgene$cell.number > 100,,drop=FALSE]
      if(nrow(cellnum.currgene)==0){
        numcheck <- FALSE
      }else{
        cellnum.currgene$virus <- "mock"
        cellnum.currgene$virus[cellnum.currgene$replicate %in% c(4,5,6)] <- "vir"
        cellnum.currgene$virus[cellnum.currgene$replicate %in% c(7,8,9)] <- "vir"
        cellnum.currgene$virus[cellnum.currgene$WellID %in% alwaysMockwells] <- "mock"
        if(!any(cellnum.currgene$virus=="mock") | !any(cellnum.currgene$treatment == "control")){
          numcheck <- FALSE
        }
      }
      if(numcheck){
        dat.currgene <- dat.currlib[dat.currlib$Sci_WellID %in% c(NCwells, expwells[i]),,drop=FALSE]
        geneinfo <- dat.currgene[dat.currgene$Sci_WellID %in% expwells[i],,drop=FALSE][1,c(1:2, 11:14),drop=FALSE]
        names(cellnum.currgene)[2] <- "Sci_WellID"
        cellnum.currgene <- cellnum.currgene[,c(2,4)]
        dat.currgene <- merge(cellnum.currgene, dat.currgene)
      }
    }
      if(!numcheck)
      {
        fit <- list(-3,-3)
      }else{
        dat.currgene$siRNA <- "knockdown"
        dat.currgene$siRNA[dat.currgene$Sci_WellID %in% NCwells] <- "control"
        dat.currgene$replicate <- as.factor(dat.currgene$replicate)
        dat.currgene$virus <- relevel(as.factor(dat.currgene$virus), "mock")
        dat.currgene$siRNA <- relevel(as.factor(dat.currgene$siRNA), "control")
        outliergroup <- split(dat.currgene, list(dat.currgene$replicate, dat.currgene$virus, dat.currgene$siRNA), drop=TRUE)
        for (k in 1:length(outliergroup)){
          curr.out <- outliergroup[[k]]
          if(length(unique(curr.out$Sci_WellID)) < 3) next
          todetect <- curr.out[,c("Sci_WellID", "number_green_dots")]
          todetect$debug <- 1 # counteract a bug in rm.outliers (drop=FALSE missing)
          todetect <- rm.outliers(todetect, parameter = 'Sci_WellID', n=3)
          curr.out <- curr.out[curr.out$Sci_WellID %in% todetect$Sci_WellID[!is.na(todetect$number_green_dots)],]
          outliergroup[[k]] <- curr.out
        }
      dat.currgene <- do.call(rbind, args = outliergroup)
      dat.currgene$Sci_WellID <- as.factor(dat.currgene$Sci_WellID)  
      dat.currgene$obs <- 1:nrow(dat.currgene)
      dat.currgene$obs <- as.factor(dat.currgene$obs)
      fit <- tryCatch(glmer(number_green_dots ~ 1 + siRNA + virus+ virus:siRNA + (1 | replicate) + (1 | replicate:Sci_WellID), data=dat.currgene, family=poisson), error=function(e){list(-1, conditionMessage(e))}, warning=function(e){list(-2, conditionMessage(e))})
      }
    flatresheader <- c("(Intercept).coef", "(Intercept).sd(coef)", "(Intercept).Z", "(Intercept).pvalue", "siRNAknockdown.coef", 
                       "siRNAknockdown.sd(coef)", "siRNAknockdown.Z", "siRNAknockdown.pvalue", "virusdna.coef","virusdna.sd(coef)",
                       "virusdna.Z",                       "virusdna.pvalue",                  "virusrna.coef",                    "virusrna.sd(coef)",                "virusrna.Z", 
                       "virusrna.pvalue",                  "siRNAknockdown:virusdna.coef",     "siRNAknockdown:virusdna.sd(coef)", "siRNAknockdown:virusdna.Z",        "siRNAknockdown:virusdna.pvalue", 
                       "siRNAknockdown:virusrna.coef",     "siRNAknockdown:virusrna.sd(coef)", "siRNAknockdown:virusrna.Z",        "siRNAknockdown:virusrna.pvalue",   "modelstatus"
    )
    if(is.list(fit)){
      flatres <- data.frame(matrix(data=NA, nrow=1, ncol= 25))
      names(flatres) <- flatresheader
	  if(fit[[1]] == -1){
	    flatres$modelstatus <- paste("error:", fit[[2]])
	  }else if(fit[[1]] == -2){
	    flatres$modelstatus <- paste("warning:", fit[[2]])
      }else if(fit[[1]] == -3){
	    flatres$modelstatus <- "checknum = FALSE; no mock or NC"
	  }
    }else{
      res <- data.frame(summary(fit)$coefficients)
      names(res) <- c("coef", "sd(coef)", "Z", "pvalue")
      flatres <- data.frame(t(as.matrix(as.vector(t(res)))))
      flatresnameframe <- expand.grid(names(res), row.names(res))
      names(flatres) <- paste(flatresnameframe[,2], flatresnameframe[,1], sep=".")
      flatres$modelstatus <- "success"
    }
    formatframe <- data.frame(matrix(nrow=0, ncol=25))
    names(formatframe) <- flatresheader
    flatres <- merge(flatres, formatframe, all=TRUE)
    flatres <- flatres[,match(flatresheader,names(flatres))]
    return(data.frame(geneinfo, flatres))
  }
write.csv(effects, "second_analysis.csv")
elapsed <- proc.time() - start
sessionInfo()