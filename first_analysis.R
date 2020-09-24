.libPaths()
#####################
### data analysis ###
#####################
require(lme4)
require(MIPHENO)
require(foreach)
require(doMC)
registerDoMC(cores=16)


files <- list.files("./csv")

### load Gene ID data ###
load("./GeneID.RData")
GeneID$Symbol <- as.character(GeneID$Symbol)
GeneID$BARCODE <- as.character(GeneID$BARCODE)
dir.create("./first_analysis")
### functions ###

# Object Plate #
Plate=function(x, y, cn, rn){
  y <- as.numeric(y)
  rgb.palette <- colorRampPalette(c("green","black","red"), space = "rgb")
  y.len <- length(y)
  color.code.y <- y - min(y, na.rm=TRUE)  
  color.code.y <- color.code.y/max(color.code.y, na.rm=TRUE)  
  color.code.y <- floor(color.code.y*(y.len - 1))+1
  color.code.y <- rgb.palette(y.len)[color.code.y]  
  return(list(x=x, y=y, cn=cn, rn=rn, color.code.y=color.code.y))
}

# Plot Plate #
plate.plot <- function(object, x.axis=NULL, y.axis=NULL, cex=10, cex.axis=1,
                       xlab="", ylab="", cex.lab=1.5, main=NULL){
  loc <- object$x
  y <- object$y
  color.code.y <- object$color.code.y
  cn <- object$cn
  rn <- object$rn
  rgb.palette <- colorRampPalette(c("green","black","red"), space = "rgb")
  par(mar=c(5, 5, 4, 2))
  plot(c(0.5,cn+0.5), c(0.5,rn+0.5), type="n", xlab=xlab, ylab=ylab,
       xaxt="n", yaxt="n", cex.lab=cex.lab,  main=main)
  points(loc[, 1], loc[, 2], col=color.code.y, pch=20, cex=cex)
  if(!is.null(x.axis)){
    axis(1, at=1:cn, lab=x.axis, cex.axis=cex.axis)
  }
  if(!is.null(y.axis)){
    axis(2, at=1:rn, lab=y.axis, cex.axis=cex.axis)
  }
}

# heatmap legend #
heatmap.legend <- function(x)
{
  par(mar <- c(1,1,1,4),cex.axis <- 1.5)
  ix <- 1
  nlevel=16
  minz <- min(x, na.rm=TRUE)
  maxz <- max(x, na.rm=TRUE)
  binwidth <- (maxz - minz)/nlevel
  midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iy <- midpoints
  iz <- matrix(iy, nrow = 1, ncol = length(iy))
  rgb.palette <- colorRampPalette(c("green","black","red"), 
                                  space = "rgb")
  image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "", col = rgb.palette(13))
  axis(4, mgp = c(3, 1, 0), las = 2)
}

# boxplot for cell number
boxplot.cellnum <- function(y, x, ...){
  plot(y~x, type="n", xaxt='n', ...)
  points(y~x, ...)
  axis(1, at=1:12, 1:12, ...)
}

for(i in files){
  print(i)
  plateid <- as.integer(substr(i, 2, 5))
  dir.create(paste("./first_analysis/", i, sep=""))
  dir.create(paste("./first_analysis/", i, "/plot", sep=""))
  dir.create(paste("./first_analysis/", i, "/result", sep=""))
  inputfiles <- list.files(paste("./csv/", i, sep=""))
  
  # check if directory contains 9 files
  if(length(inputfiles) != 9){
    stop("abnormal inputfiles number")
  }
  
  fit.df <- list()
  for(j in 1:9){
	cat("\nlibrary ", i,"replicate ", j)
    ### data input ###    
    filename <- unlist(strsplit(inputfiles[j], "\\."))[1]
    scdat <- read.csv(paste("./csv/", i, "/", inputfiles[j], sep=""), stringsAsFactors=FALSE)
    scdat$row <- substr(scdat$Sci_WellID, 1, 1)
    scdat$col <- substr(scdat$Sci_WellID, 2, 3)
    
    ### quality control ###
    
    ### positive and negetive contorl ###
    
    scdat.sub <- scdat[scdat$col %in% c("01", "12"), ]
    scdat.sub$well <- paste(scdat.sub$col, scdat.sub$row, sep="")
    jpeg(file=paste("./first_analysis/", i, "/plot/", filename, "_quality.jpeg", sep=""), width=1000, height=600)
    par(mar=c(5, 5, 2, 2))
    boxplot(log2(number_green_dots+1)~well, data=scdat.sub, 
            col=c(rep(2:3, 4), rep(3:2, 4)),notch=T,
            xlab="Colume Row ID", ylab="Log2 of dots number", 
            cex.lab=2, cex.axis=1.5)
    dev.off()
    
    ### cell nunber outlier detection ###
    
    cell.num <- aggregate(scdat$totROI_ObjectCount,
                          by=list(Row=scdat$row, Col=scdat$col), mean, na.rm=T)
    names(cell.num)[3] <- "cell.num"
    cell.num$row <- 9-match(cell.num$Row, LETTERS)
    cell.num$col <- as.integer(cell.num$Col)
    cell.num$well.id <- paste(cell.num$Row, cell.num$Col, sep="")
    loc <- cell.num[, c("col", "row")]
    p <- Plate(x=loc, y=cell.num$cell.num, 12, 8)
    
    jpeg(file=paste("./first_analysis/", i, "/plot/", filename, "_cell number_boxplot.jpeg", sep=""), 
         width=1000, height=600, quality=100)
    par(mar=c(5, 5, 2, 2))
    boxplot.cellnum(cell.num$cell.num, cell.num$col, col=cell.num$row, pch=20, cex=2, 
                    cex.axis=1.5, xlab="Column", ylab="Cell Number", cex.lab=2, 
                    xlim=c(-1, 12))
    legend("topleft", paste("Row", LETTERS[1:8]), col=8:1, bty="n", cex=1.5, pch=20)
    dev.off()
    
    
    jpeg(file=paste("./first_analysis/", i, "/plot/", filename, "_cell number.jpeg", sep=""), 
         width=1000, height=800, quality=100)
    par(mar=c(10, 10, 2, 2))
    plate.plot(p, xlab="Col", ylab="Row", x.axis=1:12, cex=19,
               y.axis=rev(c("A", "B", "C", "D", "E", "F", "G", "H")), 
               cex.axis=2, cex.lab=2.5)
    dev.off()
    
    jpeg(file=paste("./first_analysis/", i, "/plot/", filename, "_cell number_legend.jpeg", sep=""), 
         width=200, height=800)
    heatmap.legend(p$y)
    dev.off()
    
    # analysis #
    agg.mean <- aggregate(log2(scdat$number_green_dots+1), 
                          by=list(row=scdat$row, col=scdat$col),
                          mean, na.rm=T)
    loc <- cbind(rep(1:12, each=8), rep(8:1, time=12))
    p <- Plate(x=loc, y=agg.mean[, 3], 12, 8)
    jpeg(file=paste("./first_analysis/", i, "/plot/", filename, "_plate.jpeg", sep=""), 
         width=1000, height=800)
    par(mar=c(10, 10, 2, 2))
    plate.plot(p, xlab="Col", ylab="Row", x.axis=1:12, cex=19,
               y.axis=rev(c("A", "B", "C", "D", "E", "F", "G", "H")), 
               cex.axis=2, cex.lab=2.5)
    dev.off()
    
    jpeg(file=paste("./first_analysis/", i, "/plot/", filename, "_plate_legend.jpeg", sep=""), 
         width=200, height=800)
    heatmap.legend(p$y)
    dev.off()
    
    #######################################
    ### Generalized Mixed Effects model ###
    #######################################
    
    scdat.poss <- scdat
    scdat.poss$id <- as.character(scdat.poss$Sci_WellID)
    scdat.nc <- scdat.poss[scdat.poss$col == 12 & 
                             scdat.poss$row %in% c("A", "C", "E", "G"), ]
    scdat.nc <- scdat.nc[scdat.nc$totROI_ObjectCount > 100 ,,drop=FALSE]
    if(nrow(scdat.nc) == 0){
      num.status <- FALSE
      curr.fit <- list(-3, -3)
    }else{
      num.status <- TRUE
      if(length(unique(scdat.nc$Sci_WellID)) > 2){
        todetect <- scdat.nc[,c("Sci_WellID", "number_green_dots")]
        todetect$debug <- 1
        todetect <- rm.outliers(todetect, parameter = "Sci_WellID", n = 3)
        scdat.nc <- scdat.nc[scdat.nc$Sci_WellID %in% todetect$Sci_WellID[!is.na(todetect$number_green_dots)],,drop=FALSE]
      }
      scdat.nc$id <- "00000"
    }
    wells <- unique(scdat.poss$Sci_WellID)
    fit.coef <- foreach(gi = 1:length(wells), .combine = rbind) %dopar% {
      require(lme4)
      if(num.status){
        curr.scdat.poss <- rbind(scdat.poss[scdat.poss$Sci_WellID == wells[gi],], scdat.nc)
        curr.scdat.poss$obs <- 1:nrow(curr.scdat.poss)
        curr.scdat.poss$Sci_WellID <- as.factor(curr.scdat.poss$Sci_WellID)
        curr.scdat.poss$obs <- as.factor(curr.scdat.poss$obs)
        curr.scdat.poss$id <- relevel(as.factor(curr.scdat.poss$id), ref = "00000")

        curr.fit <- tryCatch(glmer(number_green_dots ~ 1 + id + (1 | Sci_WellID), data=curr.scdat.poss, family=poisson()),
                            error=function(e){list(-1, conditionMessage(e))}, warning=function(e){list(-2, conditionMessage(e))})
 
      }
      coefformat <- data.frame(matrix(data=NA, nrow=1, ncol=5))
      names(coefformat) <- c("Estimate", "Std..Error", "z.value", "Pr...z..", "well")
      if(!is.list(curr.fit)){
        curr.fit.coef <- data.frame(summary(curr.fit)$coefficient)
        curr.fit.coef <- curr.fit.coef[row.names(curr.fit.coef) != "(Intercept)", ]
        curr.fit.coef$well <- substr(row.names(curr.fit.coef), 3, 5)
        curr.fit.coef$modelstatus <- "success"
      }else if(curr.fit[[1]]==-1){
        curr.fit.coef <- coefformat
        curr.fit.coef$modelstatus <- paste("error:", curr.fit[[2]])
        curr.fit.coef$well <- wells[gi]
      }else if(curr.fit[[1]]==-2){
        curr.fit.coef <- coefformat
        curr.fit.coef$modelstatus <- paste("warning:", curr.fit[[2]])
        curr.fit.coef$well <- wells[gi]
      }else if(curr.fit[[1]]==-3){
        curr.fit.coef <- coefformat
        curr.fit.coef$modelstatus <- "noNCpassCellnum"
        curr.fit.coef$well <- wells[gi]
      }
      return(curr.fit.coef)
    }
    
    ### Gene ID mergeing ###
    
    id <- GeneID[GeneID$UTSWPlateNo == plateid, ]
    
    fit.coef <- merge(fit.coef, id, by.x="well", by.y="UTSWWellNo", 
                      all.x=T)
    fit.coef$Symbol <- as.character(fit.coef$Symbol)
    fit.coef$Symbol[fit.coef$well  %in% c("A01", "C01", "E01", "G01")] <- "Atg7_Mock"
    fit.coef$Symbol[fit.coef$well  %in% c("B01", "D01", "F01", "H01")] <- "NC_Mock"
    fit.coef$Symbol[fit.coef$well  %in% c("A12", "C12", "E12", "G12")] <- "NC"
    fit.coef$Symbol[fit.coef$well  %in% c("B12", "D12", "F12", "H12")] <- "Atg7"
    
    ### cell number ###
    cell.agg <- aggregate(scdat$totROI_ObjectCount,
                          by=list(scdat$Sci_WellID),mean, na.rm=T)
    
    fit.coef <- merge(fit.coef, cell.agg, by.x="well", by.y="Group.1", all.x=T)
    names(fit.coef)[15] <- "Total.Cell.Number"
    
    if(!all(is.na(fit.coef$Estimate))){
      fit.p <- fit.coef[, c("well", "Estimate")]
      fit.p$Row <- substr(fit.p$well, 1, 1)
      fit.p$Col <- substr(fit.p$well, 2, 3)
      fit.p$row <- 9-match(fit.p$Row, LETTERS)
      fit.p$col <- as.integer(fit.p$Col)
      loc <- fit.p[, c("col", "row")]
      p <- Plate(x=loc, y=fit.p$Estimate, 12, 8)
      
      jpeg(file=paste("./first_analysis/", i, "/plot/", filename, "_coefficients.jpeg", sep=""), 
           width=1000, height=800, quality=100)
      par(mar=c(10, 10, 2, 2))
      plate.plot(p, xlab="Col", ylab="Row", x.axis=1:12, cex=19,
                 y.axis=rev(c("A", "B", "C", "D", "E", "F", "G", "H")), 
                 cex.axis=2, cex.lab=2.5)
      dev.off()
      
      jpeg(file=paste("./first_analysis/", i, "/plot/", filename, "_coefficients_legend.jpeg", sep=""), 
           width=200, height=800)
      heatmap.legend(p$y)
      dev.off()
      
      jpeg(file=paste("./first_analysis/", i, "/plot/", filename, "_coefficients_boxplot.jpeg", sep=""), 
           width=1000, height=600, quality=100)
      par(mar=c(5, 5, 2, 2))
      boxplot.cellnum(fit.p$Estimate, fit.p$col, col=fit.p$row, pch=20, cex=2, 
                      cex.axis=1.5, xlab="Column", ylab="Coefficient", cex.lab=2, 
                      xlim=c(-1, 12))
      legend("topleft", paste("Row", LETTERS[1:8]), col=8:1, bty="n", cex=1.5, pch=20)
      dev.off()
    }   
    
    fit.df[[j]] <- fit.coef
    write.csv(fit.coef, file=paste("./first_analysis/", i, "/result/", filename, "_coefficient.csv", sep=""))
  }
  
  dat <- fit.df[[1]]
  dat <- dat[, c("well","Symbol", "Locus.ID","Accession", "Estimate", "Std..Error", "z.value", "Pr...z..", "modelstatus",  "Total.Cell.Number" 
                  )]
  for(k in 2:9){
    a <- fit.df[[k]]
    a <- a[, c("well","Symbol", "Locus.ID","Accession", "Estimate", "Std..Error", "z.value", "Pr...z..", "modelstatus",  "Total.Cell.Number" 
    )]
    dat <- cbind(dat, a)
  }
  
  if(all(apply(dat[, seq(1, 90, 10)], 1, FUN=function(x){
    return(all(as.character(x) == as.character(x[1])))
  }))){
    dat <- dat[, -c(seq(11, 90, 10), seq(12, 90, 10), 
                    seq(13, 90, 10), seq(14, 90, 10))]
    names(dat)[1:4] <- c("WellID", "Gene.Symbol", 
                         "Gene.ID", "Gene.Accession")
    names(dat)[5:58] <- paste(rep(paste("Plate", 1:9, sep="."), each=6), 
                              rep(c("Estimate", "Std.Error", "Z.Value", 
                                    "P.Value", "modelstatus" ,"Cell.Number"), times=9), 
                              sep=".")
    dat$Barcode <- id$BARCODE[1]
  }else{
    stop("merging error!")
  }
  write.csv(dat, paste("./first_analysis/", i, "/", plateid, "_all.csv", sep=""), 
            row.names=F)
}

#######################
### merge all files ###
#######################

dat <- list()
j <- 1
for(i in list.files("./first_analysis")){
  print(i)
  file.name <- list.files(paste("./first_analysis/", i, sep=""))
  file.name <- file.name[grep("all", file.name)]
  dat[[j]] <- read.csv(paste("./first_analysis/", i, "/", file.name, sep=""))
  j <- j+1
}

res <- dat[[1]]
for(i in 2:length(dat)){
  res <- rbind(res, dat[[i]])
}

write.csv(res, "first_analysis.csv")
sessionInfo()
