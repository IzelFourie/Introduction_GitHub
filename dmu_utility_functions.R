#############################################################################################
# utility function for using DMU 
#############################################################################################



# Main function for reml analyses using DMU
aireml <- function(fm=NULL,vfm=NULL,Klist=NULL, data=NULL,validate=NULL, bin=NULL) {
     tol <- 0.001
     fit <- cvfit <- model <- NULL
     model <- modelDMU(fm=fm,vfm=vfm, data=data)
     #for (i in 1:length(fm)) {
     #  model[[i]] <- modelDMU(fm=fm,vfm=vfm, data=data)
     #}
     flevels <- writeDMU(model=model,data=data,Klist=Klist)
     executeDMU(bin=bin)
     fit <- readDMU(model=model,flevels=flevels)
     if(!is.null(validate)) {
          for (i in 1:ncol(validate)) {
               #set data missing
               writeDMU(model=model,data=data,Klist=Klist)
               executeDMU(bin=bin)
               cvfit[[i]] <- readDMU(model=model,flevels=flevels)
          }
          fit$cv <- cvfit
     }
     return(fit=fit)
}



# Extract model information from fm and vfm 
modelDMU <- function(fm=NULL,vfm=NULL, data=NULL) {
     
     vtype <- sapply(data,class)
     
     ffvar <- frvar <- vvar <- yvar <- NULL 
     yvar <- as.character(fm)[2]
     
     fmsplit <- unlist(strsplit(as.character(fm)[3],split="+", fixed=TRUE))
     rwsR <- grep(")", fmsplit,fixed=TRUE)
     rwsF <- (1:length(fmsplit))[-rwsR]
     
     fvar <- fmsplit[rwsF]
     fvar <- gsub(" ", "",fvar, fixed=TRUE)
     ffvar <- fvar[vtype[fvar]=="factor"]
     frvar <- fvar[vtype[fvar]=="numeric"]
     
     vvar <- fmsplit[rwsR]
     vvar <- gsub("1 | ", "", vvar, fixed=TRUE)
     vvar <- gsub("+", "",vvar, fixed=TRUE)
     vvar <- gsub("(", "",vvar, fixed=TRUE)
     vvar <- gsub(")", "",vvar, fixed=TRUE)
     vvar <- gsub(" ", "",vvar, fixed=TRUE)
     
     vvar <- lapply(vvar, function(x) {
          x <- unlist( strsplit( as.character(x),split="~", fixed=TRUE)) })
     
     cvvar <- sapply(vvar,function(x){x[2]})
     vvar <- sapply(vvar,function(x){x[1]})
     names(cvvar) <- vvar
     cvvar <- cvvar[!is.na(cvvar)]    
     vvartype <- rep("I",length(vvar))
     names(vvartype) <- vvar
     vvartype[names(cvvar)] <- "COR"
     
     
     nreg <- length(frvar)
     nrandom <- length(vvar)
     nfixed <- length(ffvar)
     nfactors <- nrandom+nfixed 
     
     #vvartype <- rep("I",length(vvar))
     #names(vvartype) <- vvar
     #cvvar <- t(sapply(vfm, function(x) {
     #     x <- unlist( strsplit( as.character(x),split="~", fixed=TRUE))[-1]
     #}))
     #rownames(cvvar) <- cvvar[,1]
     #cvvar <- cvvar[,-1]
     #vvartype[names(cvvar)] <- "COR"
     
     variables <- list( fixed=ffvar,
                        regression=frvar,
                        random=vvar,
                        response=yvar,
                        factors=c(ffvar,vvar),
                        variables=c(ffvar,vvar,yvar,frvar))
     
     n <- as.data.frame(t(sapply(variables,length)))
     
     return(list( fixed=ffvar,
                  regression=frvar,
                  random=vvar,
                  response=yvar,
                  factors=c(ffvar,vvar),
                  variables=c(ffvar,vvar,yvar,frvar),
                  n=n,
                  covtype=vvartype,
                  covmat=cvvar)
     ) 
}

# Recode factors for DMU 
recodeDMU <- function(data=NULL) {
     flevels <- rlevels <- NULL
     for (i in 1:ncol(data)) {
          f <- data[, i]
          flevels[[i]] <- levels(f)
          names(flevels[[i]]) <- 1:nlevels(f)
          rlevels[[i]] <- 1:nlevels(f)
          names(rlevels[[i]]) <- levels(f)
          data[,i] <- names(flevels[[i]])[f]
     }
     names(flevels) <- names(rlevels) <- colnames(data)
     head(data)
     return(list(rlevels=rlevels,flevels=flevels,data=data)) 
}


# Write DIR file, data file, and cor files for DMU 
writeDMU <- function(model=NULL,data=NULL,Klist=NULL,tol=0.001){
     
     # Write DMU DIR file
     dir.file <- "gfm.DIR"
     write("$COMMENT", file = dir.file)
     write("DIR file DMU generated from R ", file = dir.file, 
           append = TRUE)
     write(" ", file = dir.file, append = TRUE)
     write(paste("$ANALYSE", 1,1,0,0), 
           file = dir.file, append = TRUE)
     write(" ", file = dir.file, append = TRUE)
     write(c("$DATA ASCII (", model$n$factors, ",", model$n$response+model$n$regression, 
             ",", -9999, ") data.txt"), file = dir.file, append = TRUE, 
           ncolumns = 12, sep = "")
     write(" ", file = dir.file, append = TRUE)
     write("$VARIABLE", file = dir.file, append = TRUE)
     write(model$variables, file = dir.file, append = TRUE, ncolumns = model$n$variables )
     write(" ", file = dir.file, append = TRUE)
     
     write("$MODEL", file = dir.file, append = TRUE)
     write(1, file = dir.file, append = TRUE)
     write(0, file = dir.file, append = TRUE)
     write(c(1,0,model$n$factors,1:model$n$factors), file = dir.file, append = TRUE, ncolumns = 3+model$n$factors  )
     write(c(model$n$random,1:model$n$random), file = dir.file, append = TRUE, ncolumns = 1+model$n$random  )
     #write(c(length(frvar),1:length(frvar)), file = dir.file, append = TRUE, ncolumns = 1+length(frvar)  )
     write(0, file = dir.file, append = TRUE)
     write(0, file = dir.file, append = TRUE)
     write(" ", file = dir.file, append = TRUE)
     
     for (i in 1:model$n$random) {
          if(model$covtype[i]=="COR") {
               vvarname <- names(model$covtype)[i]
               vvarfile <- paste(vvarname,".txt",sep="")
               
               write(c("$VAR_STR",i,"COR", "ASCII",vvarfile), file = dir.file, append = TRUE, 
                     ncolumns = 5, sep = " ")
          }
     }  
     write(" ", file = dir.file, append = TRUE)
     
     write(c("$DMUAI", format(c(10,0.0000001,0.000001,1,0,0),scientific=FALSE)), 
           file = dir.file, append = TRUE)
     write(" ", file = dir.file, append = TRUE)
     
     write("$RESIDUALS ASCII", file = dir.file, append = TRUE)
     
     
     # Write DMU data file
     data.file <- "data.txt"
     # Recode data (factors only) (check recoding this again)
     rec <- recodeDMU(data[,model$factors])
     data[,model$factors] <- rec$data
     write.table(format(data[,model$variables], scientific = FALSE), file = data.file, 
                 quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
     
     # Write DMU cor data files
     if(!is.null(Klist)) {
          for ( i in 1:length(model$covmat) ) {
               vvarname <- names(model$covmat)[i]
               vvarfile <- paste(vvarname,".txt",sep="")
               iG <- ginv(Klist[[model$covmat[i]]], tol=tol)
               colnames(iG$G) <- rownames(iG$G) <- rec$rlevels[[vvarname]][rownames(iG$G)]
               writeG(G=iG$G, filename=vvarfile, ldet=iG$ldet)
          }
     }
     return(flevels=rec$flevels)
     
}


# Remove DMU output files 
rename.and.clean.windows <- function (jobname=NULL) 
{
     ll.ff <- list.files()
     "my.system" <- function(cmd) return(system(paste(Sys.getenv("COMSPEC"), 
                                                      "/c", cmd),show.output.on.console = FALSE))
     ll.name <- c("SOL", "PAROUT", "PAROUT_STD", "PEDDATA", "RESIDUAL", 
                  "LLIK", "SOL_STD")
     ll <- ll.name[ll.name %in% ll.ff]
     for (kk in ll) my.system(paste("move ", kk, " ", jobname, 
                                    ".", kk, sep = ""))
     junk.files <- c("DMU1.log","DMUAI.log",paste("COR",1:20,sep=""),"CODE_TABLE", "DMU1.dir", "DMUAI.dir", "DMU_LOG", 
                     "DUMMY", "FSPAKWK", "Latest_parm", "LEVAL", "MODINF", 
                     "PARIN", "RCDATA_I", "RCDATA_R", "INBREED", "AINV1", 
                     "AINV2", "PEDFILE1", "PEDFILE2", "fort.81","fort.66")
     del.files <- ll.ff[ll.ff %in% junk.files]
     if (length(del.files)>0) 
          for (kk in 1:length(del.files)) my.system(paste("del ", 
                                                          del.files[kk], sep = ""))
}


# Execute DMU 
executeDMU <- function(bin=NULL) {
     
     jobname <- "gfm"
     
     dmu1 <- "dmu1"
     dmuai <- "dmuai"
     
     if(!is.null(bin)) dmu1 <- paste(bin,"dmu1",sep="/")
     if(!is.null(bin)) dmuai <- paste(bin,"dmuai",sep="/")
     
     out <- paste(jobname, ".dmuai.lst", sep = "")
     dir <- paste(jobname, ".DIR", sep = "")
     "my.system" <- function(cmd) return(system(paste(Sys.getenv("COMSPEC"), 
                                                      "/c", cmd)))
     my.system("set MKL_NUM_THREADS=1")
     test <- my.system(paste(shQuote(dmu1), " < ", dir, " > ", out, sep = ""))
     if (test == 0 & "MODINF" %in% list.files()) {
          test <- my.system(paste(shQuote(dmuai), " < ", dir, " >> ", out, sep = ""))
     }
     rename.and.clean.windows(jobname)
     
}


# Read DMU output files
readDMU <- function(model=NULL, flevels=NULL) {
     
     jobname <- "gfm"
     
     fit <- NULL
     
     llik <- scan(paste(jobname,".LLIK",sep=""), what = character(0), quiet = TRUE)
     
     fit$llik <- as.numeric(llik[12])
     cls1 <- grep("Theta",llik)
     cls2 <- grep("ASD",llik)
     fit$sigma <- as.numeric(llik[(cls1+1):(cls2-1)])
     fit$asd <- as.numeric(llik[(cls2+1):length(llik)])
     names(fit$sigma) <- names(fit$asd) <- c(model$random,"e")
     
     sol <- as.matrix(read.table(paste(jobname,".SOL",sep=""), as.is = TRUE)[-1, , drop = FALSE])
     blue <- sol[sol[,1]==2,c(2,4:6,8:9)]  # "=2" is estimates effects for fixed factors
     blup <- sol[sol[,1]==3,c(2,4:6,8:9)]  # "=3" is estimates effects for random factors
     f <- NULL
     for(i in 1:model$n$random){
          f[[i]] <- blup[blup[,2]==i,5]
          names(f[[i]]) <- flevels[[model$random[i]]][blup[blup[,2]==i,3]]
     }
     names(f) <- model$random
     fit$f <- f
     
     #fit$sigma <- as.matrix(read.table(paste(jobname,".PAROUT",sep="")))
     #fit$asd <- scan(paste(jobname,".PAROUT_STD",sep=""), what = character(0), quiet = TRUE)
     
     resi <- as.matrix(read.table(paste(jobname, ".RESIDUAL", 
                                        sep = "")))
     nc <- ncol(resi)
     nn.per.tr <- 4 #if (object$glmm) 7 else 4
     
     n.trait <- (nc - 1)/nn.per.tr
     if (n.trait != round(n.trait)) 
          stop("something wrong")
     fit$residuals <- resi[, 1 + (((nn.per.tr - 2) * n.trait + 
                                        1):((nn.per.tr - 1) * n.trait))]
     fit$fitted <- resi[, 1 + ((1 * n.trait + 1):(2 * n.trait))]
     fit$hat.matrix <- resi[, 1 + (((nn.per.tr - 1) * n.trait + 
                                         1):(nn.per.tr * n.trait))]
     
     return(fit)
}


# 
#   # Main function for reml analyses suing DMU
#   aireml <- function(fm=NULL,vfm=NULL,corlist=NULL, data=NULL,validate=NULL, tol=0.001, bin=NULL) {
#     fit <- cvfit <- NULL
#     model <- modelDMU(fm=fm,vfm=vfm, data=data)
#     flevels <- writeDMU(model=model,data=data,corlist=corlist)
#     executeDMU(bin=bin)
#     fit <- readDMU(model=model,flevels=flevels)
#     if(!is.null(validate)) {
#       for (i in 1:ncol(validate)) {
#         #set data missing
#         writeDMU(model=model,data=data,corlist=corlist)
#         executeDMU(bin=bin)
#         cvfit[[i]] <- readDMU(model=model,flevels=flevels)
#       }
#       fit$cv <- cvfit
#     }
#     return(fit=fit)
#   }
# 
# 
# 
#   # Extract model information from fm and vfm 
#   modelDMU <- function(fm=NULL,vfm=NULL, data=NULL) {
# 
#     vtype <- sapply(data,class)
#     
#     ffvar <- frvar <- vvar <- yvar <- NULL 
#     yvar <- as.character(fm)[2]
#     
#     fmsplit <- unlist(strsplit(as.character(fm)[3],split="+", fixed=TRUE))
#     rwsR <- grep(")", fmsplit,fixed=TRUE)
#     rwsF <- (1:length(fmsplit))[-rwsR]
#     
#     fvar <- fmsplit[rwsF]
#     fvar <- gsub(" ", "",fvar, fixed=TRUE)
#     ffvar <- fvar[vtype[fvar]=="factor"]
#     frvar <- fvar[vtype[fvar]=="numeric"]
#     
#     vvar <- fmsplit[rwsR]
#     vvar <- gsub("1 | ", "", vvar, fixed=TRUE)
#     vvar <- gsub("+", "",vvar, fixed=TRUE)
#     vvar <- gsub("(", "",vvar, fixed=TRUE)
#     vvar <- gsub(")", "",vvar, fixed=TRUE)
#     vvar <- gsub(" ", "",vvar, fixed=TRUE)
# 
#     vvar <- lapply(vvar, function(x) {
#          x <- unlist( strsplit( as.character(x),split="~", fixed=TRUE)) })
#     
#     cvvar <- sapply(vvar,function(x){x[2]})
#     vvar <- sapply(vvar,function(x){x[1]})
#       names(cvvar) <- vvar
#     cvvar <- cvvar[!is.na(cvvar)]    
#     vvartype <- rep("I",length(vvar))
#      names(vvartype) <- vvar
#     vvartype[names(cvvar)] <- "COR"
#     
#     nreg <- length(frvar)
#     nrandom <- length(vvar)
#     nfixed <- length(ffvar)
#     nfactors <- nrandom+nfixed 
#       
#     #vvartype <- rep("I",length(vvar))
#     # names(vvartype) <- vvar
#     #cvvar <- t(sapply(vfm, function(x) {
#      # x <- unlist( strsplit( as.character(x),split="~", fixed=TRUE))[-1]
#     #}))
#      # rownames(cvvar) <- cvvar[,1]
#      # cvvar <- cvvar[,-1]
#     #vvartype[names(cvvar)] <- "COR"
#     
#     variables <- list( fixed=ffvar,
#                        regression=frvar,
#                        random=vvar,
#                        response=yvar,
#                        factors=c(ffvar,vvar),
#                        variables=c(ffvar,vvar,yvar,frvar))
#     
#     n <- as.data.frame(t(sapply(variables,length)))
#     
#     return(list( fixed=ffvar,
#                  regression=frvar,
#                  random=vvar,
#                  response=yvar,
#                  factors=c(ffvar,vvar),
#                  variables=c(ffvar,vvar,yvar,frvar),
#                  n=n,
#                  covtype=vvartype,
#                  covmat=cvvar)
#     ) 
#   }
# 
#   # Recode factors for DMU 
#   recodeDMU <- function(data=NULL) {
#     flevels <- rlevels <- NULL
#     for (i in 1:ncol(data)) {
#       f <- data[, i]
#       flevels[[i]] <- levels(f)
#       names(flevels[[i]]) <- 1:nlevels(f)
#       rlevels[[i]] <- 1:nlevels(f)
#       names(rlevels[[i]]) <- levels(f)
#       data[,i] <- names(flevels[[i]])[f]
#     }
#     names(flevels) <- names(rlevels) <- colnames(data)
#     head(data)
#     return(list(rlevels=rlevels,flevels=flevels,data=data)) 
#   }
#   
#   
#   # Write DIR file, data file, and cor files for DMU 
#   writeDMU <- function(model=NULL,data=NULL,corlist=NULL,tol=0.001){
#     
#     # Write DMU DIR file
#     dir.file <- "gfm.DIR"
#     write("$COMMENT", file = dir.file)
#     write("DIR file DMU generated from R ", file = dir.file, 
#           append = TRUE)
#     write(" ", file = dir.file, append = TRUE)
#     write(paste("$ANALYSE", 1,1,0,0), 
#           file = dir.file, append = TRUE)
#     write(" ", file = dir.file, append = TRUE)
#     write(c("$DATA ASCII (", model$n$factors, ",", model$n$response+model$n$regression, 
#             ",", -9999, ") data.txt"), file = dir.file, append = TRUE, 
#           ncolumns = 12, sep = "")
#     write(" ", file = dir.file, append = TRUE)
#     write("$VARIABLE", file = dir.file, append = TRUE)
#     write(model$variables, file = dir.file, append = TRUE, ncolumns = model$n$variables )
#     write(" ", file = dir.file, append = TRUE)
#     
#     write("$MODEL", file = dir.file, append = TRUE)
#     write(1, file = dir.file, append = TRUE)
#     write(0, file = dir.file, append = TRUE)
#     write(c(1,0,model$n$factors,1:model$n$factors), file = dir.file, append = TRUE, ncolumns = 3+model$n$factors  )
#     write(c(model$n$random,1:model$n$random), file = dir.file, append = TRUE, ncolumns = 1+model$n$random  )
#     #write(c(length(frvar),1:length(frvar)), file = dir.file, append = TRUE, ncolumns = 1+length(frvar)  )
#     write(0, file = dir.file, append = TRUE)
#     write(0, file = dir.file, append = TRUE)
#     write(" ", file = dir.file, append = TRUE)
#     
#     for (i in 1:model$n$random) {
#       if(model$covtype[i]=="COR") {
#         vvarname <- names(model$covtype)[i]
#         vvarfile <- paste(vvarname,".txt",sep="")
#         
#         write(c("$VAR_STR",i,"COR", "ASCII",vvarfile), file = dir.file, append = TRUE, 
#               ncolumns = 5, sep = " ")
#       }
#     }  
#     write(" ", file = dir.file, append = TRUE)
#     
#     write(c("$DMUAI", format(c(10,0.0000001,0.000001,1,0,0),scientific=FALSE)), 
#           file = dir.file, append = TRUE)
#     write(" ", file = dir.file, append = TRUE)
#     
#     write("$RESIDUALS ASCII", file = dir.file, append = TRUE)
#     
#     
#     # Write DMU data file
#     data.file <- "data.txt"
#     # Recode data (factors only) (check recoding this again)
#     rec <- recodeDMU(data[,model$factors])
#     data[,model$factors] <- rec$data
#     write.table(format(data[,model$variables], scientific = FALSE), file = data.file, 
#                 quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
#     
#     # Write DMU cor data files
#     if(!is.null(corlist)) {
#       for ( i in 1:length(model$covmat) ) {
#         vvarname <- names(model$covmat)[i]
#         vvarfile <- paste(vvarname,".txt",sep="")
#         iG <- ginv(corlist[[model$covmat[i]]], tol=tol)
#         colnames(iG$G) <- rownames(iG$G) <- rec$rlevels[[vvarname]][rownames(iG$G)]
#         writeG(G=iG$G, filename=vvarfile, ldet=iG$ldet)
#       }
#     }
#     return(flevels=rec$flevels)
#     
#   }
# 
# 
#   # Remove DMU output files 
#   rename.and.clean.windows <- function (jobname=NULL) 
#   {
#     ll.ff <- list.files()
#     "my.system" <- function(cmd) return(system(paste(Sys.getenv("COMSPEC"), 
#                                                      "/c", cmd),show.output.on.console = FALSE))
#     ll.name <- c("SOL", "PAROUT", "PAROUT_STD", "PEDDATA", "RESIDUAL", 
#                  "LLIK", "SOL_STD")
#     ll <- ll.name[ll.name %in% ll.ff]
#     for (kk in ll) my.system(paste("move ", kk, " ", jobname, 
#                                    ".", kk, sep = ""))
#     junk.files <- c("DMU1.log","DMUAI.log",paste("COR",1:20,sep=""),"CODE_TABLE", "DMU1.dir", "DMUAI.dir", "DMU_LOG", 
#                     "DUMMY", "FSPAKWK", "Latest_parm", "LEVAL", "MODINF", 
#                     "PARIN", "RCDATA_I", "RCDATA_R", "INBREED", "AINV1", 
#                     "AINV2", "PEDFILE1", "PEDFILE2", "fort.81","fort.66")
#     del.files <- ll.ff[ll.ff %in% junk.files]
#     if (length(del.files)>0) 
#       for (kk in 1:length(del.files)) my.system(paste("del ", 
#                                                       del.files[kk], sep = ""))
#   }
# 
# 
#   # Execute DMU 
#   executeDMU <- function(bin=NULL) {
#     
#     jobname <- "gfm"
# 
#     dmu1 <- "dmu1"
#     dmuai <- "dmuai"
#     
#     if(!is.null(bin)) dmu1 <- paste(bin,"dmu1",sep="/")
#     if(!is.null(bin)) dmuai <- paste(bin,"dmuai",sep="/")
#     
#     out <- paste(jobname, ".dmuai.lst", sep = "")
#     dir <- paste(jobname, ".DIR", sep = "")
#     "my.system" <- function(cmd) return(system(paste(Sys.getenv("COMSPEC"), 
#                                                      "/c", cmd)))
#     my.system("set MKL_NUM_THREADS=1")
#     test <- my.system(paste(shQuote(dmu1), " < ", dir, " > ", out, sep = ""))
#     if (test == 0 & "MODINF" %in% list.files()) {
#       test <- my.system(paste(shQuote(dmuai), " < ", dir, " >> ", out, sep = ""))
#     }
#     rename.and.clean.windows(jobname)
#     
#   }
# 
# 
#   # Read DMU output files
#   readDMU <- function(model=NULL, flevels=NULL) {
#       
#     jobname <- "gfm"
#     
#     fit <- NULL
#     
#     llik <- scan(paste(jobname,".LLIK",sep=""), what = character(0), quiet = TRUE)
#     
#     fit$llik <- as.numeric(llik[12])
#       cls1 <- grep("Theta",llik)
#       cls2 <- grep("ASD",llik)
#     fit$sigmas <- as.numeric(llik[(cls1+1):(cls2-1)])
#     fit$asd <- as.numeric(llik[(cls2+1):length(llik)])
#       names(fit$sigmas) <- names(fit$asd) <- c(model$random,"e")
#     
#     sol <- as.matrix(read.table(paste(jobname,".SOL",sep=""), as.is = TRUE)[-1, , drop = FALSE])
#       blue <- sol[sol[,1]==2,c(2,4:6,8:9)]  # "=2" is estimates effects for fixed factors
#       blup <- sol[sol[,1]==3,c(2,4:6,8:9)]  # "=3" is estimates effects for random factors
#     f <- NULL
#     for(i in 1:model$n$random){
#       f[[i]] <- blup[blup[,2]==i,5]
#       names(f[[i]]) <- flevels[[model$random[i]]][blup[blup[,2]==i,3]]
#     }
#     names(f) <- model$random
#     fit$f <- f
#      
#     #fit$sigma <- as.matrix(read.table(paste(jobname,".PAROUT",sep="")))
#     #fit$asd <- scan(paste(jobname,".PAROUT_STD",sep=""), what = character(0), quiet = TRUE)
# 
#     resi <- as.matrix(read.table(paste(jobname, ".RESIDUAL", 
#                                        sep = "")))
#     nc <- ncol(resi)
#     nn.per.tr <- 4 #if (object$glmm) 7 else 4
#     
#     n.trait <- (nc - 1)/nn.per.tr
#     if (n.trait != round(n.trait)) 
#       stop("something wrong")
#     fit$residuals <- resi[, 1 + (((nn.per.tr - 2) * n.trait + 
#                                        1):((nn.per.tr - 1) * n.trait))]
#     fit$fitted <- resi[, 1 + ((1 * n.trait + 1):(2 * n.trait))]
#     #fit$hat.matrix <- resi[, 1 + (((nn.per.tr - 1) * n.trait + 
#     #                                    1):(nn.per.tr * n.trait))]
#     
#     return(fit)
#   }


computeG <- function(W=NULL,pdf=TRUE) {
  SS <- tcrossprod(W)        # compute crossproduct all SNPs
  N <- tcrossprod(!W==0)     # compute number of observation all SNPs
  G <- SS/N
  if (pdf) G <- makepdf(G)
  return(list(G=G,SS=SS,N=N))
}  



  ginv <- function(G=NULL, tol=NULL) {
    rn <- rownames(G)
    e <- eigen(G)                                    # eigen value decomposition of the matrix G
    U <- e$vectors                                   # eigen vectors
    e <- e$values                                    # eigen values
    ie <- e
    ie[e>tol] <- 1/e[e>tol]
    ie[e<tol] <- 0
    D <- diag(ie)                                   # set inverse D to 1/e
    G <- U%*%D%*%t(U)   			                        # compute inverse
    ldet <- sum(log(e[e>tol])) 
    colnames(G) <- rownames(G) <- rn
    return(list(G=G,ldet=ldet))                      # log determinant 
  } 
  
  vec2mat <- function(vec=NULL,n=NULL, rcnames=NULL) {
    X <- diag(n)
    X[lower.tri(X, diag=TRUE)] <- vec
    X <- X + t(X) - diag(diag(X))  
    if(!is.null(rcnames)) { rownames(X) <- colnames(X) <- rcnames}
    X
  }
  
  writeG <- function(G=NULL, filename=NULL, clear=TRUE, ldet=NULL) {
    if (clear) {file.remove(filename)}
    nr <- nrow(G) 
    if(!is.null(ldet)) { write.table( t(c(0,0,ldet )),filename,quote=F,sep=" ",row.names=F, col.names=F, append=TRUE)}
    for (i in 1:nr) {
      out <- data.frame(rownames(G)[i], rownames(G)[i:nr], G[i:nr,i])
      write.table( out,filename,quote=F,sep=" ",row.names=F, col.names=F, append=TRUE)
    }
  }

