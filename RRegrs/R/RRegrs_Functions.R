# ======================================================================
# RRegrs - R Regressions
# ======================================================================
# Get the best regression models for one dataset using R caret methods
# eNanoMapper.net
# -------------------------------------------------------------------------------------------------------------
# AUTHORS: 
# -------------------------------------------------------------------------------------------------------------
# Georgia Tsiliki: ChemEng - NTUA, Greece, g_tsiliki@hotmail.com
# Cristian R. Munteanu: RNASA-IMEDIR, University of A Coruna, Spain, muntisa@gmail.com
# Jose A. Seoane: Stanford Cancer Institute, USA, seoane@stanford.edu
# Carlos Fernandez-Lozano: RNASA-IMEDIR, University of A Coruna, Spain, carlos.fernandez@udc.es
# Haralambos Sarimveis: ChemEng - NTUA, Greece, hsarimv@central.ntua.gr
# Egon Willighagen: BiGCaT - Maastricht University, Netherlands, egon.willighagen@gmail.com
# -------------------------------------------------------------------------------------------------------------

#library(caret)
#======================================================================================================================
# General functions
#======================================================================================================================
r2.adj.t.funct<- function(obs,pred,num.pred){
#obs==y, pred=predicted, num.pred=number of idependent variables (predictors)
#t: traditional formula
  y.mean<- mean(obs)
  x.in<- sum((obs-pred)^2)/sum((obs-y.mean)^2)
  x.in<- 1-x.in #r squared
  
  x.in<- (1-x.in)*((length(obs)-1)/(length(obs)-num.pred-1))
  x.in<- 1 - x.in 
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------
r2.adj.funct<- function(obs,pred,num.pred){
#obs==y, pred=predicted, num.pred=number of idependent variables (predictors)
 
  x.in<- cor(obs,pred)^2
  x.in<- (1-x.in)*((length(obs)-1)/(length(obs)-num.pred-1))
  x.in<- 1 - x.in 
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------
rmse.funct<- function(obs,pred){ 
#obs==y, pred=predicted
  return(sqrt(mean((pred - obs)^2)))
}
#----------------------------------------------------------------------------------------------------------------------
r2.funct<- function(obs,pred){  
#obs==y, pred=predicted
  x.in<- cor(obs,pred)^2
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------
r2.t.funct<- function(obs,pred){  
#obs==y, pred=predicted
  y.mean<- mean(obs)
  x.in<- sum((obs-pred)^2)/sum((obs-y.mean)^2)
  x.in<- 1-x.in #r squared
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------

AppendList2CSv <- function(l,csvFile) {
  #--------------------------------------------------------------------
  # Write a LIST to CSV file
  #--------------------------------------------------------------------
  out_file <- file(csvFile, open="a")  #creates a file in append mode 
  for (i in seq_along(l)){ 
    # writes the name of the list elements ("A", "B", etc.)
    write.table(names(l)[i],file=out_file,sep=",",dec=".",quote=F,col.names=F,row.names=F) 
    write.table(l[[i]],     file=out_file,sep=",",dec=".",quote=F,col.names=NA,row.names=T) #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  
#----------------------------------------------------------------------------------------------------------------------

AppendList2txt <- function(l,csvFile) {
  #--------------------------------------------------------------------
  # Write a LIST to TXT file
  #--------------------------------------------------------------------
  out_file <- file(csvFile, open="a")  #creates a file in append mode 
  for (i in seq_along(l)){ 
    #writes the name of the list elements ("A", "B", etc) 
    write.table(names(l)[i],file=out_file,sep=" ",dec=".",quote=F,col.names=F, row.names=F)
    write.table(l[[i]],     file=out_file,sep=" ",dec=".",quote=F,col.names=NA,row.names=T) #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  

# ************************************
# RRegrs Specific functions
# ************************************

RemNear0VarCols <- function(ds,fDet=FALSE,outFile="ds3.No0Var.csv") {
  #================================================
  # Removal of near zero variance columns (Step 3)
  #================================================
  # inputs:
  # - ds = dataset frame
  # - fDet = flag for detais (TRUE/FALSE)
  # - outFileName = new file name  (it could include the path)
  
  # output = ds.Rem0NearVar  (ds without columns with near zero variance)
  # if datails = TRUE, output the new ds as a file
  # ------------------------------------------
  
  # default parameters are no details, with a CSV file name
  #library(caret)
  ds.Rem0NearVar <- ds               # default output without any modification
  ds.var <- nearZeroVar(ds)          # get the near zero columns
  if (!length(ds.var) == FALSE) {    # remove the columns only if nearZeroVar identified; if no columns to remove, ds will be the same
    ds.Rem0NearVar <- ds[,-(ds.var)] # get only the columns without this problem
    if (fDet == TRUE) {              # write as details the corrected ds file
      write.csv(ds.Rem0NearVar, outFile,row.names=F, quote=F)
    }
  }
  return(as.data.frame(ds.Rem0NearVar)) # return the new data frame without near zero variance
}
#----------------------------------------------------------------------------------------------------------------------

ScalingDS <- function(ds,s=1,c=1,fDet=FALSE,outFileName="ds4.scaled.csv") {
  #===========================
  # Scaling dataset (Step 4)
  #===========================
  # s = { 1,2,3 } - type of scaling: 1 = normalization, 2 = standardization, 3 = other
  # c = the number of column into the dataset to start scaling
  # - if c = 1: included the dependent variable
  # - if c = 2: only the features will be scaled
  # fDet = if details need to be printed    (TRUE/FALSE)
  # outFileName = new file name  (it could include the path)
  # Default scaling = NORMALIZATION !
  
  # DEFAULT scaled dataset = original
  # if other s diffent of 1,2,3 is used => no scaling!
  DataSet.scaled <- ds
  
  # if NORMALIZATION
  if (s==1) {
    # Scale all the features (from column c; column 1 is the predictor output)
    DataSet.scaled <- ((ds-min(ds))/(max(ds)-min(ds))) # normalize all the columns
  }
  # if STADARDIZATION
  if (s==2) {
    # Scale all the features (from column c; column 1 is the predictor output)
    DataSet.scaled <- scale(ds[c:ncol(ds)],center=TRUE,scale=TRUE)  
  }
  
  # if other scaling
  if (s==3) {
    # Scale all the features (from feature 2 bacause feature 1 is the predictor output)
    # TO ADD THE CODE ! 
  }
  
  # if DETAILS
  if (fDet ==TRUE) {
    # write the result into a separated file
    write.csv(DataSet.scaled, outFileName,row.names=F, quote=F)  
  }
  return (as.data.frame(DataSet.scaled)) # return the scaled data frame
}
#----------------------------------------------------------------------------------------------------------------------

RemCorrs <- function(ds,fDet,cutoff,outFile) {
  # ========================================
  # Remove the correlated columns (Step 5)
  # ========================================
  # ds = dataset frame
  # fDet = flag fro details (TRUE/FALSE)
  # cutoff = correlation cut off (ex: 0.9)
  # outFileName = new file name  (it could include the path)
  
  # Generates 5 file results:
  # - returns a dataset without the correlated columns  (1 file)
  # - generate initial correlation matrix 
  #   and the one after removing the correlated features (2 files)
  # - plots for the before and after correlation removal (2 files)
  # ------------------------------------------------------------------------
  
  # another version of this function should be implemented using
  # pairwise test between i and j descriptors- if(r2>=0.9){remove the j descriptor}
  # using findCorelations() from caret
  
  #library(corrplot) #corrplot: the library to compute correlation matrix.
  #library(caret)
  
  DataSet <- ds               # input dataset
  DataSetFiltered.scale <- ds # default results without any modification
  
  # calculate the correlation matrix for the entire file!
  # !!! NEED TO BE CORRECTED to avoid dependent variable (first column) but to report it!
  corrMat <- cor(DataSet)                                             # get corralation matrix
  
  if (fDet==TRUE) {
    CorrMatFile <- paste(outFile,".corrMAT.csv",sep='')
    # write correlation matrix as output file
    write.csv(corrMat, CorrMatFile, row.names=F, quote=F)
    
    # Plot the matrix, clustering features by correlation index
    # corrplot(corrMat, order = "hclust")
    
    # plot the correlatio plot before correlation removal
    CorrPlotFile <-  paste(outFile,".corrs.png",sep='')
    png(height=1200, width=1200, pointsize=25, filename=CorrPlotFile)
    col1 <-rainbow(100, s = 1, v = 1, start = 0, end = 0.9, alpha = 1)
    corrplot(corrMat,tl.cex=3,title="Initial feature correlation matrix",
             method="circle",is.corr=FALSE,#type="full",
             cl.lim=c(-1,1),cl.cex=2,addgrid.col="red",
             addshade="positive",col=col1,
             addCoef.col = rgb(0,0,0, alpha = 0.6), mar=c(0,0,1,0), diag= FALSE) 
    dev.off()
  }
  
  highlyCor <- findCorrelation(corrMat, cutoff) # find corralated columns
  
  # if no correlation found, return the original dataset
  if (length(highlyCor) == 0){
    return (ds)
  }

  # Apply correlation filter with the cutoff only if exists!
  # by removing all the variable correlated with more than cutoff
  DataSetFiltered.scale <- DataSet[,-highlyCor]
  
  if (fDet==TRUE) {
    corrMat <- cor(DataSetFiltered.scale)
    # plot again the rest of correlations after removing the correlated columns
    #corrplot(corrMat, order = "hclust")
    
    # plot the correlation plot AFTER correlation removal
    #CorrPlotFile2 =  paste(outFile,".afterRemCorr.png",sep='')
    #png(height=1200, width=1200, pointsize=25, file=CorrPlotFile2)
    #col1 <-rainbow(100, s = 1, v = 1, start = 0, end = 0.9, alpha = 1)
    #corrplot(corrMat,tl.cex=3,title="Correlation matrix after removing correlated features",
    #         method="circle",is.corr=FALSE,type="full",
    #         cl.lim=c(-1,1),cl.cex=2,addgrid.col="red",
    #         addshade="positive",col=col1,
    #         addCoef.col = rgb(0,0,0, alpha = 0.6), mar=c(0,0,1,0), diag= FALSE) 
    #dev.off()
    # correlation matrix for the rest of the columns after removal
    #CorrMatFile2 <- paste(outFile,".corrMAT4Selected.csv",sep='')
    # write correlation matrix as output file
    #write.csv(corrMat, CorrMatFile2, row.names=F, quote=F)
    # write the new dataset without the correlated features
    write.csv(DataSetFiltered.scale, outFile, row.names=F, quote=F)
  }
  
  return(as.data.frame(DataSetFiltered.scale))
}
#----------------------------------------------------------------------------------------------------------------------

DsSplit <- function(ds,trainFrac=3/4,fDet=FALSE,PathDataSet="",iSeed) {
  # ===============================================
  # Dataset spliting in Training and Test (Step 6)
  # ===============================================
  # Inputs
  # - ds = frame dataset object
  # - fDet = flag for detais (TRUE/FALSE)
  # - PathDataSet = pathway for results
  # Output = training and test datasets (to be used for regressions in other functions)
  
  # if datails = TRUE, output files will be created
  my.datf<- ds
  
  # create TRAIN and TEST sets to build a model
  set.seed(iSeed)
  inTrain <- createDataPartition(1:dim(my.datf)[1],p = trainFrac,list = FALSE,groups=2)
  # groups==2 forces to NOT partition
  # based on quantiles of numeric values
  
  my.datf.train<- my.datf[inTrain,]            # TRAIN dataset frame         
  my.datf.test <- my.datf[-inTrain,]           # TEST dataset frame
  
  if (fDet == TRUE) {
    # write the TRAIN and TEST set files
    # the index of each row will in the dataset will not be saved (row.names=F)
    outTrain <- file.path(PathDataSet,paste("ds.Train.split",iSeed,".csv")) # the same folder as the input
    write.csv(my.datf.train,outTrain,row.names=FALSE)
    outTest <- file.path(PathDataSet,paste("ds.Test.split",iSeed,".csv")) # the same folder as the input
    write.csv(my.datf.test,outTest,row.names=FALSE) 
  }
  MyList<- list("train"=my.datf.train, "test"=my.datf.test) 
  return(MyList)  # return a list with training and test datasets
}

# *************************************
# REGRESSION METHODS
# *************************************

LMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #==================
  # 8.1. Basic LM
  #==================

  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "lm" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=10,repeats=10,
                      summaryFunction=defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  lm.fit<- train(net.c~.,data=my.datf.train,
                 method='lm', tuneLength = 10,trControl=ctrl,
                 metric='RMSE')
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- lm.fit$results[,2]
  R2.tr    <- lm.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- lm.fit$results[,4]
    R2sd.tr   <- lm.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(lm.fit)
  lm.test.res  <- postResample(predict(lm.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(lm.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(lm.fit,my.datf.test)  # predicted Y for test
  noFeats.fit <- length(predictors(lm.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(lm.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(lm.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # -------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   "RMSE.tr"   = as.numeric(RMSE.tr),
                   "R2.tr"     = as.numeric(R2.tr),
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr),
                   "R2sd.tr"   = as.numeric(R2sd.tr),
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                   "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile, append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),       file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",         file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),        file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    write.table("Fitting Summary: ",                      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(summary(lm.fit)$coefficients), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res),file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(lm.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- lm.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
  
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Leverage / Hat values
    hat.fit <- hatvalues(fitModel)          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    # Cook's distance
    cook.dists<- cooks.distance(fitModel)
    cutoff.Cook <- 4/((nrow(my.datf.train)-length(fitModel$coefficients)-2)) # Cook's distance cutoff
    
    write.table("Cook's distances output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Cook's distance cutoff: ", cutoff.Cook), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Cook's distances: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(cook.dists), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Influence
    infl <- influence(fitModel)#produces several statistics of the kind
    
    write.table("Point influence output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Influences: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(infl), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    if(length(is.na(c(FeatImp$importance$Overall)))<=(length(c(FeatImp$importance$Overall))-3)){  
     dotchart(as.matrix(FeatImp$importance),main="Feature Importance")}                          # plot 3    
  
    # Fitted vs Residuals - plot 4
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    # Cook's distance - plot 6
    if(length(is.na(cook.dists))<=(length(cook.dists)-3)){  
    	   plot(cook.dists,
         main="Cook's Distance for Fitted Model",
         xlab="Index", ylab="Cook Distance")
    
         for (p in 1:6) {
           plot(fitModel, which=p, cook.levels=cutoff.Cook) # 6 standard fitting plots
         }
    }    
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------

  }
  return(list(stat.values=my.stats, model=lm.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

GLMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #======================================================
  # 8.2- GLM stepwise regression - based on AIC (caret)
  #======================================================
  # Inputs:
  # - my.datf.train,my.datf.test = training and test dataset frames
  # - sCV = type of cross-validation such as repeatedcv, LOOCV, etc.
  # - iSplit = index of splitalse
  # - fDet = flag for detais (True/F)
  # - outFile = output file for GLM details
  # Output:
  # - list of statistics equal with the header introduced in the main script and the full model
  #   (tr = train, ts = test, both = tr+ts = full dataset)
  # -----------------------------------------------------------------------------------------------
  
  #library(caret)
  #attach(my.datf.train)    # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "glmStepAIC" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=10,repeats=10,
                      summaryFunction=defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  glm.fit<- train(net.c~.,data=my.datf.train,
                  method='glmStepAIC', tuneLength=10, trControl=ctrl,
                  metric='RMSE')
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- glm.fit$results[,2]
  R2.tr    <- glm.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- glm.fit$results[,4]
    R2sd.tr   <- glm.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(glm.fit)
  lm.test.res  <- postResample(predict(glm.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(glm.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(glm.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(glm.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(glm.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(glm.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # -------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   "RMSE.tr"   = as.numeric(RMSE.tr),
                   "R2.tr"     = as.numeric(R2.tr),
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr),
                   "R2sd.tr"   = as.numeric(R2sd.tr),
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                   "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper",                  file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    write.table("Fitting Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(summary(glm.fit)$coefficients), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(glm.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(glm.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- glm.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Leverage / Hat values
    hat.fit <- hatvalues(fitModel)          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
    # Cook's distance
    cook.dists<- cooks.distance(fitModel)
    cutoff.Cook <- 4/((nrow(my.datf.train)-length(fitModel$coefficients)-2)) # Cook's distance cutoff
    
    write.table("Cook's distances output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Cook's distance cutoff: ", cutoff.Cook), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Cook's distances: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(cook.dists), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Influence
    infl <- influence(fitModel)#produces several statistics of the kind
    
    write.table("Point influence output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Influences: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(infl), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    if(length(is.na(c(FeatImp$importance$Overall)))<=(length(c(FeatImp$importance$Overall))-3)){  
     dotchart(as.matrix(FeatImp$importance),main="Feature Importance")}                          # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    # Cook's distance - plot 6
    if(length(is.na(cook.dists))<=(length(cook.dists)-3)){  
    	   plot(cook.dists,
         main="Cook's Distance for Fitted Model",
         xlab="Index", ylab="Cook Distance")
    
         for (p in 1:6) {
           plot(fitModel, which=p, cook.levels=cutoff.Cook) # 6 standard fitting plots
         }
    }    
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
  }
  return(list(stat.values= my.stats, model=glm.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

PLSreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #================================
  # 8.3. PLS regression (caret)
  #================================
  
  #library(caret)
  
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "pls" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  floor.param<- floor((dim(my.datf.train)[2]-1)/5)
  if(floor.param<1){floor.param <- 1}  

  pls.fit<- train(net.c~.,data=my.datf.train,
                  method = 'pls', tuneLength = 10, trControl = ctrl,
                  metric = 'RMSE',
                  tuneGrid=expand.grid(.ncomp=c(1:floor.param)))

  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- pls.fit$results[,2]
  R2.tr    <- pls.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- pls.fit$results[,4]
    R2sd.tr   <- pls.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(pls.fit)
  lm.test.res  <- postResample(predict(pls.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(pls.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(pls.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(pls.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(pls.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(pls.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # --------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
                   
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                   "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",   file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),  file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(pls.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(pls.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- pls.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # ADDED !
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)

    if(det.Traind.pls!=0){
     Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
     Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
     
     # Leverage / Hat values
     hat.fit <- Hat.test          # hat values
     hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
     hat.mean <- mean(hat.fit)               # mean hat values
     hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
     
     write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
     
     #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
     thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
     hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
     
     write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
     
     # Cook's distance ?
    }

    # Influence ?
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    if(det.Traind.pls!=0){
     plot(hat.fit, type = "h",
          main="Leverage for Fitted Model",
          xlab="Index", ylab="Hat")
     abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }
    dev.off()
    # --------------------------------------------------------------
  }
  
  return(list(stat.values=my.stats, model=pls.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

LASSOreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #================================
  # 8.4 Lasso Regression (caret)
  #================================
  
  #library(caret)
  
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "lasso.RMSE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  las.fit<- train(net.c~.,data=my.datf.train,
                  method='lasso', tuneLength = 10, trControl = ctrl,
                  metric='RMSE' ) #,tuneGrid=expand.grid(.fraction= seq(0.1,1,by=0.1)))
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- las.fit$results[,2]
  R2.tr    <- las.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- las.fit$results[,4]
    R2sd.tr   <- las.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(las.fit)
  lm.test.res  <- postResample(predict(las.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(las.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(las.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(las.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(las.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(las.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # --------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
           
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                   "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(las.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(predictors(lm.train.res),outFile)
    #write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,quote=F)
    AppendList2CSv(predictors(lm.test.res),outFile)
    #write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(las.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- las.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- abs(pred.both - ds.full[,1]) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)

    if(det.Traind.pls!=0){
     Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
     Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
     
     # Leverage / Hat values
     hat.fit <- Hat.test          # hat values
     hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
     hat.mean <- mean(hat.fit)               # mean hat values
     hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
     
     write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
     
     #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
     thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
     hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
     
     write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    
     # Cook's distance ?
    }
    # Influence ?
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(pred.both,resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    if(det.Traind.pls!=0){
     plot(hat.fit, type = "h",
          main="Leverage for Fitted Model",
          xlab="Index", ylab="Hat")
     abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }
    dev.off()
    # --------------------------------------------------------------
  }
  return(list(stat.values=my.stats, model=las.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

# RBF_DDAreg <- function(my.datf.train,my.datf.test,sCV,negThrStep=0.5,iSplit=1,fDet=F,outFile="") {
#   #============================================================
#   # 8.5. RBF network with the DDA algorithm regression (caret)
#   #============================================================
#   
#   library(caret)
#   
#   net.c = my.datf.train[,1] # dependent variable is the first column in Training set
#   RegrMethod <- "rbfDDA" # type of regression
#   
#   # Define the CV conditions
#   ctrl<- trainControl(method=sCV,number=10,repeats=10,
#                       summaryFunction=defaultSummary)
#   
#   # Train the model using only training set
#   set.seed(iSplit)
#   
#   rbf.fit<- train(net.c~.,data=my.datf.train,
#                   method='rbfDDA',trControl=ctrl,
#                   tuneGrid=expand.grid(.negativeThreshold=seq(0,1,negThrStep)))
#   
#   #------------------------------
#   # Training RESULTS
#   #------------------------------
#   RMSE.tr  <- rbf.fit$results[,2]
#   R2.tr    <- rbf.fit$results[,3]
#   if (sCV == "repeatedcv"){ # if 10-fold CV
#     RMSEsd.tr <- rbf.fit$results[,4]
#     R2sd.tr   <- rbf.fit$results[,5]
#   }
#   if (sCV == "LOOCV"){ # if LOOCV
#     RMSEsd.tr <- 0 # formulas will be added later!
#     R2sd.tr   <- 0 # formulas will be added later!
#   }
#   
#   #------------------------------------------------
#   # RMSE & R^2, for train/test respectively
#   #------------------------------------------------
#   lm.train.res <- getTrainPerf(rbf.fit)
#   lm.test.res  <- postResample(predict(rbf.fit,my.datf.test),my.datf.test[,1])
#   
#   #------------------------------------------------
#   # Adj R2, Pearson correlation
#   #------------------------------------------------
#   pred.tr     <- predict(rbf.fit,my.datf.train) # predicted Y
#   pred.ts     <- predict(rbf.fit,my.datf.test)  # predicted Y
#   noFeats.fit <- length(predictors(rbf.fit))    # no. of features from the fitted model
#   Feats.fit   <- paste(predictors(rbf.fit),collapse="+") # string with the features included in the fitted model
#   
#   ds.full     <- rbind(my.datf.train,my.datf.test)
#   pred.both   <- predict(rbf.fit,ds.full)       # predicted Y
#   adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
#   adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
#   corP.ts     <- cor(my.datf.test[,1],pred.ts)
#   
#   adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
#   RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
#   r2.both     <- r2.funct(ds.full[,1],pred.both)
#   
#   # Generate the output list with statistics for each cross-validation type
#   # --------------------------------------------------------------------
#   my.stats <- list("RegrMeth"     = RegrMethod,
#                    "Split No"     = as.numeric(iSplit),     # from function param
#                    "CVtype"       = sCV,                    # from function param
#                    "NoModelFeats" = as.numeric(noFeats.fit),
#                    "ModelFeats"   = Feats.fit,
#                    "adjR2.tr"  = as.numeric(adjR2.tr),
#                    
#                    "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
#                    "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
#                    "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
#                    "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
#                    
#                    "adjR2.ts"= as.numeric(adjR2.ts),
#                    "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
#                    "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
#                    "corP.ts" = as.numeric(corP.ts),
#                    "adjR2.both" = as.numeric(adjR2.both),
#                    "RMSE.both"  = as.numeric(RMSE.both),
#                    "R2.both"    = as.numeric(r2.both))
#   #---------------------------------------------------------------------
#   # Write to file DETAILS for GLM for each cross-validation method
#   #---------------------------------------------------------------------
#   if (fDet==T) {   # if flag for details if true, print details about any resut
#     write("RRegr package | eNanoMapper", file=outFile,append=T)
#     write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table(paste("Split no.: ", iSplit),             file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table(paste("CV type: ", sCV),  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table(summary(my.datf.train), file=outFile,append=T,sep=",",col.names=T,quote=F)
#     write.table("Test Set Summary: ",   file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table(summary(my.datf.test),  file=outFile,append=T,sep=",",col.names=T,quote=F)   
#     
#     
#     write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table(predictors(rbf.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
#     
#     write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
#     write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
#     
#     write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,row.names=F,quote=F)
#     
#     # Variable Importance (max top 20)
#     FeatImp <- varImp(rbf.fit, scale = F)
#     components = length(FeatImp)  # default plot all feature importance
#     if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
#       components = 20
#     }
#     # Append feature importance to output details
#     AppendList2CSv(FeatImp,outFile)
#     
#     fitModel <- rbf.fit$finalModel
#     
#     # =============================================================================
#     # Assessment of Applicability Domain (plot leverage)
#     # =============================================================================
#     
#     # Residuals
#     resids <- residuals(fitModel) # residuals
#     write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#     write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
#     
#     predVals.pls.ad <- pred.ts
#     Traind.pls= as.matrix(my.datf.train)
#     Testd.pls = as.matrix(my.datf.test)
#     mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
#     det.Traind.pls<- det(mat.Traind.pls)
# 
#     if(det.Traind.pls!=0){
#   
#      Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
#      Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
#      
#      # Leverage / Hat values
#      hat.fit <- Hat.test          # hat values
#      hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
#      hat.mean <- mean(hat.fit)               # mean hat values
#      hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
#      
#      write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#      write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#      write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#      write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
#      
#      #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
#      thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
#      hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
#      
#      write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#      write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
#      write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
#      
#      # Cook's distance ?
#     }
# 
#     # Influence ?
#     
#     # PDF plots
#     # --------------------------------------------------------------
#     pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
#     plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
#     plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
#     dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
#     
#     # Fitted vs Residuals
#     plot(fitted(fitModel),residuals(fitModel),
#          main="Fitted vs. Residuals for Fitted Model",
#          xlab="Fitted", ylab="Residuals")
#     abline(h = 0, lty = 2)
#     
#     # Leverage plots
#     if(det.Traind.pls!=0){
#      plot(hat.fit, type = "h",
#           main="Leverage for Fitted Model",
#           xlab="Index", ylab="Hat")
#      abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
#     }
#     dev.off()
#     # --------------------------------------------------------------
#   }
#   return(list(stat.values=my.stats, model=rbf.fit))  # return a list with statistics and the full model
# }
#----------------------------------------------------------------------------------------------------------------------

SVRMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="",cs=c(1,5,10,15,20)) {
  #====================================
  # 8.6 SVM Radial Regression (caret)
  #====================================
  
  #library(caret)
  #library(kernlab)

  cs = as.numeric(cs)

  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "svmRadial" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV,number=10,repeats=10,
                      summaryFunction=defaultSummary)

  # Train the model using only training set
  set.seed(iSplit)
  sigma = sigest (as.matrix(my.datf.train[,-1]))[2]
  svmL.fit<- train(net.c~.,data=my.datf.train,
                   method='svmRadial',tuneLength=10,trControl=ctrl,
                   metric='RMSE',
                   tuneGrid=expand.grid(.sigma=sigma,.C= cs))
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- svmL.fit$results[,which(names(svmL.fit$results)=='RMSE')]#2]
  R2.tr    <- svmL.fit$results[,which(names(svmL.fit$results)=='Rsquared')]#3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- svmL.fit$results[,which(names(svmL.fit$results)=='RMSESD')]##4]
    R2sd.tr   <- svmL.fit$results[,which(names(svmL.fit$results)=='RsquaredSD')]#5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(svmL.fit)
  lm.test.res  <- postResample(predict(svmL.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(svmL.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(svmL.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(svmL.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(svmL.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(svmL.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # ------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
                   
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                   "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",       file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(svmL.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(svmL.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- svmL.fit$finalModel
       
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals    
    resids = pred.tr-svmL.fit$trainingData$.outcome # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)

    if(det.Traind.pls!=0){
     Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
     Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
     
     # Leverage / Hat values
     hat.fit <- Hat.test          # hat values
     hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
     hat.mean <- mean(hat.fit)               # mean hat values
     hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
     
     write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
     
     #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
     thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
     hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
     
     write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
     
     # Cook's distance ?
    }

    # Influence ?
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(fitModel),resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    if(det.Traind.pls!=0){
     plot(hat.fit, type = "h",
          main="Leverage for Fitted Model",
          xlab="Index", ylab="Hat")
     abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }

    dev.off()
    # --------------------------------------------------------------
  }

  return(list(stat.values=my.stats, model=svmL.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

NNreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #========================================
  # 8.8 Neural Network Regression (caret)
  #========================================

  #library(caret)

  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "nnet" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,
                      summaryFunction = defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  
  nn.fit<- train(net.c~.,data=my.datf.train,
                 method = 'nnet',trControl = ctrl,
                 linout=T, trace = F,MaxNWts=20000,
                 #Grid of tuning parameters to try:
                 tuneGrid=expand.grid(.size=c(1,5,10,15),.decay=c(0,0.001,0.1)))
  #Grid parameters are appearing at the print out of the model
  #size==#of units in hidden layer, decay==parameter of weight decay (default:0)
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- nn.fit$results[,which(names(nn.fit$results)=='RMSE')]#2]
  R2.tr    <- nn.fit$results[,which(names(nn.fit$results)=='Rsquared')]#3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- nn.fit$results[,which(names(nn.fit$results)=='RMSESD')]#4]
    R2sd.tr   <- nn.fit$results[,which(names(nn.fit$results)=='RsquaredSD')]#5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(nn.fit)
  lm.test.res  <- postResample(predict(nn.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(nn.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(nn.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(nn.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(nn.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(nn.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # ----------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
                   
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                   "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper",                  file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(nn.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(nn.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- nn.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)

    if(det.Traind.pls!=0){
     Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
     Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
     
     # Leverage / Hat values
     hat.fit <- Hat.test          # hat values
     hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
     hat.mean <- mean(hat.fit)               # mean hat values
     hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
     
     write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
     
     #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
     thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
     hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
     
     write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
     
     # Cook's distance ?
    }
    # Influence ?
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    if(det.Traind.pls!=0){
     plot(hat.fit, type = "h",
          main="Leverage for Fitted Model",
          xlab="Index", ylab="Hat")
     abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }

    dev.off()
    # --------------------------------------------------------------
  }

  return(list(stat.values=my.stats, model=nn.fit))  # return a list with statistics and the full model
}


# **************************************
# WRAPPER METHODS
# **************************************

PLSregWSel <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #====================================================================================================
  # 8.3W. PLS regression with filter feature selection (caret)
  #====================================================================================================

  #library(caret)
  
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "pls.WSel" # type of regression
  
  # Define the CV conditions
  ctrlw <- rfeControl(method = 'boot', number = 25,saveDetails=T)#number=10,repeats=10
  ctrl  <- trainControl(method = sCV, number = 10,repeats = 1,#numebr=10,repeats=10,
                        summaryFunction = defaultSummary,savePredictions=T)
  
  subsetsx<- seq(2,dim(my.datf.train)[2]-1, by = 10)
  
  # Train the model using only training set
  set.seed(iSplit)
  floor.param<- floor((dim(my.datf.train)[2]-1)/5)
  if(floor.param<1){floor.param <- 1}   
  pls.fit<- rfe(net.c~.,data=my.datf.train, 
                method = 'pls',
                rfeControl = ctrlw, trControl=ctrl, sizes=subsetsx, importance=T,
                metric = 'RMSE',
                tuneGrid=expand.grid(.ncomp=c(1:floor.param)))#ncomp=c(1:5)
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  
  pls.fit.best <- subset(pls.fit$results, pls.fit$results$Variables == 5) # best selected fit
  
  RMSE.tr  <- pls.fit.best$RMSE
  R2.tr    <- pls.fit.best$Rsquared
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- pls.fit.best$RMSESD
    R2sd.tr   <- pls.fit.best$RsquaredSD
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- pls.fit # ??
  lm.test.res  <- postResample(predict(pls.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(pls.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(pls.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(pls.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(pls.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(pls.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # -------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"     = as.numeric(adjR2.tr),
                   
                   "RMSE.tr"   = as.numeric(RMSE.tr),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr),
                   "R2sd.tr"   = as.numeric(R2sd.tr),
                   
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                   "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if True, print details about any resut
    write("RRegr package | eNanoMapper",                  file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(pls.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("NNet variable importance: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2txt(varImp(pls.fit),outFile)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(pls.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    # PDF plots
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    par(mfrow = c(2, 2))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    #plot(FeatImp, top = components,main="Feature Importance")
    dev.off()
  }
  return(list(stat.values=my.stats, model=pls.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

Yrandom<- function(dss,trainFrac,best.reg,best.R2.ts,noYrand,ResBestF,rfe_SVM_param_c,rfe_SVM_param_eps){
  #================================================
  # Y-randomization for the best model (Step 12)
  #================================================
  #    - 1 splitting, 1 CV type, best method
  #    - best.R2.ts will be compared with Yrand.R2.ts
  #    - returns ratios DiffsR2/bestR2
  #   (ex param for rbfDDa: negThrStep)
  # --------------------------------------------------
  
  cat("-> Best model Y-Randomization ...\n")
  dss[,1] <- sample(dss[,1]) # randomize Y values for the entire dataset
  
  # splitting dataset in training and test
  #---------------------------------------
  Yrand.R2.ts <- NULL     # all values of R2 for each Y randomization
  for (i in 1:noYrand){
    iSeed=i               
    dsList  <- DsSplit(dss,trainFrac,F,PathDataSet,iSeed) # return a list with 2 datasets = dsList$train, dsList$test
    # get train and test from the resulted list
    ds.train<- dsList$train
    ds.test <- dsList$test
    
    # Run the caret function with the method from the best method
    #    for one training-test split only; no details, we need only R2 values
    if (best.reg=="lm") {
      my.stats.reg  <- LMreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run GLM for each CV and regr method
    }
    if (best.reg=="glmStepAIC") {
      my.stats.reg  <- GLMreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run GLM for each CV and regr method
    }
    if (best.reg=="pls") {
      my.stats.reg  <- PLSreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run SVRM Radial for each CV and regr method
    }
    if (best.reg=="lasso.RMSE") {
      my.stats.reg  <- LASSOreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run SVRM Radial for each CV and regr method
    }
    if (best.reg=="svmRadial") {  
      my.stats.reg  <- SVRMreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF,rfe_SVM_param_c)$stat.values # run SVRM Radial for each CV and regr method
    }
    if (best.reg=="nnet") {  
      my.stats.reg  <- NNreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run NNet for each CV and regr method
    } 
    if (best.reg=="rf") {  
      my.stats.reg  <- RFreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run NNet for each CV and regr method
    } 
    if (best.reg=="svmRFE") {  
      my.stats.reg  <- SVMRFEreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF,rfe_SVM_param_c,rfe_SVM_param_eps)$stat.values # run NNet for each CV and regr method
    } 
    if (best.reg=="glmnet") {  
      my.stats.reg  <- ENETreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run NNet for each CV and regr method
    }  
    if (best.reg=="rfRFE") {  
      my.stats.reg  <- RFRFEreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run NNet for each CV and regr method
    } 
    
    #     if (best.reg=="rbfDDA") {  
    #       my.stats.reg  <- RBF_DDAreg(ds.train,ds.test,"repeatedcv",negThrStep,i,F,ResBestF)$stat.values # run SVRM Radial for each CV and regr method
    #    }

    Yrand.R2.ts <- c(Yrand.R2.ts,my.stats.reg$R2.ts) # adding test R2 value Y randomization
  }
  
  R2diffsPerBestR2 <- NULL

  if (is.na(my.stats.reg$R2.ts)) { # check for NA values 
    cat("       --> Y-Randomization error due to NA values!\n")
    write("Y-Randomization error due to NA values!", file=ResBestF,append=T)
    }
  else{
    # get histogram for differences between best R2 and the values for each Y randomization
    R2diffs          <- abs(Yrand.R2.ts - best.R2.ts) # absolute differences between R2 values (best model vs Y randomized results)
    R2diffsPerBestR2 <- abs(R2diffs/best.R2.ts)            # the same difference in percents
  
    pdf(file=paste(ResBestF,".Yrand.Hist.pdf",sep=""))    # save histogram if ratio diffs R2 into PDF for Y random
    Yrand.hist  <- hist(R2diffsPerBestR2)         # draw histogram of the ratio diffs/Best R2 for Y random
    dev.off()
  
    write.table("Y randomization test: ",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("=====================", file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Diffs R2 (Best Model - Y rand):",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(R2diffs, ResBestF)
    write.table("Summary Difs:",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(summary(R2diffs), ResBestF)
    write.table("Ratio Diffs R2 / Best R2 (Best Model - Y rand):",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(R2diffsPerBestR2, ResBestF)
    write.table("Summary Difs %:",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(summary(R2diffsPerBestR2),ResBestF)
  }
  return(R2diffsPerBestR2) # return the ratio of diffs with the best R2 from the same 
}

# -----------------------------------------------------------------------
# svm regression function helper
# -----------------------------------------------------------------------
# jseoane
# use:
# svmFuncsGradW: RAKOTOMAMONJY gradient w
data(model.svmRadialReg)

svmFuncsW = caretFuncs    ## regular ranking using w
svmFuncsW$fit=function(x,y,first,last,...,tuneGrid){
  #cat(param$sigma,"\n")
  #library(kernlab)
  sigma = sigest(x)[2]  
  cs = tuneGrid$.C
  eps = tuneGrid$.epsilon
  tuneGrid = expand.grid(.C=cs,.sigma=sigma,.epsilon=eps)  
  train(x,y,...,tuneGrid=tuneGrid)
}
#--------------------------------------------------------------------------------
svmFuncsW$rank=function(object,x,y){
  alphas = alpha(object$finalModel)
  alpha.idxs = alphaindex(object$finalModel)
  y.sv = as.numeric(y[alpha.idxs])
  w = (y.sv * alphas) %*% xmatrix(object$finalModel)
  sig = ifelse(object$finalModel@fitted>y,yes=1,no=-1)
  
  avImp = t(w*w)
  out = data.frame(avImp)
  colnames(out) = "Overall"
  out = out[order(out$Overall, decreasing = TRUE), , drop = FALSE]
  out$var <- rownames(out)
  out
}
#--------------------------------------------------------------------------------
svmFuncsW$pred= function(object, x)
{
  tmp = predict(object, newdata=x)
  if(object$modelType == "Classification" &
       !is.null(object$modelInfo$prob))
  {
    out1 =cbind(data.frame(pred = tmp),
                as.data.frame(predict(object$finalModel, newdata=x, type = "prob")))
  } else out1 <- tmp
  out1
}
#--------------------------------------------------------------------------------

# Based on the gradient of svm coefs
svmFuncsGradW = svmFuncsW
svmFuncsGradW$rank=function(object,x,y){ # RAKOTOMAMONJY gradient w  
  alphas = alpha(object$finalModel)#[[1]]
  alpha.idxs = alphaindex(object$finalModel)#[[1]]
  y.sv = y[alpha.idxs]  
  krnFun = kernelf(object$finalModel)
  kernel = kernelMatrix(krnFun,x)
  sigma = krnFun@kpar$sigma
  xmat = xmatrix(object$finalModel)[[1]]
  kerSV = kernel[alpha.idxs,alpha.idxs]  
  nSV = length(alpha.idxs)
  nfeat = dim(x)[2]
  avImp = numeric(nfeat)
  names(avImp) = colnames(x)
  
  for(i in 1:nfeat){
    deraux =  (  x[alpha.idxs,i] %*% t(as.matrix(rep(1,nSV))) ) -  (as.matrix(rep(1,nSV)) %*% t(x[alpha.idxs,i])     )    
    kernelDeriv1 = -(deraux * kerSV) / (sigma^2)
    kernelDeriv2 =  (deraux * kerSV) / (sigma^2)
    gradMarg1= -t(y.sv*alphas) %*% kernelDeriv1 %*%  (y.sv*alphas)
    gradMarg2= -t(y.sv*alphas) %*% kernelDeriv2 %*%  (y.sv*alphas)
    avImp[i] = gradMarg1^2 + gradMarg2^2
  }
  
  out = data.frame(avImp)
  colnames(out) = "Overall"
  out = out[order(out$Overall, decreasing = TRUE), , drop = FALSE]
  out$var <- rownames(out)
  out
}
#----------------------------------------------------------------------------------------------------------------------

RFreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #======================================
  # Basic RandomForest
  #======================================
  
  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "rf" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=10,repeats=10,#number=10,repeats=10,
                      summaryFunction=defaultSummary)
  
  tuneParam = data.frame(.mtry=c(ncol(my.datf.train)/3,ncol(my.datf.train)/2,ncol(my.datf.train)))
  # Train the model using only training set
  set.seed(iSplit)
  rf.fit<- train(net.c~.,data=my.datf.train,
                 method='rf', trControl=ctrl,
                 metric='RMSE',ntree=1500,tuneGrid =tuneParam)
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- rf.fit$results[rownames(rf.fit$bestTune),2]
  R2.tr    <- rf.fit$results[rownames(rf.fit$bestTune),3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- rf.fit$results[rownames(rf.fit$bestTune),4]
    R2sd.tr   <- rf.fit$results[rownames(rf.fit$bestTune),5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later! 
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  rf.train.res <- getTrainPerf(rf.fit)
  rf.test.res  <- postResample(predict(rf.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(rf.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(rf.fit,my.datf.test)  # predicted Y for test
  noFeats.fit <- length(predictors(rf.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(rf.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(rf.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # ------------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   "RMSE.tr"   = as.numeric(RMSE.tr),
                   "R2.tr"     = as.numeric(R2.tr),
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr),
                   "R2sd.tr"   = as.numeric(R2sd.tr),
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((rf.test.res["RMSE"])),
                   "R2.ts"   = as.numeric((rf.test.res["Rsquared"])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile, append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),       file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",         file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),        file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rf.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rf.train.res),file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rf.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- importance(rf.fit$finalModel, scale = T)
    FeatImp = FeatImp[order(FeatImp,decreasing=T),]
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- rf.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- fitModel$predicted-fitModel$y # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)

    if(det.Traind.pls!=0){
     Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
     Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
     
     # Leverage / Hat values
     hat.fit <- Hat.test          # hat values
     hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
     hat.mean <- mean(hat.fit)               # mean hat values
     hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
     
     write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
     
     #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
     thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
     hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
     
     write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    } 
     # PDF with 12 plots
     # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    dotchart(as.matrix(FeatImp),main="Feature Importance")                            # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(fitModel$predicted,resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    if(det.Traind.pls!=0){
     plot(hat.fit, type = "h",
          main="Leverage for Fitted Model",
          xlab="Index", ylab="Hat")
     abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
  }

  return(list(stat.values=my.stats, model=rf.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

SVMRFEreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="",cs=c(1,5,10,15,20),eps=c(0.01,0.1,0.3)) {
#SVMRFEreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="",cs=c(1:10),eps=c(0.01,0.1,0.3),noCores=1) {
  
  #===========================================
  # SVM-RFE
  #===========================================
  
  #library(kernlab)

  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "svmRFE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=3,repeats=1,#number=10,repeats=10,
                      summaryFunction=defaultSummary,verboseIter = F)
  
  rfeCtr = rfeControl(functions=svmFuncsGradW,method="cv",number=10,repeats=10, saveDetails = T, verbose=T,rerank = T,allowParallel=T)#number=10,repeats=10,
  sigma = sigest (as.matrix(my.datf.train[,-1]))[2]
  #cs = c(0.0001,0.1,1,5,15,50)  # EXTERNAL PARAMETERS!!!
  #cs = c(1,5,15,50)
  #eps=c(0.01,0.1,0.3)   # EXTERNAL PARAMETERS!!!
  sizes = 2^(1:sqrt(ncol(my.datf.train)-1))
  tuneVars = expand.grid(.C=cs, .sigma=sigma, .epsilon=eps)    
  
  # Train the model using only training set
  set.seed(iSplit)
  
  rfesvm.fit = rfe(as.matrix(my.datf.train[,-1]),net.c,sizes = sizes,rfeControl=rfeCtr,prob.model =F,method=svmRadialReg,tuneGrid = tuneVars,trControl=ctrl ,allowParallel=T)
  
  # warning in some of the parameters is a extreme value
  if(rfesvm.fit$fit$bestTune$C %in% cs[c(1,length(cs))])
    warning("Best fitted value of C=",rfesvm.fit$fit$bestTune$C," is a extreme value in your possible c values. You may want to reset your C paramenter options", call. = FALSE)
  if(rfesvm.fit$fit$bestTune$epsilon %in% eps[c(1,length(cs))])
    warning("Best fitted value of eps=",rfesvm.fit$fit$bestTune$epsilon," is a extreme value in your possible epsilon values. You may want to reset your eps paramenter options", call. = FALSE)
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,2]
  R2.tr    <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,4]
    R2sd.tr   <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!  TODOOOOOOOOOOOOOO
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  rfesvm.train.res <- getTrainPerf(rfesvm.fit$fit)
  rfesvm.test.res  <- postResample(predict(rfesvm.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(rfesvm.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(rfesvm.fit,my.datf.test)  # predicted Y for test
  noFeats.fit <- length(predictors(rfesvm.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(rfesvm.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(rfesvm.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # -----------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   "RMSE.tr"   = as.numeric(RMSE.tr),
                   "R2.tr"     = as.numeric(R2.tr),
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr),
                   "R2sd.tr"   = as.numeric(R2sd.tr),
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((rfesvm.test.res["RMSE"])),
                   "R2.ts"   = as.numeric((rfesvm.test.res["Rsquared"])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile, append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),       file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",         file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),        file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rfesvm.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.train.res,file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.test.res, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    
    FeatImp <- svmFuncsGradW$rank(rfesvm.fit$fit,as.matrix(ds.full[,-1]),ds.full[,1])
    FeatImp = FeatImp[order(FeatImp[,1],decreasing=T),]
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- rfesvm.fit$fit
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- pred.tr-my.datf.train[,1]    # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)

    if(det.Traind.pls!=0){
     Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
     Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
     
     # Leverage / Hat values
     hat.fit <- Hat.test          # hat values
     hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
     hat.mean <- mean(hat.fit)               # mean hat values
     hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
     
     write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
     
     #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
     thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
     hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
     
     write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    } 
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    fi = as.matrix(FeatImp[,1])
    rownames(fi)=FeatImp[,2]
    dotchart(fi,main="Feature Importance") # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(fitted(fitModel),resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    if(det.Traind.pls!=0){ 
     plot(hat.fit, type = "h",
          main="Leverage for Fitted Model",
          xlab="Index", ylab="Hat")
     abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
  }

  return(list(stat.values=my.stats, model=rfesvm.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

RFRFEreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #=============================
  # Random Forest-RFE
  #=============================
  
  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "rfRFE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=5,repeats=1,#number=10,repeats=10,
                      summaryFunction=defaultSummary,verboseIter = F)
  
  rfeCtr = rfeControl(functions = rfFuncs,method="cv",number=10,repeats=10, saveDetails = T, verbose=T,rerank = F,allowParallel=T)#number=10,repeats=10,
  
  sizes = 2^(1:sqrt(ncol(my.datf.train)-1))
  
  # Train the model using only training set
  set.seed(iSplit)
  
  input = as.matrix(my.datf.train[,2:ncol(my.datf.train)])
  rferf.fit = rfe(input,net.c,sizes = sizes,rfeControl=rfeCtr,prob.model =F,trControl=ctrl ,allowParallel=T, tuneGrid=expand.grid(.mtry=c(floor(sqrt(ncol(input))),ncol(input))), metric='RMSE')
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- rferf.fit$results[rferf.fit$results$Variables ==rferf.fit$bestSubset,2]
  R2.tr    <- rferf.fit$results[rferf.fit$results$Variables ==rferf.fit$bestSubset,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- rferf.fit$results[rferf.fit$results$Variables ==rferf.fit$bestSubset,4]
    R2sd.tr   <- rferf.fit$results[rferf.fit$results$Variables ==rferf.fit$bestSubset,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!  TODOOOOOOOOOOOOOO
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  rfesvm.train.res <- rferf.fit$results[ rferf.fit$results$Variables== rferf.fit$bestSubset, c(2,3)]
  rfesvm.test.res  <- postResample(predict(rferf.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(rferf.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(rferf.fit,my.datf.test)  # predicted Y for test
  noFeats.fit <- length(predictors(rferf.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(rferf.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(rferf.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # ----------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   "RMSE.tr"   = as.numeric(RMSE.tr),
                   "R2.tr"     = as.numeric(R2.tr),
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr),
                   "R2sd.tr"   = as.numeric(R2sd.tr),
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((rfesvm.test.res["RMSE"])),
                   "R2.ts"   = as.numeric((rfesvm.test.res["Rsquared"])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile, append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),       file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",         file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),        file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rferf.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.train.res,file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.test.res, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- importance(rferf.fit$fit, scale = T)
    FeatImp = FeatImp[order(FeatImp,decreasing=T),]
    
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- rferf.fit$fit
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- abs(pred.both-ds.full[,1])    # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)

    if(det.Traind.pls!=0){

    Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
    Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
    
    # Leverage / Hat values
    hat.fit <- Hat.test          # hat values
    hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
    hat.mean <- mean(hat.fit)               # mean hat values
    hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
    
    write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
    
    #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
    thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
    hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
    
    write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    }
    
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    
    dotchart(FeatImp,main="Feature Importance") # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(pred.both,resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    if(det.Traind.pls!=0){
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
  }

  return(list(stat.values=my.stats, model=rferf.fit))  # return a list with statistics and the full model
}

findResamps.funct<- function(caret.obj){
  #=============================================================================
  # A function to find the number of re-samples for caret, rfe or sbf objects 
  # from caret package 
  #=============================================================================
  #caret.obj== caret object of class train, rfe or sbf 

	in.caret.obj<- caret.obj$control$index
	return(length(in.caret.obj))
}

impute.funct<- function(ds,FUN=mean){
  #=============================================================================
  # A function to impute missing values from columns of matrix or data frame
  # using the mean value as the default 
  #=============================================================================
  #ds== data.frame or matrix to be imputed 

	sum.na <- apply(ds,2,sum)
	ind.na <- which(is.na(sum.na)!=FALSE)

	ds.imputeV <- apply(as.matrix(ds[,ind.na]),2,function(x)FUN(x,na.rm=T))
	ds.imputeI <- apply(as.matrix(ds[,ind.na]),2,function(x)which(is.na(x)))

	if(is.list(ds.imputeI)!=TRUE){ds.imputI<- list(ds.imputeI)}

	for(i in 1:length(ds.imputI)){ds[ds.imputI[[i]],ind.na[i]]<- ds.imputeV[i]}
	return(ds)
}

ENETreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #================================
  # 8.4 ElasticNet Regression (caret/glmnet)
  #================================
  
  #library(caret)
  #library(glmnet)
  data(glmnetModel)
  
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "glmnet" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method = sCV, number = 10,repeats = 10,verboseIter=F,#number=10,repeats=10,
                      summaryFunction = defaultSummary)
  tuneGrid=expand.grid(.alpha = seq(0.1,1,length=10),.lambda=99 )
  
  # Train the model using only training set
  set.seed(iSplit)
  enet.fit<- train(net.c~.,data=my.datf.train,
                   method=glmnetModel, tuneLength = 20, trControl = ctrl,family="gaussian",
                   metric='RMSE',tuneGrid=tuneGrid)
  
  pos = which.min(abs(enet.fit$finalModel$lambda-enet.fit$finalModel$lambdaOpt))
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- enet.fit$results[,3]
  R2.tr    <- enet.fit$results[,4]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- enet.fit$results[,5]
    R2sd.tr   <- enet.fit$results[,6]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- getTrainPerf(enet.fit)
  lm.test.res  <- postResample(predict(enet.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(enet.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(enet.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(enet.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(enet.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(enet.fit,ds.full)       # predicted Y
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # ------------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
                   
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                   "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(enet.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(lm.train.res,outFile)
    #write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,quote=F)
    AppendList2CSv(lm.test.res,outFile)
    #write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(enet.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- enet.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(enet.fit) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)

    if(det.Traind.pls!=0){
     Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
     Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
     
     # Leverage / Hat values
     hat.fit <- Hat.test          # hat values
     hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
     hat.mean <- mean(hat.fit)               # mean hat values
     hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
     
     write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
     
     #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
     thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
     hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
     
     write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
     write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
     
     # Cook's distance ?
    }
    # Influence ?
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(enet.fit),residuals(enet.fit),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    if(det.Traind.pls!=0){
     plot(hat.fit, type = "h",
          main="Leverage for Fitted Model",
          xlab="Index", ylab="Hat")
     abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }
    dev.off()
    # --------------------------------------------------------------
  }

  return(list(stat.values=my.stats, model=enet.fit))  # return a list with statistics and the full model
}


###############################################################################################
# RRegrs MAIN FUNCTION 
###############################################################################################

RRegrs <- function(DataFileName="ds.House.csv",DataFileSep=",",DependentVariable=NULL,IndependentVariables=NULL,
  PathDataSet="DataResults",noCores=1,
  ResAvgs="RRegsResAvgs.csv",ResBySplits="RRegrsResAllSplits.csv",ResBest="RRegrsResBest.csv",
  fDet="T",fFilters="F",fScaling="T",fRemNear0Var="T",fRemCorr="T",
  fLM="T",fGLM="T",fPLS="T",fLASSO="T",fSVRM="T",fNN="T",fRF="T",fRFRFE="T",fSVMRFE="T",fENET="T",
  RFE_SVM_C="1;5;10;15;20",RFE_SVM_epsilon="0.01;0.1;0.3",
  cutoff=0.9,iScaling=1,iScalCol=1,trainFrac=0.75,iSplitTimes=10,noYrand=100,
  CVtypes="repeatedcv;LOOCV",NoNAValFile="ds.NoNA.csv",
  No0NearVarFile="ds.No0Var.csv",ScaledFile="ds.scaled.csv",NoCorrFile="ds.scaled.NoCorrs.csv",
  lmFile="LM.details.csv",glmFile="GLM.details.csv",plsFile="PLS.details.csv",
  lassoFile="Lasso.details.csv",svrmFile="SVMRadial.details.csv",
  nnFile="NN.details.csv",rfFile="RF.details.csv",rfrfeFile="RFRFE.details.csv",svmrfeFile="SVMRFE.details.csv",
  enetFile="ENET.details.csv",fR2rule="T") { # input = file with all parameters

  methodCount = 0;
  if (fLM=="T") methodCount = methodCount + 1;
  if (fGLM=="T") methodCount = methodCount + 1;
  if (fPLS=="T") methodCount = methodCount + 1;
  if (fLASSO=="T") methodCount = methodCount + 1;
  if (fSVRM=="T") methodCount = methodCount + 1;
  if (fNN=="T") methodCount = methodCount + 1;
  if (fRF=="T") methodCount = methodCount + 1;
  if (fRFRFE=="T") methodCount = methodCount + 1;
  if (fSVMRFE=="T") methodCount = methodCount + 1;
  if (fENET=="T") methodCount = methodCount + 1;

  if (methodCount < 2) stop("You must select at least two modelling methods to compare.");
  
  # fRBFdda="T", rbfDDAFile="RBF_DDA.details.csv",negThrStep=0.5
  
  # Minimal use:
  # RRegrs()                              # all default params
  # RRegrs(DataFileName="MyDataSet.csv")
  # RRegrs(DataFileName="MyDataSet.csv",PathDataSet="MyResultsFolder")
  # Default: all methods, no feature selection
  ptmTot <- proc.time() # total time

  # ----------------------------------
  # Parallel support
  # ----------------------------------
  if (noCores==0 | noCores>1){ # all available CPU cores or specific no of cores (if noCores = 1, no parallel support!)
    #noCoresSys=as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')) # automatically detected no. of CPU cores
    #library(parallel)
    noCoresSys=detectCores()
    
    if (noCores==0 | noCores>noCoresSys){ # all available CPU cores or the specific cores is greater than the available ones
      noCores=noCoresSys # use the available no of cores
    }

    # parallel for Linux or Mac:
    # ------------------------------------------
    if ( Sys.info() [['sysname']] == "Linux" | Sys.info() [['sysname']] == "Darwin" ){
      #library(doMC)
    }
    # ------------------------------------------
    # parallel for windows:
    # ------------------------------------------
    if (Sys.info()[['sysname']]=="Windows"){
      #library(doSNOW)
      #library(foreach)
    }  
  }     
  
  #==========================================================================================
  # (1) Load dataset and parameters
  #==========================================================================================
  
  # (1.1) PARAMETERS
  #------------------------------------------
  # Write parameter file
  #------------------------------------------
  dir.create(PathDataSet, showWarnings = FALSE)
  ParamFile <- file.path(PathDataSet, "Parameters.csv") # file to output the parameters

  # define a data frame with all parameters of the current calculation
  Params.df = data.frame(RRegrs.Parameters="DataFileName",Parameter.Value=as.character(DataFileName),Description="Input dataset file (Step 1)") # data frame with used parameters
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="PathDataSet",Parameter.Value=as.character(PathDataSet),Description="Working folder for all input and output files"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="noCores",Parameter.Value=as.character(noCores),Description="No of CPU cores (0=all available; 1=no parallel; >1 = specific no. of cores)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="ResAvgs",Parameter.Value=as.character(ResAvgs),Description="Output file averaged statistics (by splits) for each regression method"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="ResBySplits",Parameter.Value=as.character(ResBySplits),Description="Output file statistics for each splitting and each regression method"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="ResBest",Parameter.Value=as.character(ResBest),Description="Output file statistics for the best model"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fDet",Parameter.Value=as.character(fDet),Description="If calculate and print details for all the functions"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fFilters",Parameter.Value=as.character(fFilters),Description="If run Filters (Step 2)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fScaling",Parameter.Value=as.character(fScaling),Description="If Scaling dataset (Step 3)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRemNear0Var",Parameter.Value=as.character(fRemNear0Var),Description="If run Removal of near zero variance columns (Step 4)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRemCorr",Parameter.Value=as.character(fRemCorr),Description="If run Removal of correlated columns (Step 5)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fLM",Parameter.Value=as.character(fLM),Description="If run LM (Step 8.1)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fGLM",Parameter.Value=as.character(fGLM),Description="If run GLM (Step 8.2)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fPLS",Parameter.Value=as.character(fPLS),Description="If run PLS (Step 8.3)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fLASSO",Parameter.Value=as.character(fLASSO),Description="If run LASSO (Step 8.4)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fENET",Parameter.Value=as.character(fENET),Description="If run ENET (Step 8.5)"))
  # Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRBFdda",Parameter.Value=as.character(fRBFdda),Description="If run RBF DDA (Step 8.6)"))
  # Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="negThrStep",Parameter.Value=as.character(negThrStep),Description="Negative Threshold step parameter for RBF DDA (Step 8.6)"))

  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fSVRM",Parameter.Value=as.character(fSVRM),Description="If run svmRadial.RMSE (Step 8.7)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fNN",Parameter.Value=as.character(fNN),Description="If run Neural Networks (Step 8.8)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRF",Parameter.Value=as.character(fRF),Description="If run Random Forest (Step 8.9)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRFRFE",Parameter.Value=as.character(fRFRFE),Description="If run Random Forest RFE (Step 8.10)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fSVMRFE",Parameter.Value=as.character(fSVMRFE),Description="If run Random Forest (Step 8.11)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="RFE_SVM_C",Parameter.Value=as.character(RFE_SVM_C),Description="Values of C for SVM RFE"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="RFE_SVM_epsilon",Parameter.Value=as.character(RFE_SVM_epsilon),Description="Values of epsilon for SVM RFE"))
  
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="cutoff",Parameter.Value=as.character(cutoff),Description="Cut-off for correlated features (default = 0.9)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="iScaling",Parameter.Value=as.character(iScaling),Description="Type of scaling: 1 = normalization; 2 = standardization; 3 = other; any other: no scaling"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="iScalCol",Parameter.Value=as.character(iScalCol),Description="Scaling columns: 1 = including dependent variable; 2: only all the features"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="trainFrac",Parameter.Value=as.character(trainFrac),Description="Fraction of training set from the entire dataset; the rest of dataset is the test set"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="iSplitTimes",Parameter.Value=as.character(iSplitTimes),Description="Number of splitting the dataset into train and test (default  = 10)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="noYrand",Parameter.Value=as.character(noYrand),Description="Number of Y-Randomization (default = 100)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="CVtypes",Parameter.Value=as.character(CVtypes),Description="Cross-validation types: 10-CV (repeatedcv) and LOOCV"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="NoNAValFile",Parameter.Value=as.character(NoNAValFile),Description="Dataset without NA values (if fDet is True)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="No0NearVarFile",Parameter.Value=as.character(No0NearVarFile),Description="Dataset without zero near features from Step 3 (if fDet is True)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="ScaledFile",Parameter.Value=as.character(ScaledFile),Description="Scaled dataset file from Step 4 (if fDet is True)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="NoCorrFile",Parameter.Value=as.character(NoCorrFile),Description="Dataset after correction removal in Step 5 (if fDet is True)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="lmFile",Parameter.Value=as.character(lmFile),Description="LM output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="glmFile",Parameter.Value=as.character(glmFile),Description="GLM output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="plsFile",Parameter.Value=as.character(plsFile),Description="PLS output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="lassoFile",Parameter.Value=as.character(lassoFile),Description="Lasso output file with details"))
  #Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="rbfDDAFile",Parameter.Value=as.character(rbfDDAFile),Description="RBF DDA output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="svrmFile",Parameter.Value=as.character(svrmFile),Description="SVM Radial output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="nnFile",Parameter.Value=as.character(nnFile),Description="NN output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="rfFile",Parameter.Value=as.character(rfFile),Description="RF output"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="rfrfeFile",Parameter.Value=as.character(rfrfeFile),Description="RF-RFE output"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="svmrfeFile",Parameter.Value=as.character(svmrfeFile),Description="SVM-RFE output"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="enetFile",Parameter.Value=as.character(enetFile),Description="ENET output"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fR2rule",Parameter.Value=as.character(fR2rule),Description="Best model rule: R2 (default = T) or adjR2 (F)"))

  write.csv(Params.df,file=ParamFile,row.names=F,quote=F) # write parameters to a CSV in the working folder

  # Get calculation parameters 
  fDet         = as.logical(fDet)         # flag to calculate and print details for all the functions
  fFilters     = as.logical(fFilters)     # flag to apply filters                          (2)
  fScaling     = as.logical(fScaling)     # flag for dataset Scaling                       (3)
  fRemNear0Var = as.logical(fRemNear0Var) # flag for Removal of near zero variance columns (4)
  fRemCorr     = as.logical(fRemCorr)     # flag for Removal of correlated columns         (5)
  
  cutoff       = as.numeric(as.character(cutoff))  # cut off for correlated features
  fLM          = as.logical(fLM)     # flag to run LM            (8.1)
  fGLM         = as.logical(fGLM)    # flag to run GLM           (8.2)
  fPLS         = as.logical(fPLS)    # flag to run PLS           (8.3)
  fLASSO       = as.logical(fLASSO)  # flag to run LASSO         (8.4)
  fenet        = as.logical(fENET)   # flag to run ElasticNet    (8.5)
  #fRBFdda      = as.logical(fRBFdda) # flat to run RBF DDA      (8.6)
  fSVRM        = as.logical(fSVRM)   # flat to run svmRadial     (8.7)
  fNN          = as.logical(fNN)     # flat to run NN            (8.8)
  fRF          = as.logical(fRF)     # flag to run RandomForest  (8.9)
  fRFRFE       = as.logical(fRFRFE)  # flag to run RF-RFE        (8.10)
  fSVMRFE      = as.logical(fSVMRFE) # flag to run SVM RFE       (8.11)
  
  rfe_SVM_param_c   = strsplit(as.character(RFE_SVM_C),";")[[1]] # values of C for SVM RFE
  rfe_SVM_param_eps = strsplit(as.character(RFE_SVM_epsilon),";")[[1]] # values of epsilon for SVM RFE
  # negThrStep

  fR2rule = as.logical(fR2rule) # flag to decide order rule for the best model (True for R2 and False for adjR2, default = True)
  # ----------------------------------------------------------------------------------------
  trainFrac   = as.numeric(as.character(trainFrac))   # the fraction of training set from the entire dataset; trainFrac = the rest of dataset, the test set
 
  CVtypes = strsplit(as.character(CVtypes),";")[[1]] # types of cross-validation methods
  CVtypes2 = c("repeatedcv") # for complex methods we run only 10-fold CV even the user is using other parameters!
  
  # Generate path + file name = original dataset
  if (file.exists(DataFileName)) { # is it a full path already?
    inFile <- DataFileName
  } else {
    inFile <- file.path(PathDataSet, DataFileName)
  }
  
  sDescription=paste("=======================================================================================================",
    "RRegrs - R Regression Models",
    "Get the best regression models for one dataset using R caret methods", "eNanoMapper.net","AUTHORS:",
    "Georgia Tsiliki: ChemEng - NTUA, Greece, g_tsiliki@hotmail.com",
    "Cristian R. Munteanu: RNASA-IMEDIR, University of A Coruna, Spain, muntisa@gmail.com",
    "Jose A. Seoane: Stanford Cancer Institute, USA, seoane@stanford.edu",
    "Carlos Fernandez-Lozano: RNASA-IMEDIR, University of A Coruna, Spain, carlos.fernandez@udc.es",
    "Haralambos Sarimveis: ChemEng - NTUA, Greece, hsarimv@central.ntua.gr",
    "Egon Willighagen: BiGCaT - Maastricht University, The Netherlands, egon.willighagen@gmail.com",
    "=======================================================================================================",sep="\n")
  cat(sDescription) # print package header information

  # -----------------------------------
  # (1.2) Load the ORIGINAL DATASET
  # -----------------------------------
  cat("\n-> Loading original dataset ...\n")  # it can contain errors, correlations, near zero variance columns
  cat("      ---> ",inFile,"\n")
  ds.dat0 <- read.csv(inFile,header=T,sep=DataFileSep)          # original dataset frame
  
  # resolving the text to number errors for future calculations
  if (missing(IndependentVariables)) {
    ds.indx<- colnames(ds.dat0)[2:dim(ds.dat0)[2]]                    # FEATURE names (no dependent variable)
    ds.dat1<- ds.dat0[1:dim(ds.dat0)[1],2:dim(ds.dat0)[2]]            # dataset as columns
  } else {
    cat("      ---> independent variables : ", IndependentVariables,"\n")
    ds.indx<- IndependentVariables                    # FEATURE names (no dependent variable)
    ds.dat1<- ds.dat0[1:dim(ds.dat0)[1],IndependentVariables]            # dataset as columns
  }
  ds.dat1<- apply(ds.dat1,1,function(x)as.numeric(as.character(x))) # dataset as row vectors to be used with caret!!!
  
  # dependent variable
  if (missing(DependentVariable)) {
    net.c<- ds.dat0[,1]
  } else {
    cat("      ---> dependent variable : ", DependentVariable,"\n")
    net.c<- ds.dat0[,DependentVariable]
  }
  net.c<- as.numeric(as.character(net.c)) # values
  # full ds frame with training and test
  ds<- as.data.frame(cbind(net.c,t(ds.dat1)))
  cat(colnames(ds),"\n")
  cat(as.character(head(ds)))
  
  #========================================================
  # (2) FILTERS
  #     (it will be implemented in the future versions)
  #========================================================
  # 2.1 Outlier removal
  # 2.2 Custom filter (percentage threshold)
  # 2.3 Processing of missing values - use of preProcess();
  #     caret employs knnImpute algorithm to impute values from a neighborhood of k
  # TO BE IMPLEMENTED

  # -----------------------------------------------------------------------
  # (2) Remove NA values
  # -----------------------------------------------------------------------
  if(length(which(is.na(ds)==TRUE))!=0){
    cat("-> Removal of NA values ...\n")
    outFile <- file.path(PathDataSet,NoNAValFile) # the same folder as input  
    
    # get the ds without NA values-- currently use the default function (mean) 
    ds <- impute.funct(ds)
    if (fDet == TRUE){   # write as details the corrected ds file
     write.csv(ds, outFile,row.names=F, quote=F)}
  }
 
  # -----------------------------------------------------------------------
  # (3) Remove near zero variance columns
  # -----------------------------------------------------------------------
  if (fRemNear0Var==T) {
    cat("-> Removal of near zero variance columns ...\n")
    outFile <- file.path(PathDataSet,No0NearVarFile) # the same folder as input  
    
    # get the ds without near zero cols 
    ds <- cbind("net.c" = ds[,1],RemNear0VarCols(ds[,2:dim(ds)[2]],fDet,outFile))
    # use df without Y (predicted values), reconstruct the ds
    # inputs: ds, flag for details, output file
  }
  
  # -----------------------------------------------------------------------
  # (4) Scaling dataset: normalization (default), standardization, other
  # -----------------------------------------------------------------------
  if (fScaling==T) {
    cat("-> Scaling original dataset ...\n")
    outFile <- file.path(PathDataSet,ScaledFile)       # the same folder as input
    
    # run fuction for scaling input dataset file
    ds <- ScalingDS(ds,iScaling,iScalCol,fDet,outFile)
    # use df without Y (predicted values), reconstruct the ds
    # inputs: ds, type of scaling, flag for details, starting column, output file
  }
  
  # -----------------------------------------------------------------------
  # (5) Remove correlated features
  # -----------------------------------------------------------------------
  if (fRemCorr==T) {    
    cat("-> Removing correlated features ...\n") 
    outFile <- file.path(PathDataSet,NoCorrFile)    # the same folder as the input
    
    # run function to remove the correlations between the features
    ds <- cbind("net.c" = ds[,1],RemCorrs(ds[,2:dim(ds)[2]],fDet,cutoff,outFile))
  }
  
  # Check data has at least 5 columns for meaningful analysis
  if(dim(ds)[2] < 5 || dim(ds)[1] < 3){
  print(c(dim(ds)))
  stop(paste("Your corrected data set has dimensions:", paste(as.character(dim(ds)),collapse=', '),". Try repeating analysis without filtering options.",sep=''))
  }

  # print no of CPU cores used for calculation
  # noCoresSys=as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')) # automatically detected no. of CPU cores
  #library(parallel)
  noCoresSys=detectCores()
  if (noCores==0){ cat("       -> CPU Cores = ",noCoresSys,"(only complex methods)\n") }
  else{ cat("                  -> CPU Cores = ",noCores,   "(only complex methods)\n") }
  
  #=========================================================================================================
  # Steps 6 - 8 will be repeated 10 times for reporting each result and average (iSplitTimes = 10, default)
  #=========================================================================================================
  
  # Initialize the list with the statistics results; the same HEADER as the function output
  dfRes <- list("RegrMeth"     = NULL,
                "Split No"     = NULL,    
                "CVtype"       = NULL,      
                "NoModelFeats" = NULL,
                "ModelFeats"   = NULL,
                "adjR2.tr"  = NULL,
                "RMSE.tr"   = NULL,
                "R2.tr"     = NULL,
                "RMSEsd.tr" = NULL,
                "R2sd.tr"   = NULL,
                "adjR2.ts"= NULL,
                "RMSE.ts" = NULL,
                "R2.ts"   = NULL,
                "corP.ts" = NULL,
                "adjR2.both" = NULL,
                "RMSE.both"  = NULL,
                "R2.both"    = NULL)
  #-------------------------------------------------------------------------------------------------
  for (i in 1:iSplitTimes) {   # Step splitting number = i
    # -----------------------------------------------------------------------
    # (6) Dataset split: Training and Test sets
    # -----------------------------------------------------------------------
    cat("-> Splitting dataset in Training and Test sets ...\n")
    cat(paste("--> Split No.",i,"from",iSplitTimes,"\n"))

    # Initialize the list with the statistical models for all types of CV; per iSplitTimes
    dfMod <- sapply(CVtypes,function(x) NULL)
    for(cv in 1:length(CVtypes)){class(dfMod[[cv]])<- 'list'; names(dfMod)[[cv]]<- CVtypes[cv]}
    mod.ind<- rep(1,length(CVtypes)) # dummy variable to indicate the index of each new dfMod entry (per CVtype) 
    
    iSeed=i                 # to reapeat the ds splitting, different values of seed will be used
    dsList  <- DsSplit(ds,trainFrac,fDet,PathDataSet,iSeed) # return a list with 2 datasets = dsList$train, dsList$test
    
    # get train and test from the resulted list
    ds.train<- dsList$train
    ds.test <- dsList$test
    
    # -----------------------------------------------------------------------
    # (7) Feature selection
    # -----------------------------------------------------------------------
    # TO BE IMPLEMENTED
    
    # -----------------------------------------------------------------------
    # (8) REGRESSION METHODS
    # -----------------------------------------------------------------------
    
    # ----------------------------------
    # Parallel support - registed cores
    # ----------------------------------
    if (noCores!=1){
      # ------------------------------------------
      # parallel for Linux or Mac:
      # ------------------------------------------
      if (Sys.info()[['sysname']]=="Linux" | Sys.info()[['sysname']]=="Darwin"){
        registerDoMC(cores = noCores) # CPU cores
      }
      # ------------------------------------------
      # parallel for windows:
      # ------------------------------------------
      if (Sys.info()[['sysname']]=="Windows"){
       cl<-makeCluster(noCores,outfile="") 
        registerDoSNOW(cl)
      }   
    }
    
    # --------------------------------------------
    # 8.1. Basic LM : default
    # --------------------------------------------
    if (fLM==T) {   # if LM was selected, run the method
      outFile.LM <- file.path(PathDataSet,lmFile)   # the same folder as the input is used for the output
      cat("-> LM : Linear Multi-regression ...\n")
      
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) { # there is no CV but it will be implemented in the future!!!
        cat("    -->",CVtypes[cv],"\n")
        ptmLM <- proc.time()
        lm.model <- LMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.LM) # run GLM for each CV and regr method
        print(proc.time() - ptmLM) # print running time
        
        my.stats.LM <- lm.model$stat.values # stat values
        my.model.LM <- lm.model$model # model 
        
        #-------------------------------------------------------
        # Add output from GLM to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.LM, dfRes, SIMPLIFY=F)
        
        # List of models for each splitting, CV type & regression method
        names1 <- strsplit(deparse(quote(my.model.LM)),'my.model.')[[1]][2]
        dfMod[[cv]]$names1 <- my.model.LM  
        names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
        mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
      } # end CV types
    } # end LM
    
    # -----------------------------------------------------------------------------------------
    # (8.2) GLM based on AIC regression - Generalized Linear Model with Stepwise Feature Selection
    # -----------------------------------------------------------------------------------------
    if (fGLM==T) {   # if GLM was selected, run the method
      outFile.GLM <- file.path(PathDataSet,glmFile)   # the same folder as the input is used for the output
      
      cat("-> GLM : Generalized Linear Model stepwise - based on AIC ...\n")
		  # For each type of CV do all the statistics
		  # -----------------------------------------------------
		  for (cv in 1:length(CVtypes)) {
		    cat("    -->",CVtypes[cv],"\n")
		    ptmGLM <- proc.time()
		    glm.model  <- GLMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.GLM) # run GLM for each CV and regr method
		    print(proc.time() - ptmGLM) # print running time
		    
		    my.stats.GLM <- glm.model$stat.values # stat values
		    my.model.GLM <- glm.model$model # model
		    #my.stats.split <- c(my.stats.dsInfo,my.stats.GLM) # merge the ds info with statistics results for each Cv & reg method
		    
		    #-------------------------------------------------------
		    # Add output from GLM to the list of results
		    #-------------------------------------------------------
		    # List of results for each splitting, CV type & regression method
		    dfRes = mapply(c, my.stats.GLM, dfRes, SIMPLIFY=F)
		    
		    # List of models for each splitting, CV type & regression method
		    names1 <- strsplit(deparse(quote(my.model.GLM)),'my.model.')[[1]][2]
		    dfMod[[cv]]$names1 <- my.model.GLM  
		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		  } # end CV types  
    } # end GLM
    
    # --------------------------------------------
    # 8.3. PLS
    # --------------------------------------------
    if (fPLS==T) {   # if PLS was selected, run the method
		outFile.PLS <- file.path(PathDataSet,plsFile)   # the same folder as the input is used for the output
      
		cat("-> PLS : Partial Least Squares Regression ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes)) {
			cat("    -->",CVtypes[cv],"\n")
			ptmPLS <- proc.time()
			pls.model <- PLSreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.PLS) # run PLS for each CV and regr method
			print(proc.time() - ptmPLS) # print running time
        
			my.stats.PLS <- pls.model$stat.values
			my.model.PLS <- pls.model$model
			#-------------------------------------------------------
			# Add output from PLS to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.PLS, dfRes, SIMPLIFY=F)
        
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.PLS)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.PLS  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
    } # end PLS
    
    #cat("-> PLS Wrapper Feature Selection ...\n")
    ## For each type of CV do all the statistics
    ## -----------------------------------------------------
    #for (cv in 1:length(CVtypes)) {
    #  cat("    -->",CVtypes[cv],"\n")
    #  ptmPLSw <- proc.time()
    #  pls.model <- PLSregWSel(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.PLS) # run PLSw for each CV and regr method
    #  print(proc.time() - ptmPLSw) # print running time
    
    #  my.stats.PLS <- pls.model$stat.values
    #  my.model.PLS <- pls.model$model
    
    #  #-------------------------------------------------------
    #  # Add output from PLSw to the list of results
    #  #-------------------------------------------------------
    #  # List of results for each splitting, CV type & regression method
    #  dfRes = mapply(c, my.stats.PLS, dfRes, SIMPLIFY=F)
    #} # end CV types 
    
    # --------------------------------------------
    # 8.4. LASSO regression
    # --------------------------------------------
    if (fLASSO==T) {   # if LASSO was selected, run the method
		outFile.LASSO <- file.path(PathDataSet,lassoFile)   # the same folder as the input is used for the output
      
		cat("-> Lasso ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes2)) {
			cat("    -->",CVtypes2[cv],"\n")
			ptmLASSO <- proc.time()
			lasso.model <- LASSOreg(ds.train,ds.test,CVtypes2[cv],i,fDet,outFile.LASSO) # run LASSO for each CV and regr method
			print(proc.time() - ptmLASSO) # print running time

			my.stats.LASSO <- lasso.model$stat.values
			my.model.LASSO <- lasso.model$model
			#-------------------------------------------------------
			# Add output from Lasso to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.LASSO, dfRes, SIMPLIFY=F)
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.LASSO)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.LASSO  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
    } # end Lasso
    
    # --------------------------------------------
    # 8.5. Elastic Net regression
    # --------------------------------------------
    if (fenet==T) {   # if ENET was selected, run the method
		outFile.ENET <- file.path(PathDataSet,enetFile)   # the same folder as the input is used for the output
      
		cat("-> ENET : Elastic Nets ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes)) {
			cat("    -->",CVtypes[cv],"\n")
			ptmENET <- proc.time()
			enet.model <- ENETreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.ENET) # run elastic net for each CV and regr method
			print(proc.time() - ptmENET) # print running time
        
			my.stats.ENET <- enet.model$stat.values
			my.model.ENET <- enet.model$model 
        
			#-------------------------------------------------------
			# Add output from ENET to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.ENET, dfRes, SIMPLIFY=F)
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.ENET)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.ENET  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types  
    } # end enet

#     # ----------------------------------------------------------------
#     # 8.6. RBF network with the DDA algorithm regression (caret)
#     # ----------------------------------------------------------------
#     if (fRBFdda==T) {   # if RBF-DDA was selected, run the method
# 		outFile.rbfDDA <- file.path(PathDataSet,rbfDDAFile)   # the same folder as the input is used for the output
# 		  
# 		cat("-> RBF-DDA : Radial Basis Functions - Dynamic Decay Adjustment ...\n")
# 		# For each type of CV do all the statistics
# 		# -----------------------------------------------------
# 		for (cv in 1:length(CVtypes2)) {
# 		    cat("    -->",CVtypes2[cv],"\n")
# 		    ptmRBF_DDA <- proc.time()
# 		    rbfDDA.model <- RBF_DDAreg(ds.train,ds.test,CVtypes2[cv],negThrStep,i,fDet,outFile.rbfDDA) # run rbfDDA for each CV and regr method
# 		    print(proc.time() - ptmRBF_DDA) # print running time
# 		    
# 		    my.stats.rbfDDA <- rbfDDA.model$stat.values
# 		    my.model.rbfDDA <- rbfDDA.model$model 
# 		    #-------------------------------------------------------
# 		    # Add output from RBF-DDA to the list of results
# 		    #-------------------------------------------------------
# 		    # List of results for each splitting, CV type & regression method
# 		    dfRes = mapply(c, my.stats.rbfDDA, dfRes, SIMPLIFY=F)
# 		    
# 		    # List of models for each splitting, CV type & regression method
# 		    names1 <- strsplit(deparse(quote(my.model.rbfDDA)),'my.model.')[[1]][2]
# 		    dfMod[[cv]]$names1 <- my.model.rbfDDA  
# 		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
# 		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
# 		} # end CV types
#     } # end rbfDDA
    
    # --------------------------------------------
    # 8.7. SVM radial regression
    # --------------------------------------------
    if (fSVRM==T) {   # if SVM Radial was selected, run the method
		outFile.SVRM <- file.path(PathDataSet,svrmFile)   # the same folder as the input is used for the output
      
		cat("-> SVM radial : Support vector machine using radial functions ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes)) {
		    cat("    -->",CVtypes[cv],"\n")
		    ptmSVRM <- proc.time()
		    SVRM.model <- SVRMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.SVRM,rfe_SVM_param_c) # run SVRM Radial for each CV and regr method
		    print(proc.time() - ptmSVRM) # print running time
		    
		    my.stats.SVRM <- SVRM.model$stat.values
		    my.model.SVRM <- SVRM.model$model 
		    
		    #-------------------------------------------------------
		    # Add output from SVM Radial to the list of results
		    #-------------------------------------------------------
		    # List of results for each splitting, CV type & regression method
		    dfRes = mapply(c, my.stats.SVRM, dfRes, SIMPLIFY=F)
		    # List of models for each splitting, CV type & regression method
		    names1 <- strsplit(deparse(quote(my.model.SVRM)),'my.model.')[[1]][2]
		    dfMod[[cv]]$names1 <- my.model.SVRM  
		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
    } # end SVRM
    
    # --------------------------------------------
    # 8.8. Neural Networks Regression
    # --------------------------------------------
    if (fNN==T) {   # if NNet was selected, run the method
      outFile.NN <- file.path(PathDataSet,nnFile)   # the same folder as the input is used for the output
      
		 cat("-> NN : Neural Networks ...\n")
		 # For each type of CV do all the statistics
		 # -----------------------------------------------------
		 for (cv in 1:length(CVtypes)) {
			cat("    -->",CVtypes[cv],"\n")
		    ptmNN <- proc.time()
		    nn.model <- NNreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.NN) # run NNet for each CV and regr method
		    print(proc.time() - ptmNN) # print running time
		    
		    my.stats.NN <- nn.model$stat.values
		    my.model.NN <- nn.model$model
		    
		    #-------------------------------------------------------
		    # Add output from NNet to the list of results
		    #-------------------------------------------------------
		    # List of results for each splitting, CV type & regression method
		    dfRes = mapply(c, my.stats.NN, dfRes, SIMPLIFY=F)
		    # List of models for each splitting, CV type & regression method
		    names1 <- strsplit(deparse(quote(my.model.NN)),'my.model.')[[1]][2]
		    dfMod[[cv]]$names1 <- my.model.NN  
		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types  
    } # end NNet
    
    # --------------------------------------------
    # 8.9. Random Forest Regression (RF)
    # --------------------------------------------
    if (fRF==T) {   # if RF was selected, run the method
		outFile.RF <- file.path(PathDataSet,rfFile)   # the same folder as the input is used for the output
		cat("-> RF : Random Forest ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes2)) {
			cat("    -->",CVtypes2[cv],"\n")
			ptmRF <- proc.time()
			rf.model <- RFreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.RF) # run RF for each CV and regr method
			print(proc.time() - ptmRF) # print running time
		   
			my.stats.RF <- rf.model$stat.values
			my.model.RF <- rf.model$model
		  
			#-------------------------------------------------------
			# Add output from RF to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.RF, dfRes, SIMPLIFY=F)
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.RF)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.RF  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
	}
      
	# --------------------------------------------
	# 8.10. RF-RFE: Random Forest Regression Recursive Feature Elimination
	# --------------------------------------------
	if (fRFRFE==T) {   # if RF-RFE was selected, run the method
		outFile.RFRFE <- file.path(PathDataSet,rfrfeFile)   # the same folder as the input is used for the output
      
		cat("-> RF-RFE: Random Forest-Recursive Feature Elimination ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes2)) {
			cat("    -->",CVtypes2[cv],"\n")
			ptmRFRFE <- proc.time()
			rfrfe.model <- RFRFEreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.RFRFE) # run RF for each CV and regr method
			print(proc.time() - ptmRFRFE) # print running time
        
			my.stats.RFRFE <- rfrfe.model$stat.values
			my.model.RFRFE <- rfrfe.model$model
			
			#-------------------------------------------------------
			# Add output from RF to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.RFRFE, dfRes, SIMPLIFY=F)
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.RFRFE)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.RFRFE  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
    } # end RF-REF

    # --------------------------------------------
    # 8.11. SVM-RFE
    # --------------------------------------------
    if (fSVMRFE==T) {   # if SVM-RFE was selected, run the method
		outFile.SVMRFE <- file.path(PathDataSet,svmrfeFile)   # the same folder as the input is used for the output
      
		cat("-> SVM-RFE : Support Vector Machines Recursive Feature Elimination ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes2)) {
		    cat("    -->",CVtypes2[cv],"\n")
		    ptmSVMRFE <- proc.time()
		    svmrfe.model <- SVMRFEreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.SVMRFE,rfe_SVM_param_c,rfe_SVM_param_eps) # run SVM RFEet for each CV and regr method
		    print(proc.time() - ptmSVMRFE) # print running time
		    
		    my.stats.SVMRFE <- svmrfe.model$stat.values
		    my.model.SVMRFE <- svmrfe.model$model 
		    
		    #-------------------------------------------------------
		    # Add output from SVM RFE to the list of results
		    #-------------------------------------------------------
		    # List of results for each splitting, CV type & regression method
		    dfRes = mapply(c, my.stats.SVMRFE, dfRes, SIMPLIFY=F)
		    # List of models for each splitting, CV type & regression method
		    names1 <- strsplit(deparse(quote(my.model.SVMRFE)),'my.model.')[[1]][2]
		    dfMod[[cv]]$names1 <- my.model.SVMRFE  
		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		  } # end CV types  
    } # end SVM RFE
    
    
    # Kill parallel server
    # ------------------------
    if (noCores!=1){
      # ------------------------------------------
      # parallel for windows:
      # ------------------------------------------
      if (Sys.info()[['sysname']]=="Windows"){ # clean the memory!
        stopCluster(cl)
      } 
    }
    # ------------------------
    

    # END OF REGRESSION METHODS/FUNCTIONS

    # -----------------------------------------------------------------------
    # (8.final) Produce comparison plots amongst models 
    # -----------------------------------------------------------------------

    for(cv in 1:length(CVtypes)){
      dfMod.n<- dfMod[[cv]]# keep only models with the same number of resamples
      dfMod.ind<- unlist(lapply(dfMod[[cv]],findResamps.funct))
      dfMod.ind.d<- which(duplicated(dfMod.ind)=='TRUE')[1]
      dfMod.in<- which(dfMod.ind!=dfMod.ind[dfMod.ind.d])
      dfMod.flag<- 0 # flag to indicate that dfMod consists of models with different resamples  
      
      if(is.na(dfMod.ind.d)!= TRUE){if(length(dfMod.in)!=0){dfMod.n<- dfMod.n[-dfMod.in]}}
      else{dfMod.flag<- 1}
      
      if(CVtypes[cv]!='LOOCV' && length(dfMod.n)>=2 && dfMod.flag!=1){ 
        cat("-> Comparisons plots for multiple regression methods ...\n")

        resamps <- resamples(dfMod.n)#,modelNames=names(dfMod[[cv]]))  
        # calculate their differences in terms of R2 and RMSE values
        difValues <- diff(resamps)
        #summary(difValues)
        #plot different models in terms of R2 adn RMSE values in the training set 
        pdf(file=paste(PathDataSet,"/ModelsComp.","iSplits.",i,".pdf",sep=""))
        print(bwplot(resamps, layout = c(2, 1),main=paste('Resampling results on the training set',' (data split ',i,')',sep='')))
        dev.off()    

        #plot differences of models in terms of R2 adn RMSE values in the training set 
        pdf(file=paste(PathDataSet,"/DifModels.R2.","iSplits.",i,".pdf",sep=""))
        print(dotplot(difValues,metric='Rsquared',main=paste('Models` differences on the training set',' (data split ',i,')',sep='')))
        dev.off()
        
        pdf(file=paste(PathDataSet,"/DifModels.RMSE.","iSplits.",i,".pdf",sep=""))
        print(dotplot(difValues,metric='RMSE',main=paste('Models` differences on the training set',' (data split ',i,')',sep='')))

        dev.off()
      }
    }
   
  } # END SPLITTING
  
  #------------------------------------------------------------------------------
  # 9. Results for all splittings (not ordered)
  #-------------------------------------------------------------------------------
  cat("-> Results for all splitings ...\n")
  df.res <- data.frame(dfRes)
  # print(df.res) # print all results as data frame
  
  # Writing the statistics into output files: one with detailed splits, other with only averages
  # File names includin paths for the statistics outputs (only averages and split detailed; +averages [to be implemented])
  ResBySplitsF <- file.path(PathDataSet,ResBySplits)   # the main output file with statistics for each split
  write.csv(df.res, file = ResBySplitsF)    # write statistics data frame into a CSV output file
  
  #-------------------------------------------------------------------------------------
  # Averaged values of the results by each Regression Method & CV type
  #-------------------------------------------------------------------------------------
  cat("-> Averaged statistics ...\n")
  ResAvgsF <- file.path(PathDataSet,ResAvgs)           # the main output file with averaged statistics for each regression method
  #library(data.table)
  dt.res  <- data.table(df.res) # convert data frame into data table (for sorting abd averaging)
  
  # MEANS for each Regression Method & CV type
  #--------------------------------------------------------------------------------------------------------------
  # means for all CV types, not only 10CV
  dt.mean <- dt.res[,list(adjR2.tr.Avg=mean(adjR2.tr),RMSE.tr.Avg=mean(RMSE.tr),R2.tr.Avg=mean(R2.tr),
                          RMSEsd.tr.Avg=mean(RMSEsd.tr),R2sd.tr.Avg=mean(R2sd.tr),adjR2.ts.Avg=mean(adjR2.ts),
                          RMSE.ts.Avg=mean(RMSE.ts),R2.ts.Avg=mean(R2.ts),corP.ts.Avg=mean(corP.ts),
                          adjR2.both.Avg=mean(adjR2.both),RMSE.both.Avg=mean(RMSE.both),
                          R2.both.Avg=mean(R2.both),NoModelFeats.Avg=round(mean(NoModelFeats),1)),by="RegrMeth,CVtype"]
  
  dt.mean     <- dt.mean[dt.mean$CVtype=="repeatedcv",]   # keep only the 10CV results to be used to find the best model
  
  if (fR2rule==T) { # R2 rule for the best model (default)
	dt.mean.ord <- dt.mean[order(-rank(R2.ts.Avg))]     # descendent order the averages by R2.ts.Avg
  } else { # adjR2 rule for the best model
	dt.mean.ord <- dt.mean[order(-rank(adjR2.ts.Avg))]  # descendent order the averages by adjR2.ts.Avg
  }
  
  # Write averages descendent ordered by R2.ts.Avg / adjR2.ts.Avg
  #-------------------------------------------------------------------------------
  write.csv(data.frame(dt.mean.ord), file = ResAvgsF)    # write statistics data frame into a CSV output file
  
  #------------------------------------------------------------------------------
  # 10. Best model selection - detailed statistics
  #-------------------------------------------------------------------------------
  cat("-> Best model analysis ...\n")
  
  # Algorithm to verify similar R2 / adjR2 values
  # -> from the best ones (+/- 0.05 of R2 / adjR2), chose the one with less variables, after that the one with min RMSE
  
  best.dt  <- dt.mean.ord[1] # the best model (R2.ts / adjR2.ts) should be the first value in the descendent ordered results

  # best.reg <- paste(best.dt$RegrMeth,collapse="") # best regression method
  
  # Best model rule: R2 or adjR2 for ordering
  # +/- 0.05 R2 / adjR2ts --> min(RMSE)
  # best.adjR2.ts is used for R2 or adjR2 rule
  if (fR2rule==T) { # R2 rule for the best model (default)
	best.adjR2.ts <- as.numeric(data.frame(best.dt)[,10]) # best R2.ts avgs
  } else { # adjR2 rule for the best model
	best.adjR2.ts <- as.numeric(data.frame(best.dt)[,8])  # best adjR2.ts avgs
  }
  
  # best model with R2 or adjR2.ts +/- 0.05 and min of RMSE for Avgs
  if (fR2rule==T) { # R2 rule for the best model (default)
	best.dt  <- dt.mean.ord[R2.ts.Avg %between% c(best.adjR2.ts-0.05,best.adjR2.ts+0.05)][which.min(RMSE.ts.Avg)]
  } else { # adjR2 rule for the best model
	best.dt  <- dt.mean.ord[adjR2.ts.Avg %between% c(best.adjR2.ts-0.05,best.adjR2.ts+0.05)][which.min(RMSE.ts.Avg)]
  }
  
  best.reg <- paste(best.dt$RegrMeth,collapse="") # best regrression method
  cat("    -> Method:",best.reg,"\n")

  # best model non-averaged ? no. of features
  # -----------------------------------------------
  # best.method <- dt.res[CVtype == "repeatedcv"][RegrMeth == best.reg] # best modes corresponding with the avg best values
  # best.method.mean <- mean(as.numeric(data.frame(best.method)[,11])) # best adjR2.ts)
  # dt.res[CVtype == "repeatedcv"][RegrMeth == best.reg][NoModelFeats == min(NoModelFeats)][RMSE.ts == min(RMSE.ts)]
  
  #----------------------------------------------------------
  # 11. Best model detailed statistics 
  #----------------------------------------------------------
  # Write the best model statistics
  ResBestF <- file.path(PathDataSet,ResBest)
  write.table("Averaged values for all spits: ",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
  # write.csv(data.frame(best.dt), file = ResBestF)    # write statistics data frame into a CSV output file
  write.table(data.frame(best.dt), file=ResBestF,append=T,sep=",",col.names=T,quote=F) # write statistics data frame into a CSV output file
  
  # Use the last split for dataset (ds.train & ds.test) ! (or chose other one?)
  
  # Run the caret function with the method from the best method, for one training-test split only
  # and append the details in the best model output file
  
  # ----------------------------------
  # Parallel support - registed cores
  # ----------------------------------
  if (noCores!=1){
    # ------------------------------------------
    # parallel for Linux or Mac:
    # ------------------------------------------
    if (Sys.info()[['sysname']]=="Linux" | Sys.info()[['sysname']]=="Darwin"){
      registerDoMC(cores = noCores) # CPU cores
    }
    # ------------------------------------------
    # parallel for windows:
    # ------------------------------------------
    if (Sys.info()[['sysname']]=="Windows"){
     cl<-makeCluster(noCores,outfile="") 
      registerDoSNOW(cl)
    }   
  }

  if (best.reg=="lm") {
    my.stats.reg  <- LMreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run GLM for each CV and regr method
  }
  if (best.reg=="glmStepAIC") {
    my.stats.reg  <- GLMreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run GLM for each CV and regr method
  }
  if (best.reg=="pls") {
    my.stats.reg  <- PLSreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run PLS for each CV and regr method
  }
  if (best.reg=="lasso.RMSE") {
    my.stats.reg  <- LASSOreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run LASSO for each CV and regr method
  }
  if (best.reg=="svmRadial") {  
    my.stats.reg  <- SVRMreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF,rfe_SVM_param_c)$stat.values # run SVRM Radial for each CV and regr method
  }
  if (best.reg=="nnet") {  
    my.stats.reg  <- NNreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run NNet for each CV and regr method
  } 
  if (best.reg=="rf") {  
    my.stats.reg  <- RFreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run NNet for each CV and regr method
  } 
  if (best.reg=="svmRFE") {  
    my.stats.reg  <- SVMRFEreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF,rfe_SVM_param_c,rfe_SVM_param_eps)$stat.values # run NNet for each CV and regr method
  } 
  if (best.reg=="glmnet") {  
    my.stats.reg  <- ENETreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run ENET for each CV and regr method
  }
  if (best.reg=="rfRFE") {  
    my.stats.reg  <- RFRFEreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run RF RFE for each CV and regr method
  }

#   if (best.reg=="rbfDDA") {  
#     my.stats.reg  <- RBF_DDAreg(ds.train,ds.test,"repeatedcv",negThrStep,i,T,ResBestF)$stat.values # run rbfDDA for each CV and regr method
#  }

  #--------------------------------------------------------------------------------------
  # 12. Test best model with test dataset + Y randomization
  #--------------------------------------------------------------------------------------
  # ratios Yrand R2 - Best model R2 / Best model R2
  R2Diff.Yrand = NA
  if (noYrand >  0) {
    R2Diff.Yrand <- Yrandom(
      ds,trainFrac,best.reg,my.stats.reg$R2.ts,noYrand,ResBestF,rfe_SVM_param_c,rfe_SVM_param_eps
    ) # mean value of ratio (deatails are printed to output file)
  }
  
  # (ex param: negThrStep) 

  # Kill parallel server
  # ------------------------
  if (noCores!=1){
    # ------------------------------------------
    # parallel for windows:
    # ------------------------------------------
    if (Sys.info()[['sysname']]=="Windows"){ # clean the memory!
      stopCluster(cl)
    } 
  }
  # Assessment of Applicability Domain (plot leverage) was included as details in each regression function
  
  # Print total execution time
  cat("\nRRegrs total execution time\n")
  print(proc.time() - ptmTot) # print running time

  #----------------------------
  # Indicate main result files
  #----------------------------
  cat("\nMAIN RESULT FILES\n")
  cat("======================\n")
  cat("Statistics for all data set splittings/methods/CV types:", ResBySplits,"\n")
  cat("Averages by method/CV type:",ResAvgsF,"\n")
  cat("Best model statistics:",ResBestF,"\n")
  cat("Best model plots:",paste(ResBestF,".repeatedcv.split",i,".pdf",sep=""),"\n")
  if (noYrand > 0) cat("Best model Y-randomization plot:",paste(ResBestF,".Yrand.Hist.pdf",sep=""),"\n")
  cat("\n* if you choose Details, additional CSV and PDF files will be create for each method.\n")

  return(list(BestMethod=best.reg,BestStats=my.stats.reg, Models=dfMod))
  # return a list with 3 items: the name of the best method, the statistics for the best model, the list with all the fitted models (including the best one)
}
