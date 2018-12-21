\name{RRegrs}
\alias{RRegrs}
\title{
Fitting a set of ten fully utilized Predictive Models.
}
\description{
The current tool is a collection of regression tools from R that could be used to search the best regression models for any dataset. Filtering, Cross-Validation, Y-randomization are applied.
}
\usage{
RRegrs(DataFileName="ds.House.csv",DataFileSep=",",DependentVariable=NULL,IndependentVariables=NULL,
PathDataSet="DataResults", noCores=1, ResAvgs="RRegsResAvgs.csv", ResBySplits="RRegrsResAllSplits.csv", ResBest="RRegrsResBest.csv", fDet="T", fFilters="F", fScaling="T", fRemNear0Var="T", fRemCorr="T",fLM="T", fGLM="T", fPLS="T", fLASSO="T", fSVRM="T", fNN="T", fRF="T", fRFRFE="T", fSVMRFE="T", fENET="T", RFE_SVM_C="1;5;10;15;20", RFE_SVM_epsilon="0.01;0.1;0.3", cutoff=0.9, iScaling=1, iScalCol=1, trainFrac=0.75, iSplitTimes=10, noYrand=100, CVtypes="repeatedcv;LOOCV", NoNAValFile="ds.NoNA.csv",  No0NearVarFile="ds.No0Var.csv", ScaledFile="ds.scaled.csv", NoCorrFile="ds.scaled.NoCorrs.csv",
lmFile="LM.details.csv", glmFile="GLM.details.csv", plsFile="PLS.details.csv",lassoFile="Lasso.details.csv", svrmFile="SVMRadial.details.csv", nnFile="NN.details.csv", rfFile="RF.details.csv", rfrfeFile="RFRFE.details.csv", svmrfeFile="SVMRFE.details.csv",enetFile="ENET.details.csv",fR2rule="T")
}
\arguments{
 \item{DataFileName}{ a character string naming the data set file; the data set should be an object where samples are in rows and features are in columns. The first column should be a numeric or factor vector containing the outcome for each sample. Default value is the "ds.House.csv" data.}
 \item{DataFileSep}{ the column separator character in the data file. Default value is a comma (CSV file).}
 \item{DependentVariable}{ Default value is NULL, causing the first column to be used as dependent variable.}
 \item{IndependentVariables}{ Default value is NULL, causing the second to last columns to be used as dependent variables.}
 \item{PathDataSet}{ a character string naming the directory where the data set file is stored. Default value is the "DataResults" directory.}
 \item{noCores}{The number of cores used for calculations. Default value is 1.}
 \item{ResAvgs}{a character string naming the file where average  values of all models are stored. Default value is the "RRegsResAvgs.csv" file.}
 \item{ResBySplits}{a character string naming the file where statistics for all splits and all models are stored. Default value is the "RRegrsResAllSplits.csv" file.}
 \item{ResBest}{a character string naming the file where statistics for the best model are stored. Default value is the "RRegrsResBest.csv" file.}
 \item{fDet}{logical. If TRUE all details for all splits and models are stored.}
 \item{fFilters}{If run custom flters. Not implemented yet.}
 \item{fScaling}{logical. If scaling will be performed. Default value is TRUE.}
 \item{fRemNear0Var}{logical. If run Removal of near zero variance columns. Default value is TRUE.}
 \item{fRemCorr}{logical. If correlated features will be removed based on 'cutoff' value. Default value is TRUE.}
 \item{fLM}{logical. If Linear model will be applied to the data. Default value is TRUE.}
 \item{fGLM}{logical. If Generalized Linear model will be applied to the data. Default value is TRUE.}
 \item{fPLS}{logical. If Partial Leasts Square regression model will be applied to the data. Default value is TRUE.}
 \item{fLASSO}{logical. If Lasso regression model will be applied to the data. Default value is TRUE.}
 \item{fSVRM}{logical. If Support vector machine regression model (radial basis function kernel) will be applied to the data. Default value is TRUE.}
 \item{fNN}{logical. If Neural network model will be applied to the data. Default value is TRUE.}
 \item{fRF}{logical. If Random Forest model will be applied to the data. Default value is TRUE.}
 \item{fRFRFE}{logical. If random Forest with Random Feature Elimination model will be applied to the data. Default value is TRUE.}
 \item{fSVMRFE}{logical. If Support vector machines with Random Feature elimination model (radial basis function kernel) will be applied to the data. Default value is TRUE.}
 \item{fENET}{logical. If Elastic net model will be applied to the data. Default value is TRUE.}
 \item{RFE_SVM_C}{Support vector machines cost parameter. Default values are c(1,5,10,15,20).}
 \item{RFE_SVM_epsilon}{Support vector machines epsilon parameter. Default values are c(0.01,0.1,0.3).}
 \item{cutoff}{=0.9}{Cutoff for correlated features (default = 0.9).}
 \item{iScaling}{Type of scaling: 1 = normalization; 2 = standardization; 3 = other; any other: no scaling}
 \item{iScalCol}{Scaling columns: 1 = including dependent variable; 2: only all the features}
 \item{trainFrac}{ Fraction of training set from the entire dataset (default = 0.75); the rest of dataset is the test set.}
 \item{iSplitTimes}{Number of splittings the dataset into train and test (default = 10)}
 \item{noYrand}{Number of Y-Randomization (default = 100)}
 \item{CVtypes}{a character vector indicating which cross validation method will be used. Default value is c("repeatedcv","LOOCV")}
 \item{NoNAValFile}{a character string naming the file where the data set without NA values will be stored. Default value is the "ds.NoNA.csv" file.}
 \item{No0NearVarFile}{a character string naming the file where the data set without without near zero variance values will be stored. Default value is the "ds.No0Var.csv" file.}
 \item{ScaledFile}{a character string naming the file where the scaled data set will be stored. Default value is the "ds.scaled.csv" file.}
 \item{NoCorrFile}{a character string naming the file where the data set without correlated independent variables will be stored. Default value is the "ds.scaled.NoCorrs.csv" file.}
 \item{lmFile}{a character string naming the file where LM model details will be stored. Default value is the "LM.details.csv" file.}
 \item{glmFile}{a character string naming the file where GLM model details will be stored. Default value is the "GLM.details.csv" file.}
 \item{plsFile}{a character string naming the file where PLS model details will be stored. Default value is the "PLS.details.csv" file.}
 \item{lassoFile}{a character string naming the file where LASSO model details will be stored. Default value is the "Lasso.details.csv" file.}  
 \item{svrmFile}{a character string naming the file where SVMR model details will be stored. Default value is the "SVMRadial.details.csv" file.}
 \item{nnFile}{a character string naming the file where NN model details will be stored. Default value is the "NN.details.csv" file.}
 \item{rfFile}{a character string naming the file where RF model details will be stored. Default value is the "RF.details.csv" file.}
 \item{rfrfeFile}{a character string naming the file where RFRFE model details will be stored. Default value is the "RFRFE.details.csv" file.}
 \item{svmrfeFile}{a character string naming the file where SVMRRFE model details will be stored. Default value is the "SVMRFE.details.csv" file.}
 \item{enetFile}{a character string naming the file where ENET model details will be stored. Default value is the "ENET.details.csv" file.}
 \item{fR2rule}{logical. If R2 rule will be used to find the best regression model. Default value is the TRUE (R2 rule).}
}
\details{ RRegrs() function contains several sections: loading parameters and dataset, remove near zero variance features, scaling dataset, remove correlated features, dataset splitting, run the 10 regression methods, summary of statistics for all methods and splittings, averages for each method and cross-validation type for all splittings, automatic best model statistics, best model Y-randomization. Assessment of Applicability Domain was included in each method.

In order to use RRegrs function, it is needed to specify a minimum set of parameters (the others will get default values). All the parameters will be written into a CSV file (ex: Parameters.csv), in the same working folder where it should be present the input dataset file and the outputs files. The dataset needs to have CSV format, and the first column should be the dependent variable. The input and output files will be placed into a working folder. The output files are CSV statistics files, PDF and PNG plots.

If details are needed (fDet=T), several output files are generated. A CSV file with detailed statistics about the regression model is returned (Regression method, splitting number, cross-validation type, Training set summary, Test set summary, Fitting summary, List of predictors, Training predictors, Test predictors, Full statistics, Feature importance, Residuals of the fitted model, Assessment of applicability domain / leverage analysis). If the determinant is not zero then the following are reported: mean of hat values, hat values with warnings (X3 and X2 for values 3 and 2 times than hat mean), leverage threshold, list of points with leverage greater than threshold, Cook's distance cutoff, Cook's distances, points influence. Also, 5-12 plots for fitting statistics as a PDF file for each splitting and cross-validation method are returned. Along with the above, plots for Training Yobs versus Ypred are generated, plots for Test Yobs versus Ypred, Feature Importance, Fitted vs. Residuals for Fitted Model, Leverage for Fitted Model, Cook's Distance for Fitted Model, 6 standard fitting plots using plot function with cutoff.Cook.

The outputs can differ depending on the regression method used.

After filtering the dataset for correlated variables, near-zero variance features and splitting the dataset into training and test sets, the user's selected regression methods will be executed for each splitting and cross-validation type. 

Averaged statistics across all splittings, for each regression method and cross-validation type are calculated and used to chose the best regression model. Particularly, this is chosen based on the following criteria: from the best test R-squared (+/- 0.05), the model with minimum RMSE is the final one. If fR2rule is FALSE, adjusted R2 values will be used to choose the best model. For the best model, an additional CSV file is generated containing detailed statistics as well as PDF plots for important statistics. Finally, the best model (last data split) is tested with Y-randomization runs (100 by default).  

Please note that the Cross-validation types available are the resampling methods available from trainControl() and rfeControl() caret.
}
\value{
RRegrs() returns a list with the following 3 items:
\item{BestMethod}{the name of the best method}
\item{BestStats}{the statistics for the best model}
\item{Models}{the list with all the fitted models based on caret functions (including the best model).}
}
\examples{
\dontrun{
library(data.table)
library(corrplot)
source("../RRegrs/R/RRegrs_Functions.R")

#Run RRegrs with all default parameters
# (use default dataset file and working folder,
#  run all regression methods, without wrappers,
#  10 splitings, 100 times Y-randomization,
#  no parallel calculation = 1 CPU core)

#Run RRegrs for a specific dataset file and the rest default #parameters
RRegrsResults = RRegrs(DataFileName="MyDataSet.csv")
}
}
\author{
Georgia Tsiliki, Cristian R. Munteanu, Jose A. Seoane, Carlos Fernandez-Lozano
}
