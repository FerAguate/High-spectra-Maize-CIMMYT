## USE OF HIGH-RESOLUTION IMAGE DATA OUTPERFORMS VEGETATION INDICES IN PREDICTION OF MAIZE YIELD
### Supplementary methods

#### Number of LVs for PLS model
```R
# Load pls library
library(pls)

# Loading data from heat and Drought Stress trials only
load('objectsDS.rdata')

# Create an object for scaled grain yield
y = scale(Y[,1])

# Initialize the list that will contain the number of LVs
nLV.List=list()

#### The next will run for each of the single and multi-time points
for (w in 1:9){
  if(w < 6)  {X = eval(parse(text=paste0("X",w)))}
  if(w == 6) {X = cbind(X4,X5)}
  if(w == 7) {X = cbind(X3,X4,X5)}
  if(w == 8) {X = cbind(X2,X3,X4,X5)}
  if(w == 9) {X = cbind(X1,X2,X3,X4,X5)}
  
  nFolds = length(unique(Y[,2]))
  nLV = seq(from=2,to=15,by=1)
  folds = as.numeric(Y[,2])
  COR.TST = rep(NA,nFolds)
  nLV.BEST = rep(NA,nFolds)
  
  for(i in unique(folds)){
    trn = which(folds!=i)
    COR = matrix(nrow=11,ncol=length(nLV))
    
    for(j in 1:11){
      trn_trn = which(folds[trn]!=j)
      for(k in 1:length(nLV)){
        COR[j,k] = cor(y[-trn_trn],predict(plsr(formula=y[trn_trn]~X[trn_trn,],ncomp=nLV[k]),ncomp=nLV[k],newdata=X[-trn_trn,]))
      }
      print(paste("subFold=",j,"Fold=",i,"Timepoint=",w))
    }
    
    COR = colMeans(COR)
    nLV.BEST[i] = nLV[which(COR==min(COR))]
    COR.TST[i] = cor(y[-trn],predict(plsr(formula=y[trn]~X[trn,],ncomp=nLV.BEST[i]),ncomp=nLV.BEST[i],newdata=X[-trn,]))
  }
  nLV.List[[w]] = nLV.BEST
}
```
#### The following script performs leave-one-trial out cross validation to obtain predictions for each model with data from heat and drought trials. Predictions of the irrigated trial were obtained separately with an adaptated version of this code.

```R
# Load libraries
library(BGLR)

# Create an object for identification of trials
trials = as.integer(factor(Y[,2]))

# Obtain predictions from single time point data
for (k in 1:5){
  # initialize the object that will contain predictions
  YHat = data.frame(matrix(nrow=length(y),ncol=9))
  colnames(YHat)=c('NDVI','CWMI','mND','PRI','MTCI','MCARI2','OLS','PLS','BayesB')
  X = eval(parse(text=paste0('X',k)))
  X = X/sqrt(ncol(X))
  # Cross Validation leave one trial out
  for (i in unique(trials)){
    print(i)
    trn = which(trials!=i); tst = which(trials==i)
    yTRN = y[trn] ; yTST = y[tst]
    trialTRN = trials[trn] ; trialTST = trials[tst] 
    X.TRN = X[trn,]   ; X.TST = X[tst,]
    NDVI.TRN = NDVI[trn,k] ; NDVI.TST = NDVI[tst,k]
    CWMI.TRN = CWMI[trn,k] ; CWMI.TST = CWMI[tst,k]   
    mND.TRN = mND[trn,k] ; mND.TST = mND[tst,k] 
    PRI.TRN = PRI[trn,k] ; PRI.TST = PRI[tst,k]
    MTCI.TRN = MTCI[trn,k] ; MTCI.TST = MTCI[tst,k]
    MCARI2.TRN=MCARI2[trn,k] ; MCARI2.TST = MCARI2[tst,k]    
    # NDVI
    Model1 = lm(yTRN ~ NDVI.TRN)
    pY1 = as.matrix(NDVI.TST)%*%(coef(Model1)[-1]) + coef(Model1)[1]
    YHat[tst,'NDVI'] = pY1
    # CWMI
    Model2 = lm(yTRN ~ CWMI.TRN)
    pY2 = as.matrix(CWMI.TST)%*%(coef(Model2)[-1]) + coef(Model2)[1]
    YHat[tst,'CWMI'] = pY2
    # mND
    Model3 = lm(yTRN ~ mND.TRN)
    pY3 = as.matrix(mND.TST)%*%(coef(Model3)[-1]) + coef(Model3)[1]
    YHat[tst,'mND'] = pY3
    # PRI
    Model4 = lm(yTRN ~  PRI.TRN)
    pY4 = as.matrix(PRI.TST)%*%(coef(Model4)[-1]) + coef(Model4)[1]
    YHat[tst,'PRI'] = pY4
    # MTCI
    Model5 = lm(yTRN ~ MTCI.TRN)
    pY5 = as.matrix(MTCI.TST)%*%(coef(Model5)[-1]) + coef(Model5)[1]
    YHat[tst,'MTCI'] = pY5
    # MCARI2
    Model6 = lm(yTRN ~ MCARI2.TRN)
    pY6 = as.matrix(MCARI2.TST)%*%(coef(Model6)[-1]) + coef(Model6)[1]
    YHat[tst,'MCARI2'] = pY6
    # OLS
    Model7 = lm(yTRN ~ X.TRN) 
    pY7 = as.matrix(X.TST)%*%(coef(Model7)[-1]) + coef(Model7)[1]
    YHat[tst,'OLS'] = pY7
    # PLS
    Model8 = plsr(formula = yTRN~X.TRN, ncomp = nLV.List[[k]][i])
    pY8 = predict(Model8,ncomp=nLV.List[[k]][i],newdata=X.TST)
    YHat[tst,'PLS'] = pY8
    # BayesB
    Model9 = BGLR(y=yTRN,ETA=list(list(X=X.TRN,model='BayesB')),nIter=120000,burnIn=5000,      verbose=F,rmExistingFiles=T,saveAt="~/tmp_uavs",groups=trials[trn])
    pY9 = as.matrix(X.TST)%*%Model9$ETA[[1]]$b + Model9$mu
    YHat[tst,'BayesB'] = pY9
  }
  # Save the predictions in files separated by time point  
  save(y,YHat,trials,file=paste0('YHATDSCV',timePoint,'.RData'))
}
```
#### Predictions with multi-time points data
```R
for (z in 1:4){
  timePoint = c(z:5)
  YHat = data.frame(matrix(nrow=length(y),ncol=9))
  colnames(YHat) = c('NDVI','CWMI','mND','PRI','MTCI','MCARI2','OLS','PLS','BayesB')
  X = list()
  for (i in timePoint) {
    X[[which(timePoint==i)]] = eval(parse(text=paste0('X',i)))
    X[[which(timePoint==i)]] = X[[which(timePoint==i)]]/sqrt(ncol(X[[which(timePoint==i)]]))
  }
  # Cross Validation leave one trial out
  for (i in unique(trials)){
    print(i)
    trn = which(trials!=i) ; tst = which(trials==i)
    yTRN = y[trn] ; yTST = y[tst]
    trialTRN = trials[trn] ; trialTST = trials[tst] 
    X.TRN = do.call('cbind',X)[trn,] ; X.TST = do.call('cbind',X)[tst,]
    NDVI.TRN = NDVI[trn,timePoint] ; NDVI.TST = NDVI[tst,timePoint]
    CWMI.TRN = CWMI[trn,timePoint] ; CWMI.TST = CWMI[tst,timePoint]   
    mND.TRN = mND[trn,timePoint] ; mND.TST = mND[tst,timePoint] 
    PRI.TRN = PRI[trn,timePoint] ; PRI.TST = PRI[tst,timePoint]
    MTCI.TRN = MTCI[trn,timePoint] ; MTCI.TST = MTCI[tst,timePoint]
    MCARI2.TRN = MCARI2[trn,timePoint] ; MCARI2.TST = MCARI2[tst,timePoint]
    # NDVI
    Model = lm(yTRN ~ NDVI.TRN)
    pY = as.matrix(NDVI.TST)%*%(coef(Model)[-1]) + coef(Model)[1]
    YHat[tst,'NDVI'] = pY
    # CWMI
    Model2 = lm(yTRN ~ CWMI.TRN)
    pY2 = as.matrix(CWMI.TST)%*%(coef(Model2)[-1]) + coef(Model2)[1]
    YHat[tst,'CWMI'] = pY2
    # mND
    Model3 = lm(yTRN ~ mND.TRN)
    pY3 = as.matrix(mND.TST)%*%(coef(Model3)[-1]) + coef(Model3)[1]
    YHat[tst,'mND'] = pY3
    # PRI
    Model4 = lm(yTRN ~ PRI.TRN)
    pY4 = as.matrix(PRI.TST)%*%(coef(Model4)[-1]) + coef(Model4)[1]
    YHat[tst,'PRI'] = pY4
    # MTCI
    Model5 = lm(yTRN ~ MTCI.TRN)
    pY5 = as.matrix(MTCI.TST)%*%(coef(Model5)[-1]) + coef(Model5)[1]
    YHat[tst,'MTCI'] = pY5
    # MCARI2
    Model6 = lm(yTRN ~ MCARI2.TRN)
    pY6 = as.matrix(MCARI2.TST)%*%(coef(Model6)[-1]) + coef(Model6)[1]
    YHat[tst,'MCARI2'] = pY6
    # OLS
    Model7 = lm(yTRN ~ X.TRN) 
    pY7 = as.matrix(X.TST)%*%(coef(Model7)[-1]) + coef(Model7)[1]
    YHat[tst,'OLS'] = pY7
    # PLS
    Model8 = plsr(formula=yTRN~X.TRN,ncomp=nLV.List[[10-z]][i])
    pY8 = predict(Model8,ncomp=nLV.List[[10-z]][i],newdata=X.TST)
    YHat[tst,'PLS'] = pY8
    # BayesB
    ETA = list()
    for (L in order(timePoint)){
      ETA[[L]] = list(X=X[[L]][trn,],model='BayesB')
    }
    Model9 = BGLR(y=yTRN,ETA=ETA,nIter=120000,burnIn=5000,verbose=F,rmExistingFiles=T,saveAt="~/tmp_uavs")
    pY9 = as.matrix(X.TST)%*%do.call('c',lapply(order(timePoint),function(x)Model9$ETA[[x]]$b)) + Model9$mu
    YHat[tst,'BayesB'] = pY9
  }
  #Save YHat separately
  save(y,YHat,trials,file=paste0('AnalyzingData/output/YHATDSCV',paste0(timePoint,collapse = '-'),'.RData'))
}
```
#### Bootstrap correlations and SE

```R
#number of replicates
nRep= 5000

#Load and bind predictions with heat and drought stress data and  irrigated data for T5.
load('YHATDSCV5.rdata')
load('YHATWWCV5.rdata')
YHAT = rbind(YHat,YHatWW)
y = rbind(y,yWW)
trials = rbind(trials,trialsWW) #trial 12 is the irrigated one

#Initialize objects for weighted within-trial correlations
WW.WWTC = matrix(nrow=nRep,ncol=ncol(YHAT)) #Well-watered conditions
DS.WWTC = matrix(nrow=nRep,ncol=ncol(YHAT)) #Drought stress

#Run bootstrap correlations for each model
for (index in 1:ncol(YHAT)){
  COR = matrix(nrow=nRep,ncol=length(unique(trials)))
  SV = matrix(nrow=nRep,ncol=length(unique(trials)))
  for(trial in 1:length(unique(trials))){
    n = sum(trials==trial)
    print(paste('index=',index,'trial=',trial))
    for(rep in 1:nRep){
      tmp = sample(1:n,size=n,replace=TRUE)
      COR[rep,trial] = cor(y[trials==trial][tmp],YHAT[trials==trial,index][tmp])
      SV[rep,trial] = (1-COR[rep,trial]^2)/(n-2)
    }
  }
  for(rep in 1:nRep){
    DS.WWTC[rep,index] = sum(COR[rep,1:11]/SV[rep,1:11])/sum(1/SV[rep,1:11])
    WW.WWTC[rep,index] = sum(COR[rep,12]/SV[rep,12])/sum(1/SV[rep,12])
  }
}
colnames(WW.WWTC) <- colnames(DS.WWTC) <- colnames(YHAT)
```
#### Table with results
```R
PropMatrix = matrix(ncol=ncol(YHAT),nrow=ncol(YHAT))
for(i in 1:ncol(YHAT)){ 
  for(j in 1:ncol(YHAT)){
    if(i>j) PropMatrix[i,j] = round( mean(DS.WWTC[,j]<DS.WWTC[,i]),2)
    if(i<j) PropMatrix[i,j] = round( mean(WW.WWTC[,j]<WW.WWTC[,i]),2)
  }
}
colnames(PropMatrix) <- rownames(PropMatrix) <- colnames(YHAT)
cbind('Drought'=round(apply(DS.WWTC,2,mean),2),
      'SEDro'=round(apply(DS.WWTC,2,sd),3),
      'Irrigated'=round(apply(WW.WWTC,2,mean),2),
      'SEIrr'=round(apply(WW.WWTC,2,sd),3),
      PropMatrix)
```
