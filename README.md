# USE OF HIGH-RESOLUTION IMAGE DATA OUTPERFORMS VEGETATION INDICES IN PREDICTION OF MAIZE YIELD

## Supplementary methods

To load the R objects that contains the data, follows:
```R
load('objects.rdata')
ls()
```
######[1] "CWMI" "MCARI2" "mND" "MTCI" "NDVI" "PRI"  "X1"   "X2"   "X3"   "X4"   "X5"   "Y" 
Where CWMI, MCARI2, mND, MTCI, NDVI and PRI are matrix class objects for the indices with time point 1 to 5 in columns. X1 to X5 are matrices with the 62 reflectance bands in columns, and the 1231 plot observations in rows. Object Y is a data frame that contains two columns, the first has grain yield observations and the second ID for trials.
To create the observations vector and a vector for trial identifications, execute the following:
```R
y = scale(Y[,1])
trial = as.numeric(Y[,2])
```
To define a time point to start working with, create a "tp" object as:
```R
tP = 5
```
The next script Fits models on training data and obtain predictions matrix (YHat). For PLS model, the number of latent variables was previously defined in the object "nLV".
```R
YHat = matrix(nrow = length(y), ncol = 4)
colnames(YHat) = c('NDVI', 'OLS', 'PLS', 'BB')
X = eval(parse(text = paste0('X', tP)))

model_NDVI = lm(y ~ NDVI[,tP])
YHat[,'NDVI'] = cbind(1,NDVI[,tP]) %*% coef(model_NDVI)
model_OLS = lm(y ~ X) 
YHat[,'OLS'] = cbind(1,X) %*% coef(model_OLS)
library(pls)  
YHat[,'PLS'] = predict(plsr(formula = y~X, ncomp = nLV), ncomp = nLV)
library(BGLR)
model_BB = BGLR(y = y, ETA = list(list(X = X, model = 'BayesB')), nIter = 100000, burnIn = 3000)
YHat[,'BB'] = X %*% model_BB$ETA[[1]]$b + model_BB$mu
```
The next stands to obtain across-trial correlations in training data.
```R
cor(y, YHat)
```
For within-trial correlations in training data.
```R
do.call('rbind', lapply(by(cbind(y, YHat), Y[,2], I), function(x) cor(x[,1], x[,2:5])))
```
Leave-one-trial out predictions are obtain by fitting models with training trials and predicting in each of the selected testing trial.
```R
YHat_cv = matrix(nrow = length(y), ncol=4)
colnames(YHat_cv) = c('NDVI', 'OLS', 'PLS', 'BB')

for (i in unique(trial)){
  trn = which(trial != i) ; tst = which(trial == i)
  yTRN = y[trn]           ; yTST = y[tst]
  X.TRN = X[trn,]         ; X.TST = X[tst,]
  NDVI.TRN = NDVI[trn,tP] ; NDVI.TST=NDVI[tst,tP]

  model_NDVI_cv = lm(yTRN ~ NDVI.TRN)
  YHat_cv[tst,'NDVI'] = cbind(1, NDVI.TST) %*% coef(model_NDVI_cv)
  model_OLS_cv = lm(yTRN ~ X.TRN) 
  YHat_cv[tst,'OLS'] = cbind(1, X.TST) %*% coef(model_OLS_cv)
  model_PLS_cv = plsr(formula = yTRN ~ X.TRN, ncomp = cvnLV[i]) 
  YHat_cv[tst,'PLS'] = predict(model_PLS_cv, ncomp = cvnLV[i], newdata = X.TST)
  model_BB_cv = BGLR(y = yTRN, ETA = list(list(X = X.TRN, model = 'BayesB')), nIter = 100000, burnIn = 3000)
  YHat_cv[tst,'BB'] = X.TST %*% model_BB_cv$ETA[[1]]$b + model_BB_cv$mu
}
```
Set multi-time points by creating a numeric vector:
```R
mtP=4:5
```
This code will Fit models with multi-time point information on the training data to evaluate godness of fit.
```R
YHat_mtp = data.frame(matrix(nrow=length(y),ncol=4))
colnames(YHat_mtp)=c('NDVI','OLS','PLS','BB')

Xtp = list()
for (tp in mtP){Xtp[[which(tp == mtP)]] = eval(parse(text = paste0('X', tp)))}
Xmtp = do.call('cbind',Xtp)

model_NDVI_mtP = lm(y ~ NDVI[,mtP])
YHat_mtp[,'NDVI'] = cbind(1,NDVI[,mtP]) %*% coef(model_NDVI_mtP)
model_OLS_mtP = lm(y ~ Xmtp)
YHat_mtp[,'OLS'] = cbind(1, Xmtp) %*% coef(model_OLS_mtP)
YHat_mtp[,'PLS'] = predict(plsr(formula = y~Xmtp, ncomp = nLVmtp), ncomp = nLVmtp)
  
ETA=list()
for (L in order(timePoint)){ETA[[L]]= list(X=Xtp[[L]],model=model)}
model_BB_mtp = BGLR(y = y, ETA = ETA, nIter = 100000, burnIn = 3000)
YHat[,'BB'] = cbind(1, Xmtp) %*% do.call('c', lapply(order(mtP), function(x) model_BB_mtp$ETA[[x]]$b)) + model_BB_mtp$mu
```

