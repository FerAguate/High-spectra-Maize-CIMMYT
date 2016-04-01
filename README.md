# USE OF HIGH-RESOLUTION IMAGE DATA OUTPERFORMS VEGETATION INDICES IN PREDICTION OF MAIZE YIELD

## Supplementary methods

To load the R objects that contains the data, follows:
```R
load('objects.rdata')
ls()
```
######[1] "CWMI" "mND"  "NDVI" "PRI"  "X1"   "X2"   "X3"   "X4"   "X5"   "Y" 
Where CWMI, mND, NDVI and PRI are matrix class objects for the indices by time points (in columns). X1 to X5 are matrices with 62 reflectance bands (in columns) and 1231 observations. Object Y is a data frame that contains grain yield observations and the ID for trials.
To create observations an trial vectors, run the next script:
```R
y = scale(Y[,1])
trial = as.numeric(Y[,2])
```
Define time point as:
```R
tP = 5
```
Fitting models on training data and obtain predictions matrix (YHat).
```R
YHat = matrix(nrow = length(y), ncol = 4)
colnames(YHat) = c('NDVI', 'OLS', 'PC', 'BB')
X = eval(parse(text = paste0('X', tP)))
SVD = svd(X, nu = ncol(X), nv = 0)
PC = SVD$u %*% diag(SVD$d)[,1:5]
  
model_NDVI = lm(y ~ NDVI[,tP])
YHat[,'NDVI'] = cbind(1,NDVI[,tP]) %*% coef(model_NDVI)
model_OLS = lm(y ~ X) 
YHat[,'OLS'] = cbind(1,X) %*% coef(model_OLS)
model_PC = lm(y ~ PC) 
YHat[,'PC'] = cbind(1, PC) %*% coef(model_PC)
library(BGLR)
model_BB = BGLR(y = y, ETA = list(list(X = X, model = 'BayesB')), nIter = 100000, burnIn = 3000)
YHat[,'BB'] = X %*% model_BB$ETA[[1]]$b + model_BB$mu
```
Calculate across-trial correlations in training data.
```R
cor(y, YHat)
```
Calculate within-trial correlations in training data.
```R
do.call('rbind', lapply(by(cbind(y, YHat), Y[,2], I), function(x) cor(x[,1], x[,2:5])))
```
Obtain leave-one-trial out predictions.
```R
YHat_cv = matrix(nrow = length(y), ncol=4)
colnames(YHat_cv) = c('NDVI', 'OLS', 'PC', 'BB')

for (i in unique(trial)){
  trn = which(trial != i) ; tst = which(trial == i)
  yTRN = y[trn]           ; yTST = y[tst]
  PC.TRN = PC[trn,]       ; PC.TST = PC[tst,] 
  X.TRN = X[trn,]         ; X.TST = X[tst,]
  NDVI.TRN = NDVI[trn,tP] ; NDVI.TST=NDVI[tst,tP]

  model_NDVI_cv = lm(yTRN ~ NDVI.TRN)
  YHat_cv[tst,'NDVI'] = cbind(1, NDVI.TST) %*% coef(model_NDVI_cv)
  model_OLS_cv = lm(yTRN ~ X.TRN) 
  YHat_cv[tst,'OLS'] = cbind(1, X.TST) %*% coef(model_OLS_cv)
  model_PC_cv = lm(yTRN ~ PC.TRN) 
  YHat_cv[tst,'PC'] = cbind(1, PC.TST) %*% coef(model_PC_cv)
  model_BB_cv = BGLR(y = yTRN, ETA = list(list(X = X.TRN, model = 'BayesB')), nIter = 100000, burnIn = 3000)
  YHat_cv[tst,'BB'] = X.TST %*% model_BB_cv$ETA[[1]]$b + model_BB_cv$mu
}
```
Define multi-time points:
```R
mtP=4:5
```
Fitting models with multi-time points on the training data.
```R
YHat_mtp = data.frame(matrix(nrow=length(y),ncol=4))
colnames(YHat_mtp)=c('NDVI','OLS','PC','BB')

Xtp = list()
for (tp in mtP){Xtp[[which(tp == mtP)]] = eval(parse(text = paste0('X', tp)))}
Xmtp = do.call('cbind',Xtp)
SVD = svd(Xjoin, nu = ncol(Xjoin), nv = 0)
PCmtp = SVD$u %*% diag(SVD$d)[,1:5]
  
model_NDVI_mtP = lm(y ~ NDVI[,mtP])
YHat_mtp[,'NDVI'] = cbind(1,NDVI[,mtP]) %*% coef(model_NDVI_mtP)
model_OLS_mtP = lm(y ~ Xmtp)
YHat_mtp[,'OLS'] = cbind(1, Xmtp) %*% coef(model_OLS_mtP)
model_PC_mtP = lm(y ~ PCmtp)
YHat_mtp[,'PC'] = cbind(1, PCmtp) %*% coef(model_PC_mtP)
  
ETA=list()
for (L in order(timePoint)){ETA[[L]]= list(X=Xtp[[L]],model=model)}
model_BB_mtp = BGLR(y = y, ETA = ETA, nIter = 100000, burnIn = 3000)
YHat[,'BB'] = cbind(1, Xmtp) %*% do.call('c', lapply(order(mtP), function(x) model_BB_mtp$ETA[[x]]$b)) + model_BB_mtp$mu
```

