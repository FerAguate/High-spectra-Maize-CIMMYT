# USE OF HIGH-RESOLUTION IMAGE DATA OUTPERFORMS VEGETATION INDICES IN PREDICTION OF MAIZE YIELD

## Supplementary methods

Load and prepare data
```R
load('objects.rdata')
y = scale(Y[,1])
trial = as.numeric(Y[,2])
```
Define time point
```R
tP = 5
```
Fitting models on training data and obtain predictions matrix
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
Calculate Across trial correlations in training data
```R
cor(y, YHat)
```
Calculate within trial correlations in training data
```R
do.call('rbind', lapply(by(cbind(y, YHat), Y[,2], I), function(x) cor(x[,1], x[,2:5])))
```
Obtain leave-one-trial out predictions
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
