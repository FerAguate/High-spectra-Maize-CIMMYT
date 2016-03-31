# USE OF HIGH-RESOLUTION IMAGE DATA OUTPERFORMS VEGETATION INDICES IN PREDICTION OF MAIZE YIELD

## Supplementary methods

Load and prepare data
```R
load('objects.rdata')
y = scale(Y[,1])
```
Define time point
```R
timePoint = 5
```
Fitting models on training data and obtain predictions matrix
```R
YHat = matrix(nrow = length(y), ncol = 4)
colnames(YHat) = c('NDVI', 'OLS', 'PC', 'BB')
X = eval(parse(text = paste0('X', timePoint)))
SVD = svd(X, nu = ncol(X), nv = 0)
PC = SVD$u %*% diag(SVD$d)
  
Model_NDVI = lm(y ~ NDVI[,timePoint])
YHat[,'NDVI'] = cbind(1,NDVI[,timePoint]) %*% coef(Model_NDVI)
Model_OLS = lm(y ~ X) 
YHat[,'OLS'] = cbind(1,X) %*% coef(Model_OLS)
Model_PC = lm(y ~ PC[,1:5]) 
YHat[,'PC'] = cbind(1,PC[,1:5]) %*% coef(Model_PC)
library(BGLR)
Model_BB = BGLR(y = y, ETA = list(list(X = X, model = 'BayesB')), nIter = 103000, burnIn = 3000)
YHat[,'BB'] = X %*% Model_BB$ETA[[1]]$b + Model_BB$mu
```
Calculate Across trial correlations in training data
```R
cor(y,YHat)
```
Calculate within trial correlations in training data
```R
do.call('rbind', lapply(by(cbind(y,YHat), Y[,2], I), function(x) cor(x[,1], x[,2:5])))
```
