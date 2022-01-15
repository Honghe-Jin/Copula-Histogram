## This is the implementation of the proposed Copula method on histograms
## Honghe Jin
## 01/10/2022
####
library(mvtnorm)
library(tidyverse)
library(mvQuad)
library(gplots)
library(copula)
#Source functions
source('~/Documents/researches/Symbolic_new/Rcode/Copula_functions.R')

####Test Copulas####
##Gaussian
set.seed(123)
n = 30
list_hist = list()
rho0 = .1
for(i in 1:n){ #generate histograms
  x = rmvnorm(100, mean = c(0,0), sigma = matrix(c(1,rho0,rho0,1),2))
  h = hist2d(x, nbins = c(10,10), show = F)
  list_hist[[i]] = h
}

result = copula_MLE(list_hist = list_hist, Copula_name = 'Gaussian', para_range = c(-.9,.9))
result #result$mle is the estimator of rho

covariance_copula(h_list = list_hist, Copula_name = 'Gaussian', theta = result$mle)

##Clayton
n = 30
theta0 = 2

list_hist = list()
for(i in 1:n){
  cl3 <- claytonCopula(param = theta0, dim = 2)
  x <- rCopula(100, cl3)
  h = hist2d(x, nbins = c(10,10), show = F)
  list_hist[[i]] = h
}

result = copula_MLE(list_hist = list_hist, Copula_name = 'Clayton', para_range = c(0, 5))
result #result$mle is the estimator of theta

result = copula_MLE_list(list_hist = list_hist, Copula_name = 'Clayton', para_range = c(0, 5))
result #result$mle is the estimator of theta

covariance_copula(h_list = list_hist, Copula_name = 'Clayton', theta = result$mle)

##Frank
n = 30
theta0 = 3

list_hist = list()
for(i in 1:n){
  cl3 <- frankCopula(param = theta0, dim = 2)
  x <- rCopula(100, cl3)
  h = hist2d(x, nbins = c(10,10), show = F)
  list_hist[[i]] = h
}

result = copula_MLE_list(list_hist = list_hist, Copula_name = 'Frank', para_range = c(0, 5))
result

covariance_copula(h_list = list_hist, Copula_name = 'Frank', theta = result$mle)

##Gumbel
n = 20
theta0 = 3

list_hist = list()
for(i in 1:n){
  cl3 <- gumbelCopula(param = theta0, dim = 2)
  x <- rCopula(1000, cl3)
  h = hist2d(x, nbins = c(10,10), show = F)
  list_hist[[i]] = h
}

result = copula_MLE_list(list_hist = list_hist, Copula_name = 'Gumbel', para_range = c(1, 5))
result

covariance_copula(h_list = list_hist, Copula_name = 'Gumbel', theta = result$mle)

##AMH
n = 20
theta0 = .3

list_hist = list()
for(i in 1:n){
  cl3 <- amhCopula(param = theta0, dim = 2)
  x <- rCopula(1000, cl3)
  h = hist2d(x, nbins = c(10,10), show = F)
  list_hist[[i]] = h
}

result = copula_MLE(list_hist = list_hist, Copula_name = 'AMH', para_range = c(0,.9))
result

covariance_copula(h_list = list_hist, Copula_name = 'AMH', theta = result$mle)


####Correlation####
#use the Gaussian Distribution to compare the correlation
set.seed(123)
n = 30
list_hist = list()
rho0 = .1
for(i in 1:n){ #generate histograms
  x = rmvnorm(100, mean = c(0,0), sigma = matrix(c(1,rho0,rho0,1),2))
  h = hist2d(x, nbins = c(10,10), show = F)
  list_hist[[i]] = h
}

#Gaussian
#result = copula_MLE_list(list_hist = list_hist, Copula_name = 'Gaussian', para_range = c(-.9,.9))
result = copula_MLE(list_hist = list_hist, Copula_name = 'Gaussian', para_range = c(-.9,.9))
print(paste0("The estimated parameter for Gaussian copula is: ", result$mle))

bic = log(n) -2*result$max_likelihood
print(paste0("BIC for the Gaussian copula is: ", bic))

cov_result = covariance_copula(h_list = list_hist, Copula_name = 'Gaussian', theta = result$mle)
print(paste0("correlation estimation for the Gaussian copula is: ", cov_result$correlation))

#Clayton
result = copula_MLE(list_hist = list_hist, Copula_name = 'Clayton', para_range = c(0, 5))
print(paste0("The estimated parameter for Clayton copula is: ", result$mle))

bic = log(n) -2*result$max_likelihood
print(paste0("BIC for the Clayton copula is: ", bic))

cov_result = covariance_copula(h_list = list_hist, Copula_name = 'Clayton', theta = result$mle)
print(paste0("correlation estimation for the Clayton copula is: ", cov_result$correlation))

#Frank
result = copula_MLE(list_hist = list_hist, Copula_name = 'Frank', para_range = c(0, 5))
print(paste0("The estimated parameter for Frank copula is: ", result$mle))

bic = log(n) -2*result$max_likelihood
print(paste0("BIC for the Frank copula is: ", bic))

cov_result = covariance_copula(h_list = list_hist, Copula_name = 'Frank', theta = result$mle)
print(paste0("correlation estimation for the Frank copula is: ", cov_result$correlation))

#Gumbel
result = copula_MLE(list_hist = list_hist, Copula_name = 'Gumbel', para_range = c(1, 5))
print(paste0("The estimated parameter for Gumbel copula is: ", result$mle))

bic = log(n) -2*result$max_likelihood
print(paste0("BIC for the Gumbel copula is: ", bic))

cov_result = covariance_copula(h_list = list_hist, Copula_name = 'Gumbel', theta = result$mle)
print(paste0("correlation estimation for the Gumbel copula is: ", cov_result$correlation))

#AMH
result = copula_MLE(list_hist = list_hist, Copula_name = 'AMH', para_range = c(0,.9))
print(paste0("The estimated parameter for AMH copula is: ", result$mle))

bic = log(n) -2*result$max_likelihood
print(paste0("BIC for the AMH copula is: ", bic))

cov_result = covariance_copula(h_list = list_hist, Copula_name = 'AMH', theta = result$mle)
print(paste0("correlation estimation for the AMH copula is: ", cov_result$correlation))


