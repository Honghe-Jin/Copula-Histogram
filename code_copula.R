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
n = 150
list_hist = list()
rho0 = .1
for(i in 1:n){ #generate histograms
  x = rmvnorm(1000, mean = c(0,0), sigma = matrix(c(1,rho0,rho0,1),2))
  h = hist2d(x, nbins = c(20,20), show = F)
  list_hist[[i]] = h
}

result = copula_MLE(list_hist = list_hist, Copula_name = 'Gaussian', para_range = c(-.9,.9))
result #result$mle is the estimator of rho


##Clayton
n = 30
theta0 = 2

list_hist = list()
for(i in 1:n){
  cl3 <- claytonCopula(param = theta0, dim = 2)
  x <- rCopula(1000, cl3)
  h = hist2d(x, nbins = c(10,10), show = F)
  list_hist[[i]] = h
}

result = copula_MLE(list_hist = list_hist, Copula_name = 'Clayton', para_range = c(0, 5))
result #result$mle is the estimator of theta

##Frank
n = 30
theta0 = 3

list_hist = list()
for(i in 1:n){
  cl3 <- frankCopula(param = theta0, dim = 2)
  x <- rCopula(1000, cl3)
  h = hist2d(x, nbins = c(10,10), show = F)
  list_hist[[i]] = h
}

result = copula_MLE(list_hist = list_hist, Copula_name = 'Frank', para_range = c(0, 5))
result

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

result = copula_MLE(list_hist = list_hist, Copula_name = 'Gumbel', para_range = c(0, 5))
result

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

