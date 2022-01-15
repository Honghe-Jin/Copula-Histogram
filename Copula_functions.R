###
#this includes the necessary functions to fit copulas on histograms
#Honghe Jin
#01/10/2022
##
####Copulas####
# p-dimensional Gaussian Copula 
copula_Gaussian <- function(u = u, rho = NULL, Sigma = NULL){
  #Gaussian copula function
  #rho is special for bivariate Gaussian copula
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  if(is.null(rho) & is.null(Sigma)) stop('Must include rho or Sigma')
  if(is.null(Sigma)) Sigma <- matrix(c(1,rho,rho,1), nrow = 2) #bivariate
  pmvnorm(upper = qnorm(u), sigma = Sigma)[1]
}

dcopula_Gaussian <- function(u = u, rho = NULL, Sigma = NULL){
  #Gaussian copula density 
  p <- length(u)
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  if(is.null(rho) & is.null(Sigma)) stop('Must include rho or Sigma')
  if(max(u)==1|min(u)==0) return(0)
  if(is.null(Sigma)) Sigma <- matrix(c(1,rho,rho,1), nrow = 2) #bivariate
  
  exp(-.5*qnorm(u)%*%(solve(Sigma)-diag(p))%*%qnorm(u) )/sqrt(det(Sigma))
}

# Student t Copula
copula_t <- function(u = u, Sigma = Sigma, df = df){
  #t copula function
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  return(pmvt(qt(u, df = df), sigma = Sigma))
}

dcopula_t <- function(u = u, Sigma = Sigma, df = df){
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  p = length(u)
  gamma((df + p)/2)*gamma(df/2)^(p-1)*(1 + t(u)%*%solve(Sigma)%*%u)^(-(df + p)/2)/
    sqrt(det(Sigma))/gamma((df + 1)/2)^p/(prod((1 + qt(u, df = df)^2/df )))^(-(df + p)/2)
}

# Clayton Copula
copula_clayton <- function(u, theta){
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  p = length(u)
  cl_c = claytonCopula(theta, dim = p)
  pCopula(u, cl_c)
}

dcopula_clayton <- function(u, theta){
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  p = length(u)
  cl_c = claytonCopula(theta, dim = p)
  dCopula(u, cl_c)
}

# Frank Copula
copula_frank <- function(u, theta){
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  p = length(u)
  fr_c = frankCopula(theta, dim = p)
  pCopula(u, fr_c)
}

dcopula_frank <- function(u, theta){
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  p = length(u)
  fr_c = frankCopula(theta, dim = p)
  dCopula(u, fr_c)
}

# Gumbel Copula 
copula_gumbel <- function(u, theta){
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  p = length(u)
  gu_c = gumbelCopula(theta, dim = p)
  pCopula(u, gu_c)
}

dcopula_gumbel <- function(u, theta){
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  p = length(u)
  gu_c = gumbelCopula(theta, dim = p)
  dCopula(u, gu_c)
}

# AMH copula
copula_amh <- function(u, theta){
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  p = length(u)
  amh_c = amhCopula(theta, dim = p)
  pCopula(u, amh_c)
}

dcopula_amh <- function(u, theta){
  if(max(u)>1|min(u)<0) stop('element not between 0 and 1')
  p = length(u)
  amh_c = amhCopula(theta, dim = p)
  dCopula(u, amh_c)
}


####2d histogram to a function ####
h2dtofunction = function(h){
  if(class(h) != 'hist2d') stop('input must be a hist2d object')
  x_break = length(h$x.breaks)
  y_break = length(h$y.breaks)
  f_hist = function(x){
    #the input should be a 2-d vector 
    if(x[1] <= h$x.breaks[1]| x[1] >= h$x.breaks[x_break]|x[2] <= h$y.breaks[1]| x[2]  >= h$y.breaks[x_break]){
      return(0)
    }
    x_point = which(order(c(x[1], h$x.breaks)) == 1)
    y_point = which(order(c(x[2], h$y.breaks)) == 1)
    return(h$counts[x_point - 1, y_point - 1]/sum(h$counts))
  }
  return(f_hist)
}

####1d empirical CDF####
cdf_hist = function(h){
  if(class(h) != 'hist2d') stop('input must be a hist2d object')
  x_counts = rowSums(h$counts)
  y_counts = colSums(h$counts)
  x_cdf = function(x){
    x_point = which(order(c(x, h$x.breaks)) == 1)
    if(x_point == 1) return(0)
    if(x_point > length(h$x.breaks)) return(1)
    
    if(x_point == 2){
      before = 0
    }else{before = sum(x_counts[1:(x_point-2)])/sum(x_counts)}
    
    cur = (x - h$x.breaks[x_point - 1])/(h$x.breaks[x_point] - h$x.breaks[x_point - 1])*x_counts[x_point - 1]/sum(x_counts)
    res  = before + cur
    names(res) = NULL
    return(res)
  }
  
  y_cdf = function(y){
    y_point = which(order(c(y, h$y.breaks)) == 1)
    if(y_point == 1) return(0)
    if(y_point > length(h$y.breaks)) return(1)
    
    if(y_point == 2){
      before = 0
    }else{before = sum(y_counts[1:(y_point-2)])/sum(y_counts)}
    
    cur = (y - h$y.breaks[y_point - 1])/(h$y.breaks[y_point] - h$y.breaks[y_point - 1])*y_counts[y_point - 1]/sum(y_counts)
    res  = before + cur
    names(res) = NULL
    return(res)
  }
  return(c(x_cdf, y_cdf))
}

cdf_hist_x = function(h){
  if(class(h) != 'hist2d') stop('input must be a hist2d object')
  x_counts = rowSums(h$counts)
  x_cdf = function(x){
    x_point = which(order(c(x, h$x.breaks)) == 1)
    if(x_point == 1) return(0)
    if(x_point > length(h$x.breaks)) return(1)
    
    if(x_point == 2){
      before = 0
    }else{before = sum(x_counts[1:(x_point-2)])/sum(x_counts)}
    
    cur = (x - h$x.breaks[x_point - 1])/(h$x.breaks[x_point] - h$x.breaks[x_point - 1])*x_counts[x_point - 1]/sum(x_counts)
    res  = before + cur
    names(res) = NULL
    return(res)
  }
  
  return(x_cdf)
}

cdf_hist_y = function(h){
  if(class(h) != 'hist2d') stop('input must be a hist2d object')
  y_counts = colSums(h$counts)
  
  y_cdf = function(y){
    y_point = which(order(c(y, h$y.breaks)) == 1)
    if(y_point == 1) return(0)
    if(y_point > length(h$y.breaks)) return(1)
    
    if(y_point == 2){
      before = 0
    }else{before = sum(y_counts[1:(y_point-2)])/sum(y_counts)}
    
    cur = (y - h$y.breaks[y_point - 1])/(h$y.breaks[y_point] - h$y.breaks[y_point - 1])*y_counts[y_point - 1]/sum(y_counts)
    res  = before + cur
    names(res) = NULL
    return(res)
  }
  return(y_cdf)
}
####1d empirical pdf####
pdf_hist = function(h){
  if(class(h) != 'hist2d') stop('input must be a hist2d object')
  x_counts = rowSums(h$counts)
  y_counts = colSums(h$counts)
  x_pdf = function(x){
    x_point = which(order(c(x, h$x.breaks)) == 1)
    if(x_point == 1) return(0)
    if(x_point > length(h$x.breaks)) return(0)
    res = x_counts[x_point - 1]/sum(x_counts)
    names(res) = NULL
    return(res)
    
  }
  
  y_pdf = function(y){
    y_point = which(order(c(y, h$y.breaks)) == 1)
    if(y_point == 1) return(0)
    if(y_point > length(h$y.breaks)) return(0)
    res = y_counts[y_point - 1]/sum(y_counts)
    names(res) = NULL
    return(res)
  }
  return(c(x_pdf, y_pdf))
}

pdf_hist_x = function(h){
  if(class(h) != 'hist2d') stop('input must be a hist2d object')
  x_counts = rowSums(h$counts)
  x_pdf = function(x){
    x_point = which(order(c(x, h$x.breaks)) == 1)
    if(x_point == 1) return(0)
    if(x_point > length(h$x.breaks)) return(0)
    res = x_counts[x_point - 1]/sum(x_counts)
    names(res) = NULL
    return(res)
    
  }
  
  return(x_pdf)
}

pdf_hist_y = function(h){
  if(class(h) != 'hist2d') stop('input must be a hist2d object')
  y_counts = colSums(h$counts)
  
  y_pdf = function(y){
    y_point = which(order(c(y, h$y.breaks)) == 1)
    if(y_point == 1) return(0)
    if(y_point > length(h$y.breaks)) return(0)
    res = y_counts[y_point - 1]/sum(y_counts)
    names(res) = NULL
    return(res)
  }
  return( y_pdf)
}

####likelihood####
likelihood_2d_copula = function(h,f_copula, int_level = 15){
  #f_copula is a 2-d function with specified parameter
  #the input h should be a hist2d object
  #int_level is the number of breaks in each dimension while doing the integral
  if(class(h) != 'hist2d') stop('input must be a hist2d object')
  f_hist = h2dtofunction(h) ## convert it to a function
  x_break = length(h$x.breaks)
  y_break = length(h$y.breaks)
  ##following does the integral
  nw <- createNIGrid(dim=2, type="GLe", level=int_level)
  rescale(nw, domain = matrix(c(h$x.breaks[1], h$y.breaks[1],
                                h$x.breaks[x_break], h$y.breaks[y_break]),nrow = 2),
          dec.type = 0)
  f_int = function(x) sqrt(f_hist(x)*f_copula(x))
  result = sum(apply(getNodes(nw), 1, f_int)*getWeights(nw))
  return(result)
}

####return a likelihood function####
likelihoodfunction = function(list_hist, Copula_name){
  if(Copula_name == 'Gaussian'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        log_c = function(x){
          dcopula_Gaussian(c(x_cdf(x[1]), y_cdf(x[2])), rho = theta)
        }
        sum_likelihood = sum_likelihood + log(likelihood_2d_copula(h, log_c))
      }
      return(sum_likelihood)
    }
  }else if(Copula_name == 'Clayton'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        pdfs = pdf_hist((h))
        x_pdf = pdfs[1][[1]]
        y_pdf = pdfs[2][[1]]
        log_c = function(x){
          dcopula_clayton(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
        }
        sum_likelihood = sum_likelihood + log(likelihood_2d_copula(h, log_c))
      }
      return(sum_likelihood)
    }
  }else if(Copula_name == 'Frank'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        pdfs = pdf_hist((h))
        x_pdf = pdfs[1][[1]]
        y_pdf = pdfs[2][[1]]
        log_c = function(x){
          dcopula_frank(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
        }
        sum_likelihood = sum_likelihood + likelihood_2d_copula(h, log_c)
      }
      return(sum_likelihood)
    }
  }else if(Copula_name == 'Gumbel'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        pdfs = pdf_hist((h))
        x_pdf = pdfs[1][[1]]
        y_pdf = pdfs[2][[1]]
        log_c = function(x){
          dcopula_gumbel(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
        }
        sum_likelihood = sum_likelihood + likelihood_2d_copula(h, log_c)
      }
      return(sum_likelihood)
    }
  }else if(Copula_name == 'AMH'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        pdfs = pdf_hist((h))
        x_pdf = pdfs[1][[1]]
        y_pdf = pdfs[2][[1]]
        log_c = function(x){
          dcopula_amh(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)#*x_pdf(x[1])*y_pdf(x[2])
        }
        sum_likelihood = sum_likelihood + likelihood_2d_copula(h, log_c)
      }
      return(sum_likelihood)
    }
  }
  return(likelihood2max)
}

####Copula_MLE####
copula_MLE = function(list_hist, Copula_name, para_range = c(0,1)){
  ##must specify the range
  if(!Copula_name%in%c('Gaussian', 't', 'Clayton', 'Gumbel', 'Frank', 'AMH')){
    stop('Copula_name must be in {Gaussian, t, Clayton, Gumbel, Frank, AMH}')
  }
  if(Copula_name == 'Gaussian'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        log_c = function(x){
          dcopula_Gaussian(c(x_cdf(x[1]), y_cdf(x[2])), rho = theta)
        }
        sum_likelihood = sum_likelihood + log(likelihood_2d_copula(h, log_c))
      }
      return(sum_likelihood)
    }
  }else if(Copula_name == 'Clayton'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        pdfs = pdf_hist((h))
        x_pdf = pdfs[1][[1]]
        y_pdf = pdfs[2][[1]]
        log_c = function(x){
          dcopula_clayton(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
        }
        sum_likelihood = sum_likelihood + log(likelihood_2d_copula(h, log_c))
      }
      return(sum_likelihood)
    }
  }else if(Copula_name == 'Frank'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        pdfs = pdf_hist((h))
        x_pdf = pdfs[1][[1]]
        y_pdf = pdfs[2][[1]]
        log_c = function(x){
          dcopula_frank(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
        }
        sum_likelihood = sum_likelihood + likelihood_2d_copula(h, log_c)
      }
      return(sum_likelihood)
    }
  }else if(Copula_name == 'Gumbel'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        pdfs = pdf_hist((h))
        x_pdf = pdfs[1][[1]]
        y_pdf = pdfs[2][[1]]
        log_c = function(x){
          dcopula_gumbel(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
        }
        sum_likelihood = sum_likelihood + likelihood_2d_copula(h, log_c)
      }
      return(sum_likelihood)
    }
  }else if(Copula_name == 'AMH'){
    likelihood2max = function(theta){
      sum_likelihood = 0
      for(i in 1:n){
        h = list_hist[[i]]
        cdfs = cdf_hist((h))
        x_cdf = cdfs[1][[1]]
        y_cdf = cdfs[2][[1]]
        pdfs = pdf_hist((h))
        x_pdf = pdfs[1][[1]]
        y_pdf = pdfs[2][[1]]
        log_c = function(x){
          dcopula_amh(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)#*x_pdf(x[1])*y_pdf(x[2])
        }
        sum_likelihood = sum_likelihood + likelihood_2d_copula(h, log_c)
      }
      return(sum_likelihood)
    }
  }
  
  
  mle_result = optimize(f = likelihood2max, interval = para_range,maximum = T)
  result = NULL
  result$mle = mle_result$maximum
  result$max_likelihood = mle_result$objective
  return(result)
}

####histogram list to estimate one pdf function####
hlist2dtofunction = function(h_list){
  f_list = function(x){
    result = 0
    for(h in h_list){
      f = h2dtofunction(h)
      result = result + f(x)
    }
    return(result/length(h_list))
  }
  return(f_list)
}

####combined marginal cdf and pdf for list####
cdf_list_hist_x = function(h_list){
  cdfx = function(x){
    result = 0
    for(h in h_list){
      f = cdf_hist_x(h)
      result = result + f(x)
    }
    return(result/length(h_list))
  }
  return(cdfx)
}

cdf_list_hist_y = function(h_list){
  cdfy = function(x){
    result = 0
    for(h in h_list){
      f = cdf_hist_y(h)
      result = result + f(x)
    }
    return(result/length(h_list))
  }
  return(cdfy)
}

pdf_list_hist_x = function(h_list){
  pdfx = function(x){
    result = 0
    for(h in h_list){
      f = pdf_hist_x(h)
      result = result + f(x)
    }
    return(result/length(h_list))
  }
  return(pdfx)
}

pdf_list_hist_y = function(h_list){
  pdfy = function(x){
    result = 0
    for(h in h_list){
      f = pdf_hist_y(h)
      result = result + f(x)
    }
    return(result/length(h_list))
  }
  return(pdfy)
}

#### likelihood for list####
likelihood_2d_multiple = function(h_list,f_copula, int_level = 20){
  #f_copula is a 2-d function with specified parameter
  #int_level is the number of breaks in each dimension while doing the integral
  f_hist = hlist2dtofunction(h_list)
  x_min = Inf
  x_max = -Inf
  y_min = Inf
  y_max = -Inf
  for(h in h_list){  
    x_break = length(h$x.breaks)
    y_break = length(h$y.breaks)
    x_min = min(x_min, h$x.breaks[1])
    y_min = min(y_min, h$y.breaks[1])
    x_max = max(x_max, h$x.breaks[x_break])
    y_max = max(y_max, h$y.breaks[y_break])
  }
  ##following does the integral
  nw <- createNIGrid(dim=2, type="GLe", level=int_level)
  rescale(nw, domain = matrix(c(x_min, y_min,x_max, y_max),nrow = 2),
          dec.type = 0)
  f_int = function(x) sqrt(f_hist(x)*f_copula(x))
  result = sum(apply(getNodes(nw), 1, f_int)*getWeights(nw))
  return(result)
}

####likelihood for list####
likelihoodfunctionlist = function(list_hist, Copula_name){
  x_cdf = cdf_list_hist_x(list_hist)
  y_cdf = cdf_list_hist_y(list_hist)
  x_pdf = pdf_list_hist_x(list_hist)
  y_pdf = pdf_list_hist_y(list_hist)
  if(Copula_name == 'Gaussian'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_Gaussian(c(x_cdf(x[1]), y_cdf(x[2])), rho = theta)
      }
      return(log(likelihood_2d_multiple(h, log_c)))
    }
  }else if(Copula_name == 'Clayton'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_clayton(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
      }
      return(log(likelihood_2d_multiple(h, log_c)))
    }
  }else if(Copula_name == 'Frank'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_frank(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
      }
      return(likelihood_2d_multiple(h, log_c))
    }
  }else if(Copula_name == 'Gumbel'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_gumbel(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
      }
      return(likelihood_2d_multiple(h, log_c))
    }
  }else if(Copula_name == 'AMH'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_amh(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
      }
      return(likelihood_2d_multiple(h, log_c))
    }
  }
  return(likelihood2max)
}


####Copula_MLE list####
copula_MLE_list = function(list_hist, Copula_name, para_range = c(0,1)){
  ##must specify the range
  if(!Copula_name%in%c('Gaussian', 't', 'Clayton', 'Gumbel', 'Frank', 'AMH')){
    stop('Copula_name must be in {Gaussian, t, Clayton, Gumbel, Frank, AMH}')
  }
  x_cdf = cdf_list_hist_x(list_hist)
  y_cdf = cdf_list_hist_y(list_hist)
  x_pdf = pdf_list_hist_x(list_hist)
  y_pdf = pdf_list_hist_y(list_hist)
  if(Copula_name == 'Gaussian'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_Gaussian(c(x_cdf(x[1]), y_cdf(x[2])), rho = theta)
      }
      return(log(likelihood_2d_multiple(list_hist, log_c)))
    }
  }else if(Copula_name == 'Clayton'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_clayton(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
      }
      return(log(likelihood_2d_multiple(list_hist, log_c)))
    }
  }else if(Copula_name == 'Frank'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_frank(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
      }
      return(likelihood_2d_multiple(list_hist, log_c))
    }
  }else if(Copula_name == 'Gumbel'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_gumbel(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
      }
      return(likelihood_2d_multiple(list_hist, log_c))
    }
  }else if(Copula_name == 'AMH'){
    likelihood2max = function(theta){
      log_c = function(x){
        dcopula_amh(c(x_cdf(x[1]), y_cdf(x[2])), theta = theta)*x_pdf(x[1])*y_pdf(x[2])
      }
      return(likelihood_2d_multiple(list_hist, log_c))
    }
  }
  
  
  mle_result = optimize(f = likelihood2max, interval = para_range,maximum = T)
  result = NULL
  result$mle = mle_result$maximum
  result$max_likelihood = mle_result$objective
  return(result)
}

####Variance####
variance_copula = function(h_list){
  result = NULL
  #set up the integral
  x_min = Inf
  x_max = -Inf
  y_min = Inf
  y_max = -Inf
  for(h in h_list){  
    x_break = length(h$x.breaks)
    y_break = length(h$y.breaks)
    x_min = min(x_min, h$x.breaks[1])
    y_min = min(y_min, h$y.breaks[1])
    x_max = max(x_max, h$x.breaks[x_break])
    y_max = max(y_max, h$y.breaks[y_break])
  }
  x_pdf = pdf_list_hist_x(list_hist)
  y_pdf = pdf_list_hist_y(list_hist)
  
  #x
  nw <- createNIGrid(dim=1, type="GLe", level=40)
  rescale(nw, domain = c(x_min, x_max),dec.type = 0)
  f_int = function(x) x_pdf(x)*x
  x_mean = sum(sapply(getNodes(nw), f_int)*getWeights(nw))/sum(sapply(getNodes(nw), x_pdf)*getWeights(nw))
  
  f_int = function(x) x_pdf(x)*x^2
  x2_mean = sum(sapply(getNodes(nw), f_int)*getWeights(nw))/sum(sapply(getNodes(nw), x_pdf)*getWeights(nw))
  x_var = x2_mean - x_mean^2
  
  #y
  rescale(nw, domain = c(y_min, y_max),dec.type = 0)
  f_int = function(x) y_pdf(x)*x
  y_mean = sum(sapply(getNodes(nw), f_int)*getWeights(nw))/sum(sapply(getNodes(nw), y_pdf)*getWeights(nw))
  
  f_int = function(x) y_pdf(x)*x^2
  y2_mean = sum(sapply(getNodes(nw), f_int)*getWeights(nw))/sum(sapply(getNodes(nw), y_pdf)*getWeights(nw))
  y_var = y2_mean - y_mean^2
  
  result$x_mean = x_mean
  result$y_mean = y_mean
  result$x_var = x_var
  result$y_var = y_var
  
  return(result)
}

####Covariance####
covariance_copula = function(h_list, Copula_name = 'Gaussian', theta, int_level = 10){
  ##this is to calculate the covariance with the estiamted/given theta
  #set up the integral
  x_min = Inf
  x_max = -Inf
  y_min = Inf
  y_max = -Inf
  for(h in h_list){  
    x_break = length(h$x.breaks)
    y_break = length(h$y.breaks)
    x_min = min(x_min, h$x.breaks[1])
    y_min = min(y_min, h$y.breaks[1])
    x_max = max(x_max, h$x.breaks[x_break])
    y_max = max(y_max, h$y.breaks[y_break])
  }
  nw <- createNIGrid(dim=2, type="GLe", level=int_level)
  rescale(nw, domain = matrix(c(x_min, y_min,x_max, y_max),nrow = 2),
          dec.type = 0)
  x_cdf = cdf_list_hist_x(list_hist)
  y_cdf = cdf_list_hist_y(list_hist)
  
  #decide the copula
  if(Copula_name == 'Gaussian'){
    copula2use = copula_Gaussian
  }else if(Copula_name == 'Clayton'){
    copula2use = copula_clayton
  }else if(Copula_name == 'Frank'){
    copula2use = copula_frank
  }else if(Copula_name == 'Gumbel'){
    copula2use = copula_gumbel
  }else if(Copula_name == 'AMH'){
    copula2use = copula_amh
  }
  #the function to integral
  f_int = function(x){
    copula2use(c(x_cdf(x[1]),y_cdf(x[2])), theta) - x_cdf(x[1])*y_cdf(x[2])
  }
  result = NULL
  result$covariance = sum(apply(getNodes(nw), 1, f_int)*getWeights(nw))
  
  vars = variance_copula(h_list)
  result$correlation = result$covariance/sqrt(vars$x_var*vars$y_var)
  return(result)
}
