# Minimum Volatility method for optimal window  size
# simulation for master thesis
# author: Xiaotong Meng
# E-mail: xiaotong.meng@rwth-aachen.de


set.seed(2022)
library("zoo")
######################functions#######################

rollapply_naive <- function(x, m, mM){
  T = length(x) 
  TM <- T-mM+1 
  S_jm <- sapply(1: TM, function(y){sum(x[y: (y+m-1)])})
  #print(paste("length(S_jm)" , length(S_jm)))
  return(S_jm)
}


get_gamma_naive <- function(m,x,mM){
  T <- length(x) 
  S_jm  <-rollapply_naive(x, m, mM)
  Sn <- sum(x)
  tmp <- ( S_jm - m/T *Sn)^2 / (m*(T-m+1))
  tmp1 <-  cumsum(tmp) 
  return(tmp1)
  
}

func2 <- function(x, k){
  temp =   rollapply(x, k, sum)
  return(temp/k)
  
}

func3 <- function(i,k, gm1,gm2){
  j = i + k -1    
  temp <- gm1[i:j,] - gm2[i,]
  temp <- temp^2
  temp <- apply(temp, 2, sum)
  temp <- sqrt(temp)
  return(temp)
}

optimize_window <- function(x, m1, mM, mm, k){
  gc()  
  T <- length(x)
  m_all <- seq(m1, mM, by = mm)
  mM_real <- tail(m_all, n=1)
  #cat("# window candidate", length(m_all))
  gamma_all <- t(sapply(m_all, get_gamma_naive, x = x, mM = mM_real))
  gamma_bar <- apply(gamma_all, 2, func2, k = k)
  jn <- nrow(gamma_bar)
  j<- seq(1,jn, by = 1)
  se <- t(sapply(j, func3, k = k, gm1 =gamma_all,gm2 = gamma_bar))
  max_se <- apply(se,1, max)
  min_max = min(max_se)
  m_MV_index <- which(max_se == min_max)
  m_MV <- m_all[m_MV_index + (k-1)/2]
  m <-list("window" = m_MV, "min_max_Se" = min_max)
  return(m)
}

########################initialization########################
# m <- optimize_window_naive(x, 4, 2340, 1, 7)
# m