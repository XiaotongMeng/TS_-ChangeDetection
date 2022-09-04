# classic CUSUM test with lagged window estimator
# simulation for master thesis
# author: Xiaotong Meng
# E-mail: xiaotong.meng@rwth-aachen.de


set.seed(2022)
library("zoo")
library('e1071')


######################functions#######################

sim_tvAR1 <-function(a, TT, mu){
  ww <- rnorm(TT)
  tvAR1 <- rep(0, TT)
  tvAR1[1] <- ww[1]
  for(t in 2:TT) { 
    tvAR1[t] <- a[t-1]*tvAR1[t-1]  + ww[t] 
  }
  tvAR1 <- tvAR1+ mu
  return(tvAR1)
}

# Optimize window with AR1 method
get_window <- function(x){
  TT = length(x)
  ar = acf(x, plot=FALSE)
  a = ar$acf[2]
  m = floor( 1.1447*( 4*a^2 *TT / (1- a^2) )^(1/3) )
  if( m==0) {m =1}
  return(m)
} 

#test statistics with lagged window estimator
Tn_LW <- function(i, a, TT, mu){
  j = 1: TT
  j = j/TT
  x <- sim_tvAR1(a,TT,mu)
  #m <- optimize_window (x, 1, 0.1 *TT, 1, 7)
  #m <- m$window
  m <- get_window(x)
  #cat("Window size:",m )
  mvsum = rollapply(x, m, sum)
  LW = (mvsum/m - sum(x)/TT)^2
  LW <- m/(TT-m+1)* sum(LW) 
  #cat("lagged window estimator", LW )
  Tn <-  max(abs(cumsum(x) - j*sum(x))) / (sqrt(TT* LW) )
  return(Tn) 
}


#maximum of Brownian bridge
get_mbb <- function(i,f){
  tn = rbridge(end = 1, frequency = f)
  tn = max(abs(tn))
  return(tn)
}



get_p <- function(Tn, tn, Nmbb){
  p = max(which(tn <= Tn))
  #print(p)
  p = 1- p/Nmbb
  return(p)
}



get_phist <- function(Nmbb, a, TT, wB, mu){
  j <- rep(1, Nmbb)
  tn_all <- sapply(j, get_mbb, f=1000)
  tn_s <- sort(tn_all)
  i <- rep(1, wB)
  Tn_all <- sapply(i, Tn_LW, a = a, TT=TT, mu=mu)
  #print(Tn_all[1:3])
  p_all <- sapply(Tn_all, get_p, tn = tn_s, Nmbb= Nmbb)
  hist(p_all, main=NULL, xlab =" ")
  hist(p_all, main=NULL, xlim =c(0,1),freq = FALSE, xlab = "")
  return(p_all)
}


#type I Error
get_t1 <- function(p_all, nl){
  t1 <- sum(p_all < nl)
  t1 <- t1/ length(p_all)*100
  return(t1)
}



########################initialization########################

TT <- 200
wB <- 2000
Nmbb<- 1000

a0 <- rep(0.6, TT)
a1 <- c(rep(0.6, TT/2), rep(-0.6, TT/2))
t <-seq(0,1, length.out =TT) * 2 * pi
a2 <- 0.75*cos(t)
t <-seq(0,1, length.out =TT/2) * 2 * pi
a3_1 <- 0.75*cos(t)
a3_2 <- seq(0.5, 1, length.out = TT/2)
a3 <- c(a3_1, 0.6-a3_2)

mu <- rep(0, TT)
#mu <- c(rep(1, TT/4),rep(-1, TT/4), rep(1, TT/4),rep(-1, TT/4))



########################test###################################

p_all0 <- get_phist(Nmbb, a0,TT, wB,mu)
a <- 0.05
t0 <-get_t1(p_all0, a)