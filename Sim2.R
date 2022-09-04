# CUSUM test with robust bootstrap method
# simulation for master thesis
# author: Xiaotong Meng
# E-mail: xiaotong.meng@rwth-aachen.de

set.seed(2022)
library("zoo")
#library('e1071')


######################functions#######################

sim_tvAR1 <-function(a, TT, mu){
  ww <- rnorm(TT)
  tvAR1 <- rep(0, TT)
  tvAR1[1] <- ww[1]
  for(t in 2:TT) { 
    tvAR1[t] <- a[t]*tvAR1[t-1]  + ww[t] 
  }
  tvAR1 <- tvAR1+ mu
  return(tvAR1)
}

robust_Tn <- function(x){
  gc()
  TT <- length(x)
  t = seq(0, 1, length.out = TT)
  Tn =  max(abs(cumsum(x) - t*sum(x))/ sqrt(TT))
  #plot(Tn, pch =20,col="blue", xlab="Second")
  #which(Tn < 0.3
  #print( paste("Tn = ", Tn))
  return(Tn)
}

roll_sum <- function(x, m){
  TT = length(x) 
  Tm <- TT-m+1 
  j <- seq(1, Tm, by =1)
  S_jm <- sapply(j, function(y){sum( x[y:(y+m-1)])} )
  #print(paste("length(S_jm)" , length(S_jm)))
  return(S_jm)
}

robust_tn <- function(i,m,x,S_jm){
  TT <- length(x)
  gN <-TT-m+1
  gj <- rnorm(gN, 0 ,1)
  temp <- (S_jm - m*sum(x)/TT )* gj/ sqrt( m*(TT-m+1))
  phi_im <- cumsum(temp)
  j <- seq(0,1,length.out =gN)
  tn <- phi_im - j*sum(temp)
  tn1 <- tn[-m]
  Mr <- max(abs(tn1))
  return(Mr)
}    


robust_get_p = function(i,B, a,TT, mu){
  gc()
  Mr <- rep(1, B)
  x <- sim_tvAR1(a,TT,mu)
  #m <- get_window(x)
  m <- optimize_window(x, 1, 0.1*TT, 1, 7)
  m <- m$window
  #print(m)
  S_jm <- roll_sum(x =x, m =m)
  Mr <- sapply(Mr, robust_tn, m= m, x = x, S_jm =S_jm)
  #hist(Mr)
  Mr <- sort(Mr)
  Tn <- robust_Tn(x)
  p <- max(which(Mr <= Tn))
  p <- 1- p/B
  return(p)  
}


robust_get_phist <- function(B,Np, a,TT,mu){
  i <- rep(1,Np) 
  p_all <- sapply(i, robust_get_p, B = B, a=a, TT=TT, mu = mu) 
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
B <- 2000
Np<- 1000

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

r_p_all0 <- robust_get_phist(B,Np, a0, TT,mu)
a <- 0.05
t0 <-get_t1(p_all0, a)