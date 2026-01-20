surval.x <- function(t, x1, x2, x3){        
  surval.x <- exp(-Lambda0(t)-beta1.true*x1*t-beta2.true*x2*t-beta3.true*x3*t)
  return(surval.x)
} 

fail <- function(ran.prop, x1, x2, x3){
  tt<- uniroot(function(t){1- surval.x(t,x1,x2,x3)-ran.prop},c(0,100))$root
  return(tt)
}

get_data <- function(n, l1=0.02, u1=1.28){
  x1 <- runif(n, min = 0, max = 1)
  x2 <- runif(n, min = 0, max = 1)
  x3 <- rbinom(n, size = 1, prob = 0.5)
  
  xx<-cbind(x1,x2,x3)
  ran.prop<- runif(n, min = 0, max = 1)
  
  fail.time <- obse.time <- mapply(fail, ran.prop, x1, x2, x3)
  cens.time <- runif(n, min = l1, max = u1)
  delta <- as.numeric(I(fail.time <= cens.time))
  #cens.rate <- 1 - mean(delta)
  obse.time[!delta] <- cens.time[!delta]
  
  data <- data.frame(label_time = obse.time, delta= delta, x1=x1, x2=x2, x3=x3)
  
  return(data)
}

# Additive risk model
# Local estimator
AR.fit <- function(yobs, delta, X){
  N <- length(yobs)
  p <- ncol(X)
  X <- X %>% mutate_if(is.factor, ~ ifelse(. == levels(.)[1], 0, 1))
  X <- as.matrix(X)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  
  X.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){
    Y <- ifelse((yobs >= y.sort[j])==TRUE, 1, 0)
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- diff( c(0,y.sort) )
  
  A0 <- B0 <- array(0, dim=c(p, p))
  d0 <- rep(0,p)
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    Rdi <- Ri * sqrt(di)
    I2i <- t(Rdi) %*% Rdi # I2i <- t(Ri) %*% di %*% X[i,]
    I1i <- ( X[i,] - X.bar.sort[y.rank[i],] ) * delta[i]
    A0 <- A0 + I2i
    d0 <- d0 + I1i
    B0 <- B0 + delta[i]*( Ri[Ki,]%*%t(Ri[Ki,]) )
  }
  A <- A0/N; B <- B0/N; d <- d0/N
  # calculate the estimate of beta and SE
  Est    <- solve(A,d)
  Sigma  <- solve(A) %*% B %*% solve(A) # asymptotic var-cov matrix
  SE     <- sqrt( diag(Sigma)/N )
  zvalue <- Est/SE
  pvalue <- 2*(1-pnorm(abs(zvalue)))
  coef <- data.frame(Est=Est, SE=SE, zvalue=zvalue, pvalue=pvalue,
                     row.names=colnames(X))
  # output
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=X),
    coef=coef
  )
  return(out)
}

V_k <- function(yobs, delta, X){
  N <- length(yobs)
  p <- ncol(X)
  X <- X %>% mutate_if(is.factor, ~ ifelse(. == levels(.)[1], 0, 1))
  X <- as.matrix(X)
  y_sort <- sort( yobs )
  y_rank <- rank( yobs, ties.method = 'min')
  X_bar_sort <- array(0, dim=c(N, p))
  for( j in 1:N){
    Y <- ifelse((yobs >= y_sort[j])==TRUE, 1, 0)
    X_bar_sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  
  y_sort_diff <- diff( c(0,y_sort) )
  # calculate A, B and d
  A0 <- B0 <- array(0, dim=c(p, p))
  d0 <- rep(0,p)
  for( i in 1:N ){ # i
    Ki <- y_rank[i]
    Xi_aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi_aug - X_bar_sort[1:Ki,]
    di <- y_sort_diff[1:Ki]
    Rdi <- Ri * sqrt(di)
    I2i <- t(Rdi) %*% Rdi 
    I1i <- ( X[i,] - X_bar_sort[y_rank[i],] ) * delta[i]
    A0 <- A0 + I2i
    d0 <- d0 + I1i
    B0 <- B0 + delta[i]*( Ri[Ki,]%*%t(Ri[Ki,]) )
  }
  out <- list(N = N, A_k = A0, coef=names(d0), d_k = as.vector(d0), B_k = B0 )
  return(out)
} 


write_site_V <- function(sitename, infor) {
  mylist = list(site = sitename, N=infor$N, 
                A_k = infor$A_k,  d_k = infor$d_k, B_k = infor$B_k)
  all_site_V <- rjson::fromJSON(file = "site_infor.json")
  K = length(all_site_V)
  found = FALSE
  for (i in 1:K) {
    if (all_site_V[[i]]$site == sitename) {
      all_site_V[[i]] = mylist
      found = TRUE
    }
  }
  if (!found) {
    all_site_V[[K+1]] = mylist
  }
  json_dat = rjson::toJSON(all_site_V, 2)
  #print(json_data)
  write(json_dat, "site_infor.json")
}

read_site_V <- function(sitename) {
  all_infor<- rjson::fromJSON(file = "site_infor.json")
  K = length(all_infor)
  for (i in 1:K) {
    if (all_infor[[i]]$site == sitename) {
      return(all_infor[[i]])
    }
  }
  return()
}

get_Fedbeta <- function(site_list, var_p){

  NN <- 0
  d_all <- rep(0, var_p)
  A_all <- B_all <- matrix(0, var_p, var_p)
  for(j in 1:length(site_list)){
    V_site <- read_site_V(sitename = site_list[j])
    A_all <- A_all + matrix(V_site$A_k, nrow = var_p) 
    d_all <- d_all + V_site$d_k
    B_all <- B_all + matrix(V_site$B_k, nrow = var_p)
    NN <- NN + V_site$N
  }
  A_all <- A_all/NN; d_all <- d_all/NN; B_all <- B_all/NN
  Est    <- solve(A_all, d_all)
  Sigma  <- solve(A_all) %*% B_all %*% solve(A_all) # asymptotic var-cov matrix
  SE     <- sqrt( diag(Sigma)/NN )
  zvalue <- Est/SE
  pvalue <- 2*(1-pnorm(abs(zvalue)))
  coef <- data.frame(Est=Est, SE=SE, zvalue=zvalue, pvalue=pvalue)
  return(coef)
  
}

# Calculate Lambda0(t) at time points y_obs (time=NULL)
AR.Ht <- function(data, bet, time=NULL){
  
  yobs <- data$label_time
  delta <- data$label_status
  X <- data[, !(names(data) %in% c("label_time","delta"))]
  N <- length(yobs); p <- ncol(X)
  
  num_Risk  <- sapply( yobs, function(yi){sum( yobs >= yi )} )
  y_sort <- sort( yobs )
  y_sort_diff <- diff( c(0, y_sort) )
  X_bar_sort <- array(0, dim = c(N, p))
  cumhaz0_yobs <- rep(NA, N)
  for( j in 1:N){
    # calculate X.bar.sort till this time point yobsj
    Y <- ifelse((yobs >= y_sort[j])==TRUE, 1, 0)
    X_bar_sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
    # calculate the C(yobsj) at this time point
    Ct <- t(X_bar_sort[1:j,,drop=FALSE]) %*% y_sort_diff[1:j]
    # calculate the first part
    Ft <- sum( (yobs <= y_sort[j]) * delta / num_Risk  )
    # Final Result
    cumhaz0_yobs[j] <- Ft - as.numeric( t(bet) %*% Ct )
  }
  
  ### to make it momtone
  cumhaz0_yobs_mon <- sapply(1:length(cumhaz0_yobs),
                             function(icum){max(c(0,cumhaz0_yobs[1:icum]))}
                             )
  
  ### calculate Lambda0(t) at time points tm ###
  if(is.null(time)){
    time <- y_sort
    cumhaz0 <- cumhaz0_yobs_mon
  }else{
    cumhaz0 <- c(0,cumhaz0_yobs_mon)[
      sapply(time, function(tmi){sum( c(0, y_sort) <= tmi )}) ]
  }
  
  ### Output ###
  out <- data.frame(time=time, cumhaz0=cumhaz0)
  return(out)
}

# S(Y|X) if time=NULL, x=NULL
AR.Stx <- function(data, bet, time=NULL, x=NULL){
  
  ### some preparation ###
  yobs <- data$label_time
  delta <- data$label_status
  X <- data[, !(names(data) %in% c("label_time","label_status"))]
  
  if(is.null(time)){time <- yobs}
  if(is.null(x)){x <- X}
  if(is.null(dim(x))){x <- t(x)}
  l <- nrow(x) 
  m <- length(time)
  
  ### Obtain the Baseline Cumhaz by using AR.Ht ###
  cumhaz0 <- AR.Ht(data=data, bet = bet, time = time)[,2]
  ### Give estimates for S(t|x) - l*m ###
  Stx <- matrix(nrow=l, ncol=m, 
                dimnames=list(paste("x", 1:l, sep=""),
                              paste("time", 1:m, sep="")))
  for(im in 1:m){
    for(il in 1:l){
      Stx[il, im] <- exp(-cumhaz0[im]-sum(x[il,]*bet)*time[im])
    }
  }
  
  ### Output ###
  out <- Stx #list(time=time, x=x, Stx=Stx)
  return(out)
}


#======For unstratified model ===================================

Xk_bar <- function(yobs, X, y_sort){
  n <- length(yobs)
  N <- length(y_sort)
  p <- ncol(X)

  Xk_bar_fenmu <-NULL
  Xk_bar_fenzi <- array(0, dim=c(N, p))
  for( j in 1:N){
    Yj <- (yobs >= y_sort[j]) #N*1
    Xk_bar_fenzi[j,] <- apply( X * Yj, 2, sum ) 
    Xk_bar_fenmu[j] <- sum(Yj)
  }
  return(list(Xk_bar_fenzi=Xk_bar_fenzi, Xk_bar_fenmu=Xk_bar_fenmu))
}


#N*p for site k
V_k_un <- function(yobs, delta, X, y_sort, X_bar){
  N <- length(y_sort)
  n <- length(yobs)
  p <- ncol(X)
  X <- as.matrix(X)
  positions <- sapply(yobs, function(x) which(y_sort == x))
  
  y_sort_diff <- diff( c(0,y_sort) )
  # calculate A, B and d
  A0 <- B0 <- array(0, dim=c(p, p))
  d0 <- rep(0,p)
  for( i in 1:n ){ # i
    
    Ki <- max(unlist(positions[i]))
    Xi_aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi_aug - X_bar[1:Ki,]
    di <- y_sort_diff[1:Ki]
    Rdi <- Ri * sqrt(di) #N*p
    I2i <- t(Rdi) %*% Rdi # I2i <- t(Ri) %*% di %*% X[i,]#fenmu   sum at risk p*p
    I1i <- ( X[i,] - X_bar[Ki,] ) * delta[i] #fenzi #p*1
    A0 <- A0 + I2i #p*p
    d0 <- d0 + I1i  #p*1
    B0 <- B0 + delta[i]*( Ri[Ki,]%*%t(Ri[Ki,]) )
  }
  # output
  out <- list(N = N, A_k = A0, coef=names(d0), d_k = as.vector(d0), B_k = B0 )
  return(out)
} 

