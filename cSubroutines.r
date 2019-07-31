#Fraga Alve's estimate of the second order parameter
rho.fraga.alves <-function(x)
{
  n<-length(x)
  x<-sort(x,decreasing = TRUE)
  U<-numeric(n)
  for (i in 1:(n-1))   U[i]=i*log(x[i]/x[i+1])
  
  k1=min((n-1),round(n^0.995))
  Uk=U[1:k1]
  
  M<-numeric(3)
  T<-numeric(3)
  
  xk <- x[1:k1]/x[k1+1]
  M[1]<-sum(log(xk))/k1
  M[2]<-sum(log(xk)^2)/k1
  M[3]<-sum(log(xk)^3)/k1
  T[1]=(log(M[1])-1/2*log(M[2]/2))/(1/2*log(M[2]/2)-1/3*log(M[3]/6))
  T[2]=((M[1])^(1/2)-(M[2]/2)^(1/4))/((M[2]/2)^(1/4)-(M[3]/6)^(1/6))
  T[3]=((M[1])^(1)-(M[2]/2)^(1/2))/((M[2]/2)^(1/2)-(M[3]/6)^(1/3))
  r.tau.tilde=-abs((3*(T[1]-1))/(T[1]-3))
  
  list(rho=r.tau.tilde)
}

# Proportion of non-censored observations
ip.hat.fun <- function(delta.n){
  
  n <- length(delta.n)
  K <- 1:(n-1)
  
  p.hat <- cumsum(delta.n[n-K+1])/K
  #   for (k in 1:(n-1)) {
  #     p.hat[k] <- mean(delta.n[n:(n-k+1)])
  #   }
  return(p.hat)
}

# E_{Z,k}(s)
iEs <- function(s, Z.n){
  
  n <- length(Z.n)
  output <- numeric(n-1)
  
  for (k in 1:(n-1)) {
    v <- (Z.n[n-(1:k)+1]/Z.n[n-k])^(s[k])
    output[k] <- mean(v)
  }
  return(output)
}


# E^(c)_{Z,k}(s)
icEs <- function(s, Z.n, delta.n){
  
  n <- length(Z.n)
  output <- numeric(n-1)
  
  for (k in 1:(n-1)) {
    v <- delta.n[n-(1:k)+1] * (Z.n[n-(1:k)+1]/Z.n[n-k])^(s[k])
    output[k] <- mean(v)
  }
  return(output)
}

# Censored EPD estimator for gamma_1 and kappa_1
#
# rho is a parameter for the rho_estimator of Fraga Alves et al. (2003)
# when strictly positive or (a) choice(s) for rho if negative
cpEPD <- function(data, censored, rho = -1,H_kn, w=1,pen=FALSE,beta = NULL, logk = FALSE, plot = FALSE, add = FALSE, 
                  main = "EPD estimates of EVI", ...){
  
  n <- length(data)
  k <- n-1 
  K<-1:k
  s <- sort(data, index.return = TRUE)
  Z.n <- s$x
  # delta is 1 if a data point is not censored, we also sort it with X
  delta.n <- !(censored[s$ix])
  
  
  phat <- ip.hat.fun(delta.n)
  HillZ <- Hill(Z.n)$gamma
  
  nrho <- length(rho)
  
  if (is.null(beta)) {
    # determine beta using rho
    
    beta <- matrix(0,n-1,nrho)
    
    if (all(rho>0) & nrho==1) {
      rho <- .rhoEst(data,alpha=1,tau=rho)$rho
      
      # Estimates for rho of Fraga Alves et al. (2003) used 
      # and hence a different value of beta for each k
      beta[,1] <- -rho/HillZ
      
      for(j in 1:length(rho)) {
        beta[,j] <- -rho[j]/HillZ
      }
      
    } else if(all(rho<0)) {
      
      # rho is provided => beta is constant over k
      # (but differs with rho)
      for(j in 1:nrho) {
        beta[,j] <- -rho[j]/HillZ
      }
      
      
    } else {
      stop("rho should be a single positive number or a vector (of length >=1) of negative numbers.")
    }
    
  } else {
    nrho <- length(beta)
    
    # use provided value(s) for beta
    if (length(beta)==1) {
      beta <- matrix(beta,n-1,nrho)
      #     } else if(length(beta)==n-1){
      #       beta <- as.matrix(rep(beta,nrho),ncol=length(rho),byrow=FALSE)
      #     } else {
      #       stop(paste0("beta should have length 1 or n-1 = ",n-1,"."))
      #     }
    } else {
      beta <- matrix(rep(beta,n-1),ncol=length(beta),byrow=TRUE)
    }
  }
  
  gamma1 <- matrix(0,n-1,nrho)
  kappa1 <- matrix(0,n-1,nrho)
  Delta <- matrix(0,n-1,nrho)
  sigma.k<-matrix(0,n-1,nrho)
  
  penalty<-rep(0,n-1)
  
  for(j in 1:nrho) {
    sigma.k[K,j]<-(K/n)^(2)#*rho[j]
    if(pen==TRUE) penalty<-(w*H_kn[K])/(K*sigma.k[K,j])
    D <- - (beta[K,j]^4 * HillZ[K]^3) / ( (1+HillZ[K]*beta[K,j])^2 * (1+2*HillZ[K]*beta[K,j]))-penalty
    
    Es <- iEs(-beta[,j],Z.n)
    cEs <- icEs(-beta[,j],Z.n,delta.n)
    
    # Estimates for kappa1
    kappa1[,j] <- (1 - Es[K] - beta[K,j] * (H_kn[K]) * cEs[K]) / D[K]
    kappa1[,j] <- pmax(kappa1[,j], pmax(-1,-1/beta[,j])+0.001)
    
    # Estimates for Delta
    Delta[,j] <- (kappa1[,j]*(1-Es))/phat 
    
    # Has to be strictly negative otherwise the bias is increased,
    # hence we put it to 0 such that gamma1 = cHill then
    Delta[,j] <- pmin(0, Delta[,j])
    
    # gamma1 = cHill + Delta
    gamma1[,j] <- H_kn - Delta[,j]*(rho[j]/(1-rho[j]))
    gamma1[gamma1[,j]<=0,j] <- 0.001
    
  }
  # Transform to vectors if rho is a single value
  if(nrho==1) {
    gamma1 <- as.vector(gamma1)
    kappa1 <- as.vector(kappa1)
    beta <- as.vector(beta)
    Delta <- as.vector(Delta)
  } else if(plot | add) {
    # Add lines
    for(j in 2:nrho) {
      lines(K,gamma1[,j],lty=j)
    }
  }
  return(list(k=K, gamma1=gamma1, kappa1=kappa1, beta=beta, Delta=Delta))
}

BRW <- function(data, censored, rho = -1,L_kn, w=1,pen=T,beta = NULL, logk = FALSE, plot = FALSE, add = FALSE, 
                   main = "EPD estimates of EVI", ...){
  #data<-Z;censored<-1-del1;rho = rhos;L_kn=Worms2; w=1;pen=F;beta = NULL
  n <- length(data)
  k <- n-1 
  K<-1:k
  s <- sort(data, index.return = TRUE)
  Z.n <- s$x
  # delta is 1 if a data point is not censored, we also sort it with X
  delta.n <- !(censored[s$ix])
  rho <- matrix(rep(rho,n-1),ncol=length(rho),byrow=TRUE)
  nrho<-ncol(rho)
  gamma1 <- matrix(0,n-1,nrho)
  delta.c <- matrix(0,n-1,nrho)
  kappa1<- matrix(0,n-1,nrho)

  
  sigma.k<-matrix(0,n-1,nrho)
  penalty<-rep(0,n-1)
  
  for(j in 1:nrho) {
    sigma.k[K,j]<-(K/n)^(2)#*rho[K,j]
    if(pen==TRUE) penalty<-(w*L_kn[K])/(K*sigma.k[K,j])
  
  for (k in K) {
      prods <- numeric(k)
      
      # Compute products for i=1,...,k
      for (i in 1:k) {
        ind <- which(Z.n<Z.n[n-i+1])
        prods[i] <- prod((1-delta.n[ind]/(n-ind+1)))
      }
      delta.c[k,j] <-1+sum(prods/prods[k] * (Z.n[n-(1:k)+1]^(rho[k,j]/L_kn[k])-Z.n[n-(1:k)]^(rho[k,j]/L_kn[k]))/(Z.n[n-k]^(rho[k,j]/L_kn[k]))) 
    }
    kappa1[K,j] <- (rho[K,j]^4/(L_kn[K] * (1-2*rho[K,j]) * (1-rho[K,j])^2)+penalty)^(-1)* (delta.c[K,j] - 1 / (1-rho[K,j]))
    #kappa1[,j] <- pmin(kappa1[K,j], 0)
    # gamma1 = cHill + Delta
    gamma1[,j] <- L_kn[K] -  kappa1[K,j]* rho[K,j] #/ (1 - rho[K,j])
  }
  # Transform to vectors if rho is a single value
  if(nrho==1) {
    gamma1 <- as.vector(gamma1)
    kappa1 <- as.vector(kappa1)
    beta <- as.vector(beta)
    delta.c <- as.vector(delta.c)
  }
  
  return(list(k=K, gamma1=gamma1, kappa1=kappa1, beta=beta, Delta=delta.c))
}

CR_Worms <- function(Ytilde, Delta_Y) {
  n <- length(Ytilde)
  if (length(Ytilde)!=n) {
    stop("Xtilde and Ytilde need to be of the same length.")
  }
  
  if (length(Delta_Y)!=n) {
    stop("Ytilde and Delta_Y need to be of the same length.")
  }
  
  # Convert to logical vectors
  Delta_Y <- as.logical(Delta_Y)
  
  
  # Sort data according to Ytilde values
  s <- sort(Ytilde, index.return=TRUE)
  Ytilde_sort <- s$x
  Delta_Y_sort <- Delta_Y[s$ix]
  
  K <- 1:(n-1)
  H <- numeric(n)
  kappa <- numeric(n)
  
  for (k in 1:(n-1)) {
    
    prods <- numeric(k)
    # Compute products for i = 1,...,k
    for (i in 1:k) {
      #ind <- 1:(n-i)
      ind <- which(Ytilde_sort<Ytilde_sort[n-i+1])
      prods[i] <- prod((1-1/(n-ind+1))^(Delta_Y_sort[ind]-1))
    }
    i <- 1:k
    H[k] <- sum(i/n * prods * (log(Ytilde_sort[n-(1:k)+1])-log(Ytilde_sort[n-(1:k)]))) / (prods[k]*k/n)
  }
  return(list(k=K, H=H[K]))
}