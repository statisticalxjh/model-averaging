DGP1_C1 <- function(n,p,R2,rho){
  M1 <- matrix(rep(1:p,p),ncol=p,byrow=F)
  corM <- rho^(abs(M1-t(M1)))
  X <- mvrnorm(n,rep(0,p),corM)
  
   # beta <- 1/(1:p)
    beta<- sort(1/(1:p))
 
  ############################################case1
  epsilon=rnorm(n,0,1)
  aa<- 1
   bb<- t(beta)%*%corM%*%beta
  theta<- sqrt(aa*R2/((1-R2)*bb))
  
  Y <-as.vector(theta)*X%*%as.vector(beta)+epsilon
  delta=cbind(matrix(1,n,3),outer(1*(X[,1]< 0),c(1,1,1),"*"),
              outer(1*(X[,2]< 0),c(1,1,1),"*"),outer(1*(X[,3]< 0),c(1,1,1),"*"))
  
  res <- list(Y=Y,X=X,delta=delta)
  return(res)
}
 

DGP1_C2<- function(n,p,R2,rho){
  M1 <- matrix(rep(1:p,p),ncol=p,byrow=F)
  corM <- rho^(abs(M1-t(M1)))
  X <- mvrnorm(n,rep(0,p),corM)
  X <- as.matrix(X)
  # beta <- 1/(1:p)
 beta<- sort(1/(1:p))
 
  
  
  ############################################Case2
  
  epsilon=rt(n,3)
  aa<- 3
  bb<- t(beta)%*%corM%*%beta
  theta<- sqrt(aa*R2/((1-R2)*bb))
  
  ###################################################
  
  Y <-as.vector(theta)*X%*%as.vector(beta)+epsilon
  delta=cbind(matrix(1,n,3),outer(1*(X[,1]< 0),c(1,1,1),"*"),
              outer(1*(X[,2]< 0),c(1,1,1),"*"),outer(1*(X[,3]< 0),c(1,1,1),"*"))
  
  res <- list(Y=Y,X=X,delta=delta)
  return(res)
}


DGP1_C3<- function(n,p,R2,rho){
  M1 <- matrix(rep(1:p,p),ncol=p,byrow=F)
  corM <- rho^(abs(M1-t(M1)))
  X <- mvrnorm(n,rep(0,p),corM)
  X <- as.matrix(X)
   #beta <- 1/(1:p)
 beta<- sort(1/(1:p))
 
  ############################################Case3
  epsilon=rlaplace(n,0,1)
  aa<- 2
   bb<- t(beta)%*%corM%*%beta
  theta<- sqrt(aa*R2/((1-R2)*bb))
  
  

  
  Y <-as.vector(theta)*X%*%as.vector(beta)+epsilon
   delta=cbind(matrix(1,n,3),outer(1*(X[,1]< 0),c(1,1,1),"*"),
              outer(1*(X[,2]< 0),c(1,1,1),"*"),outer(1*(X[,3]< 0),c(1,1,1),"*"))
  
  res <- list(Y=Y,X=X,delta=delta)
  return(res)
}


DGP1_C4 <- function(n,p,R2,rho){
  M1 <- matrix(rep(1:p,p),ncol=p,byrow=F)
  corM <- rho^(abs(M1-t(M1)))
  X <- mvrnorm(n,rep(0,p),corM)
  X <- as.matrix(X)
   #beta <- 1/(1:p)
 beta<- sort(1/(1:p))
 
  
  ############################################Case4
  
  epsilon=(rchisq(n,3)-3)/sqrt(6)
  aa<- 1
  bb<- t(beta)%*%corM%*%beta
  theta<- sqrt(aa*R2/((1-R2)*bb))
  ###################################################
  
  Y <-as.vector(theta)*X%*%as.vector(beta)+epsilon
  delta=cbind(matrix(1,n,3),outer(1*(X[,1]< 0),c(1,1,1),"*"),
              outer(1*(X[,2]< 0),c(1,1,1),"*"),outer(1*(X[,3]< 0),c(1,1,1),"*"))
  
  res <- list(Y=Y,X=X,delta=delta)
  return(res)
}


Checkloss <- function (u,tau){
  return(u*(tau-(u<=0)))
}


 

a0=c(1,1,1)
a1=c(0,0,0)
a2=c(0,0,1)
a3=c(0,1,0)
a4=c(0,1,1)
a5=c(1,0,0)
a6=c(1,0,1)
a7=c(1,1,0)
pattern1=rbind(rep(a0,each=3),rep(a1,each=3),rep(a2,each=3),rep(a3,each=3),
               rep(a4,each=3),rep(a5,each=3),rep(a6,each=3),rep(a7,each=3))
pattern=cbind(matrix(1,8,3),pattern1)
