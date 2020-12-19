rm(list=ls())
library(quantreg)
library(cqrReg)
library("Rglpk")
library("truncdist")
library("rmutil")
library("MASS")
source("fc.R")

set.seed(123)
n <- 300
p <- 12
M<-8
rho=0.5  
R2_all=(1:9)/10

outloop<-length(R2_all)
res_final<-matrix(0,outloop,7)

for (km in 1:outloop){
  
  loop<-200
  result<-matrix(0,loop,7)
  for (jjj in 1:loop){
    
    data<-DGP1_C1(n,p,R2_all[km],rho)
    Y<-data$Y
    X<-data$X
    delta<-data$delta
    
    SS_m=list()
    S_m=list()
    for (m in 1:M){
      Delta_m=pattern[m,]
      Index=c()
      for (i in 1:n){
        if (all(delta[i,]==Delta_m))
        {Index=union(Index,i)}
      }
      SS_m[[m]]=Index
      S_m[[m]]=union(SS_m[[m]],SS_m[[1]])
    }
    
    tauseq=(1:5)/6
    
    mu<-list()
    for (m in 1:M){
      index_m<-which(pattern[m,]==1)
      XX_m=X[S_m[[1]],index_m]
      YY_m=Y[S_m[[1]]]
      
      aa=matrix(0,length(YY_m),length(tauseq))
      for (i in 1:length(YY_m)){
        temp_regressor <- XX_m
        tempone <- temp_regressor[i,]
        temp_regressor <- temp_regressor[-i,]
        tempY <- YY_m;
        tempY <- tempY[-i]
        thetahat_CV <- cqr.fit(temp_regressor,tempY,tau <- tauseq,method="mm")
        beta_hat <- thetahat_CV$beta
        b_hat <- thetahat_CV$b
        aa[i,]=rep(tempone%*%beta_hat,length=length(tauseq))+b_hat
      }
      mu[[m]]=aa
    }
    
    
    matc<-cbind(matrix(mu[[1]],length(tauseq)*length(YY_m),1),
                matrix(mu[[2]],length(tauseq)*length(YY_m),1),
                matrix(mu[[3]],length(tauseq)*length(YY_m),1),
                matrix(mu[[4]],length(tauseq)*length(YY_m),1),
                matrix(mu[[5]],length(tauseq)*length(YY_m),1),
                matrix(mu[[6]],length(tauseq)*length(YY_m),1),
                matrix(mu[[7]],length(tauseq)*length(YY_m),1),
                matrix(mu[[8]],length(tauseq)*length(YY_m),1))
    
    obj<-c(rep(0,length=M),rep(tauseq,each=length(YY_m)),
           rep(1-tauseq,each=length(YY_m)))
    
    mat1<-cbind(matc,diag(1,length(tauseq)*length(YY_m)),
                diag(-1,length(tauseq)*length(YY_m)))
    mat2<-cbind(diag(1,M),matrix(0,M,2*length(tauseq)*length(YY_m)))
    mat3<-mat2
    mat4<-cbind(matrix(0,length(tauseq)*length(YY_m),M),
                diag(1,length(tauseq)*length(YY_m)),
                diag(0,length(tauseq)*length(YY_m)))
    mat5<-cbind(matrix(0,length(tauseq)*length(YY_m),M),
                diag(0,length(tauseq)*length(YY_m)),
                diag(1,length(tauseq)*length(YY_m)))
    mat6<-c(rep(1,length=M),rep(0,length=2*length(tauseq)*length(YY_m)))
    
    mat<-rbind(mat1,mat2,mat3,mat4,mat5,mat6)
    
    dir <- c(rep("==",length=length(tauseq)*length(YY_m)),rep(">=",length <- M),
             rep("<=",length<- M),rep(">=",length <-length(tauseq)*length(YY_m)),
             rep(">=",length <-length(tauseq)*length(YY_m)),"==")
    
    rhs <- c(rep(YY_m,length <- length(tauseq)),rep(0,length <- M),rep(1,length <- M),
             rep(0,length=2*length(tauseq)*length(YY_m)),1)
    
    est <- Rglpk_solve_LP(obj, mat, dir, rhs)$solution
    w<-est[1:M]
    
    
    
    ###################################################################################
    ############## AIC ######################################################
    
    X_AIC<-X[S_m[[1]],]
    Y_AIC<-Y[S_m[[1]]]
    n_AIC<-length(Y_AIC)
    
    AICm<-rep(0,length<- M)
    BICm<-rep(0,length<- M)
    for (mA in 1:M){
      index_AIC<-which(pattern[mA,]==1)
      X_AIC_m <- X[S_m[[1]],index_AIC]
      aa_AIC <- cqr.lasso.mm(X_AIC_m,Y_AIC,tauseq)
      beta_hat_AIC <- aa_AIC$beta
      b_hat_AIC <-  aa_AIC$b
      
      Loss_AIC <-rep(0,length <- length(tauseq))
      for (kc in 1:length(tauseq)){
        Loss_AIC[kc] <-
          mean(Checkloss(Y_AIC-X_AIC_m%*%beta_hat_AIC-b_hat_AIC[kc],tauseq[kc]))
      }
      AICm[mA]=2*n_AIC*length(tauseq)*log(mean(Loss_AIC))+2*(mA+length(tauseq))
      BICm[mA]=2*n_AIC*length(tauseq)*log(mean(Loss_AIC))+log(n_AIC)*(mA+length(tauseq))
    }
    w_AIC<-exp(-0.5*AICm)/sum(exp(-0.5*AICm))
    w_BIC<-exp(-0.5*BICm)/sum(exp(-0.5*BICm))
    
    
    ##################################################################### 
    ########################### Prediction ##############################
    
    n_T=10000
    data_test<-DGP1_C1(n_T,p,R2_all[km],rho)
    Y_test1<-data_test$Y
    X_test1<-data_test$X
    delta_test<-data_test$delta
    index_cc<-which(apply(delta_test,1,sum)==p)
    Y_test=Y_test1[index_cc]
    X_test=X_test1[index_cc,]
    
    MatTemp<-matrix(0,length(Y_test),length(tauseq))
    MatTemp_AIC<-matrix(0,length(Y_test),length(tauseq))
    MatTemp_BIC<-matrix(0,length(Y_test),length(tauseq))
    MatTemp_1<-matrix(0,length(Y_test),length(tauseq))
    MatTemp_equal<-matrix(0,length(Y_test),length(tauseq))
    for (mm in 1:M){
      index_mm<-which(pattern[mm,]==1)
      X_m<-X[S_m[[mm]],index_mm]
      Y_m <- Y[S_m[[mm]]]
      resTemp<-cqr.fit(X_m,Y_m,tau <- tauseq,method="mm")
      beta_hat_m<-resTemp$beta
      b_hat_m<-resTemp$b
      
      aaa=X_test[,index_mm]%*%beta_hat_m
      aaa1=outer(as.vector(aaa),rep(1,length=length(tauseq)),"*")
      aaa2=t(outer(b_hat_m,rep(1,nrow(aaa)),"*"))
      MatTemp1<-(aaa1+aaa2)*w[mm]
      MatTemp <- MatTemp+MatTemp1
      
      MatTemp2<-(aaa1+aaa2)*w_AIC[mm]
      MatTemp_AIC <- MatTemp_AIC+MatTemp2
      
      MatTemp3<-(aaa1+aaa2)*w_BIC[mm]
      MatTemp_BIC <- MatTemp_BIC+MatTemp3
      
      MatTemp4<-(aaa1+aaa2)*(1/M)
      MatTemp_equal<- MatTemp_equal+MatTemp4
    }
    
    
    Loss<-c(0,length=length(tauseq))
    Loss_AIC<-c(0,length=length(tauseq))
    Loss_BIC<-c(0,length=length(tauseq))
    Loss_equal<-c(0,length=length(tauseq))
    for (kk in 1:length(tauseq)){
      Loss[kk]=mean(Checkloss(Y_test-MatTemp[,kk],tauseq[kk]))
      Loss_AIC[kk]=mean(Checkloss(Y_test-MatTemp_AIC[,kk],tauseq[kk]))
      Loss_BIC[kk]=mean(Checkloss(Y_test-MatTemp_BIC[,kk],tauseq[kk]))
      Loss_equal[kk]=mean(Checkloss(Y_test-MatTemp_equal[,kk],tauseq[kk]))
    }
    
    pre1<-mean(Loss)
    pre5<-mean(Loss_AIC)
    pre6<-mean(Loss_BIC)
    pre7<-mean(Loss_equal)
    
    ######################################################################
    ######################### CC ####################################
    
    XX_cc=X[S_m[[1]],]
    YY_cc=Y[S_m[[1]]]
    resTemp_cc<-cqr.fit(XX_cc,YY_cc,tau <- tauseq,method="mm")
    beta_hat_m_cc<-resTemp_cc$beta
    b_hat_m_cc<-resTemp_cc$b
    
    Loss_cc<-c(0,length=length(tauseq))
    for (kc in 1:length(tauseq)){
      Loss_cc[kc]=mean(Checkloss(Y_test-X_test%*%beta_hat_m_cc-b_hat_m_cc[kc]
                                 ,tauseq[kc]))
    }
    
    pre2<-mean(Loss_cc)
    
    ##########################################################
    ################# G1 #############################
    
    XX_G1=X[,1:3]
    YY_G1=Y
    resTemp_G1<-cqr.fit(XX_G1,YY_G1,tau <- tauseq,method="mm")
    beta_hat_m_G1<-resTemp_G1$beta
    
    b_hat_m_G1<-resTemp_G1$b
    Loss_G1<-c(0,length=length(tauseq))
    for (kc in 1:length(tauseq)){
      Loss_G1[kc]=mean(Checkloss(Y_test-X_test[,1:3]%*%beta_hat_m_G1-b_hat_m_G1[kc]
                                 ,tauseq[kc]))
    }
    pre3<-mean(Loss_G1)
    
    
    ##########################################################################
    ################## LASSO ###################################################
    X_lasso<-X[S_m[[1]],]
    Y_lasso<-Y[S_m[[1]]]
    aa_lasso <- cqr.lasso.mm(X_lasso,Y_lasso,tauseq)
    beta_hat_lasso <- aa_lasso$beta
    b_hat_lasso <-  aa_lasso$b
    
    Loss_lasso<-c(0,length=length(tauseq))
    for (kc in 1:length(tauseq)){
      Loss_lasso[kc]=mean(Checkloss(Y_test-X_test%*%beta_hat_lasso-b_hat_lasso[kc]
                                    ,tauseq[kc]))
    }
    pre4<-mean(Loss_lasso)
    ######################################################################
    ######################################################################
    
    
    result[jjj,]<-c(pre1,pre5,pre6,pre2,pre3,pre4,pre7)
    print(jjj)
  }
  
  res_final[km,]=apply(result,2,median)
  print(km)
}
 

result.final<-data.frame(res_final)
colnames(result.final)<-c("Proposed","AIC","BIC","CC","G1","LASSO","equal")



 