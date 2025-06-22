
generate_simulationData<-function(n=100,
                                  p_x=5,
                                  rho=0,
                                  p_m=100,
                                  h=2, #0 1 2 3
                                  source=F,
                                  transferable=T,
                                  seed=NULL){
  set.seed(seed)
  gamma<-rep(0.3,p_x)
  alpha0<-0
  #alpha1<-c(rep(0.5,p_m_true+1),rep(0,p_m-p_m_true-1)) # D
  alpha1<-c(c(0.6,0.55,0.65,0.5,0.35,0.3),0,0,rep(0,p_m-8)) # D
  alpha2<-sapply(1:p_m,function(x){rep(0.3,p_x)}) # X
  
  beta0<-0
  beta1<-1 # D
  #beta2<-c(rep(0.2,p_m_true),0,0.2,rep(0,p_m-p_m_true-2)) # M
  beta2<-c(c(0.4,0.45,0.45,0.5),0,0,c(0.35,0.25),rep(0,p_m-8))
  beta3<-rep(0,p_m) # DM
  beta4<-rep(0.3,p_x) # X
  
  if(source){
    if(transferable){
      set.seed(123)
      alpha1<-c(c(0.6,0.55,0.65,0.5,0.35,0.3),0,0,sample(c(rep(0.3,h),rep(0,p_m-8-h)))) 
      set.seed(123)
      beta2<-c(c(0.4,0.45,0.45,0.5),0,0,c(0.35,0.25),sample(c(rep(0.3,h),rep(0,p_m-8-h))))
      set.seed(NULL)
    }else{
      set.seed(123)
      alpha1<-c(sample(c(rep(0.5,10+2),0,0,rep(0,p_m-14))))
      set.seed(123)
      beta2<-c(sample(c(rep(0.5,10),0,0,rep(0.5,2),rep(0,p_m-14))))
      set.seed(NULL)
    }
  }
  
  ### confounder: X
  X.Sigma=diag(rep(1,p_x))
  for(i in 1:p_x){for(j in 1:p_x){X.Sigma[i,j]=0.5^abs(i-j)}}
  X<-mvrnorm(n, c(rep(0,p_x)), X.Sigma)
  
  ### treatment: D
  #D <- ifelse((X%*%gamma+rnorm(n,0,1))>0.5,1,0)
  #D <- rbinom(n, 1, plogis(X%*%gamma))
  D <- X%*%gamma + rnorm(n,0,1)
  
  ### mediator
  M.Sigma=diag(rep(1,p_m))
  for(i in 1:p_m){for(j in 1:p_m){M.Sigma[i,j]=rho^abs(i-j)}}
  M.se<-mvrnorm(n,rep(0,p_m), M.Sigma)
  #M.se<-M.se[,sample(1:ncol(M.se))]
  #set.seed(123)
  M<-matrix(data=NA,nrow=n,ncol=p_m)
  for(i in 1:p_m){
    M[,i] <- alpha1[i]*D+X%*%alpha2[,i]+M.se[,i]
  }
  
  ### outcome
  #Y<-rbinom(n,1,plogis(beta0 + beta1*D + M%*%beta2 + X%*%beta4))
  Y<-beta0 + beta1*D + M%*%beta2 + X%*%beta4 + rnorm(n,0,1)
  
  data<-as.data.frame(cbind(Y=Y,D=D,M=M,X=X))
  colnames(data)<-c('Y','D',paste0('M',1:p_m),paste0('X',1:p_x))
  
  ### true effect
  coef<-list(alpha0=unname(alpha0),
             alpha1=unname(alpha1),
             alpha2=unname(alpha2),
             beta0=unname(beta0),
             beta1=unname(beta1),
             beta2=unname(beta2),
             beta3=unname(beta3),
             beta4=unname(beta4))
  
  set.seed(NULL)
  
  return(list(data=data,coef=coef))
}

