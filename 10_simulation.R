source('~/M')

library(caret)
library(MASS)
library(glmnet)
library(parallel)
library(doParallel)
library(foreach)
library(qvalue)
library(doParallel)

cl<-makeCluster(50)
registerDoParallel(cl)

p_x=5
n_simu=200
n_s=200

for(p_m in c(1000,2000)){
  for(n in c(100)){
    for(rho in c(0,0.3,0.6)){
      
      message('rho: ',rho)
      message('n: ',n)
      message('p_m: ',p_m)
      
      true_coef<-generate_simulationData(n=n,p_x=p_x,p_m=p_m,rho=rho)$coef
      true_effect<-true_coef$beta2*true_coef$alpha1
      
      system.time(simu_HDMA<-foreach(s = 1:n_simu, .combine = "rbind",.packages=c("glmnet","MASS","foreach","caret","qvalue",'HDMT')) %dopar% {
        target_data<-generate_simulationData(n=n,p_x=p_x,p_m=p_m,rho=rho,seed=s)$data
        HDMA.fit<-TransHDMA(target_data,verbose = T)
        c(HDMA.fit$DE_est,HDMA.fit$IDE_est)
      })
      
      simu<-list(simu_HDMA=simu_HDMA,true_effect=true_effect)
      save(simu,file = paste0('simu_HDMA_rho_',rho,'_n_',n,'_p_m_',p_m,'.RData'))
      message('simu_HDMA end')
      
      for(s_n in c(1,2,3)){
        system.time(simu_TransHDMA<-foreach(s = 1:n_simu, .combine = "rbind",.packages=c("glmnet","MASS","foreach","caret","qvalue",'HDMT')) %dopar% {
          target_data<-generate_simulationData(n=n,p_x=p_x,p_m=p_m,rho=rho,seed=s)$data
          source_data<-foreach(k= 1:s_n, .combine = "rbind") %do% {
            generate_simulationData(n=n_s,p_x=p_x,p_m=p_m,rho=rho,source = T,transferable = T,h=h)$data
          }

          TransHDMA.fit<-TransHDMA(target_data=target_data,source_data,T,verbose = T,ncore =1)
          c(TransHDMA.fit$DE_est,TransHDMA.fit$IDE_est)
        })
        
        simu<-list(simu_TransHDMA=simu_TransHDMA,true_effect=true_effect)
        save(simu,file = paste0('simu_TransHDMA_rho_',rho,'_n_',n,'_p_m_',p_m,'_s_n_',s_n,'.RData'))
        
        message('s_n: ',s_n)
        message('simu_TransHDMA end')
      }
    }
  }
}

