simu_index<-c()

for(rho in c(0,0.3,0.6)){
  for(n in c(100)){
    for(p_m in c(1000,2000)){
      for(s_n in c(NA,200,400,600)){
        if(is.na(s_n)){
          method='simu_HDMA'
          load(paste0('simu_HDMA_rho_',rho,'_n_',n,'_p_m_',p_m,'.RData'))
        }else{
          method='simu_TransHDMA'
          load(paste0('simu_TransHDMA_rho_',rho,'_n_',n,'_p_m_',p_m,'_n_s_',s_n,'.RData'))
        }
        result<-simu[[method]]
        n_simu<-dim(result)[1]
        true_effect<-simu$true_effect
        
        # true_effect
        DE_true<-1
        IDE_true<-sum(true_effect[1:4])
        PE_true<-IDE_true/(DE_true+IDE_true)
        
        IDE_perM_true<-c(mean(true_effect[1:4]),rep(0,3))
        PE_perM_true<-IDE_perM_true/(DE_true+IDE_true)
        
        # direct effect
        DE_est<-result[,'D']
        DE_bias<-abs(mean(DE_est)-DE_true)/DE_true
        DE_MSE<-sqrt(mean((DE_est-DE_true)^2))
        DE_sd<-sd(DE_est)
        
        # indirect effect
        IDE_est<-rowSums(result[,-1])
        IDE_bias<-abs(mean(IDE_est)-IDE_true)/IDE_true
        IDE_MSE<-sqrt(mean((IDE_est-IDE_true)^2))
        IDE_SD<-sd(IDE_est)
        
        # indirect effect per M type
        ### M 1-4: alpha beta !=0 
        ### M 5-6: alpha =0  beta !=0
        ### M 7-8: alpha !=0 beta =0
        ### M 9-p_m: alpha beta =0
        IDE_perM_est<-cbind(rowMeans(result[,paste0('M',1:4)]),
                            rowMeans(result[,paste0('M',5:6)]),
                            rowMeans(result[,paste0('M',7:8)]),
                            rowMeans(result[,paste0('M',9:p_m)]))
        IDE_perM_bias<-abs(colMeans(IDE_perM_est)-IDE_perM_true)
        IDE_perM_MSE<-sqrt(colMeans(t(t(IDE_perM_est)-IDE_perM_true)^2))
        IDE_perM_SD<-apply(IDE_perM_est,2,sd)
        
        # proportion effect
        PE_est<-IDE_est/(IDE_est+DE_est) 
        PE_bias<-abs(mean(PE_est)-PE_true)/PE_true
        PE_MSE<-sqrt(mean((PE_est-PE_true)^2))
        PE_SD<-sd(PE_est)
        
        # proportion effect per M type
        PE_perM_est<-cbind(rowMeans(result[,paste0('M',1:4)]),
                           rowMeans(result[,paste0('M',5:6)]),
                           rowMeans(result[,paste0('M',7:8)]),
                           rowMeans(result[,paste0('M',9:p_m)]))/(IDE_est+DE_est)
        PE_perM_bias<-abs(colMeans(PE_perM_est)-PE_perM_true)
        PE_perM_MSE<-sqrt(colMeans(t(t(PE_perM_est)-PE_perM_true)^2))
        PE_perM_SD<-apply(PE_perM_est,2,sd)
        
        # Power for M1-M10
        Power<-length(which(result[,paste0('M',1:4)]>0))/length(result[,paste0('M',1:4)])
        
        # T1E per 2-4 M type ：被错误纳入的无关变量数/总无关变量数
        T1E_perM<-c(
          sum(rowSums(result[,paste0('M',5:6)]>0))/n_simu/2,
          sum(rowSums(result[,paste0('M',7:8)]>0))/n_simu/2,
          sum(rowSums(result[,paste0('M',9:p_m)]>0))/n_simu/(p_m-4-4)
        )
        
        # FWER
        FWER<-length(which(rowSums(result[,paste0('M',5:p_m)])!=0))/n_simu
        
        index<-c(rho,n,p_m,s_n,method,
                 DE_bias,DE_MSE,DE_sd,
                 IDE_bias,IDE_MSE,IDE_SD,
                 PE_bias,PE_MSE,PE_SD,
                 Power,FWER,
                 IDE_perM_bias,IDE_perM_MSE,IDE_perM_SD,
                 PE_perM_bias,PE_perM_MSE,PE_perM_SD,
                 T1E_perM)
        names(index)<-c('rho','n','p_m','s_n','method',
                        'DE_bias','DE_MSE','DE_SD',
                        'IDE_bias','IDE_MSE','IDE_SD',
                        'PE_bias','PE_MSE','PE_SD',
                        'Power','FWER',
                        paste0('IDE_perM_bias_',1:4),paste0('IDE_perM_MSE_',1:4),paste0('IDE_perM_SD_',1:4),
                        paste0('PE_perM_bias_',1:4),paste0('PE_perM_MSE_',1:4),paste0('PE_perM_SD_',1:4),
                        paste0('T1E_perM_',2:4))
        
        simu_index<-rbind(simu_index,index)
        
      }
    }
  }
}
