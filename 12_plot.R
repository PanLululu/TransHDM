library(ggplot2)
###################################################################################################
########################################## DE IDE PE ##############################################
###################################################################################################

for(type in c('bias','MSE','SD')){
  for(rho in c(0,0.3,0.6)){
    for(p_m in c(1000,2000)){
      simu_index<-data.frame(simu_index)
      plotdata<-simu_index[simu_index$rho==rho & simu_index$p_m==p_m,
                           c(paste0('DE_',type),paste0('IDE_',type),paste0('PE_',type))]
      plotdata$x<-c(1:4)
      colnames(plotdata)<-c('DE','IDE','PE','x')
      plotdata<-sapply(plotdata, as.numeric)
      min<-min(plotdata[,1:3])
      max<-max(plotdata[,1:3])
      min_re<-min-(max-min)/10
      max_re<-max+(max-min)/10
      
      p<-ggplot(plotdata)+
        geom_line(aes(x,DE),color='#ab6380',size=0.6,linetype='dashed')+
        geom_point(aes(x,DE),color='#ab6380',shape=17,size=1.8)+
        
        geom_line(aes(x,IDE),color='#388dc2',size=0.6,linetype='dotdash')+
        geom_point(aes(x,IDE),color='#388dc2',shape=19,size=1.8)+
        
        geom_line(aes(x,PE),color='#5aa167',size=0.6,linetype='dotted')+
        geom_point(aes(x,PE),color='#5aa167',shape=18,size=2.6)+
        
        labs(title=paste0(type,'_rho_',rho,'_p_m_',p_m),x=" ", y=" ")+
        theme_grey()+
        theme(#panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          title = element_text(size=8,face='bold'),
          axis.text = element_text(size=9),
          axis.title.x = element_text(size=10,face='bold'),
          axis.title.y = element_text(size=10,face='bold'),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          #legend.box.margin=margin(-10,-10,-10,-10),
          legend.position="right",
          text=element_text(family="TimesNewRomanPSMT"))+
        scale_x_continuous(limits = c(0.8,4.2),breaks = seq(1,4,1))+
        scale_y_continuous(limits = c(min_re,max_re))
      
      tiff(paste0(type,'_rho_',rho,'_p_m_',p_m,".tiff"), 
           width = 7, height = 4, units = "cm", res = 500)
      print(p)
      dev.off()
    }
  }
}



for(type in c('Power','FWER')){
  for(p_m in c(1000,2000)){
    if(type=='Power'){min=0.83;max=1}else{min=0;max=0.31}
    plotdata<-data.frame(cbind(simu_index[simu_index$p_m==p_m & simu_index$rho==0,c(type)],
                               simu_index[simu_index$p_m==p_m & simu_index$rho==0.3,c(type)],
                               simu_index[simu_index$p_m==p_m & simu_index$rho==0.6,c(type)]))
    plotdata$x<-c(1:4)
    colnames(plotdata)<-c('rho0','rho3','rho6','x')
    plotdata<-sapply(plotdata, as.numeric)
    
    p<-ggplot(plotdata)+
      geom_line(aes(x,rho0),color='#f5944e',size=0.6,linetype='dashed')+
      geom_point(aes(x,rho0),color='#f5944e',shape=17,size=1.8)+
      
      geom_line(aes(x,rho3),color='#b777bd',size=0.6,linetype='dotdash')+
      geom_point(aes(x,rho3),color='#b777bd',shape=19,size=1.8)+
      
      geom_line(aes(x,rho6),color='#6ec4b3',size=0.6,linetype='dotted')+
      geom_point(aes(x,rho6),color='#6ec4b3',shape=18,size=2.6)+
      
      labs(title=paste0(type,'_p_m_',p_m),x=" ", y=" ")+
      theme_grey()+
      theme(#panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        title = element_text(size=8,face='bold'),
        axis.text = element_text(size=9),
        axis.title.x = element_text(size=10,face='bold'),
        axis.title.y = element_text(size=10,face='bold'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),  
        #legend.box.margin=margin(-10,-10,-10,-10),
        legend.position="right",
        text=element_text(family="TimesNewRomanPSMT"))+
      scale_x_continuous(limits = c(0.8,4.2),breaks = seq(1,4,1))+
      scale_y_continuous(limits = c(min,max))
    
    tiff(paste0(type,'_p_m_',p_m,".tiff"), 
         width = 6, height = 4, units = "cm", res = 600)
    print(p)
    dev.off()
  }
}


