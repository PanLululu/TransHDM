library(dplyr)
library(ggplot2)

###########################################################################################
###############################     Covariate      ########################################
###########################################################################################

# Import baseline demographic data including the following variables: 
# 'RID', 'AGE', 'B_DATE', 'research_group', 'PTGENDER', 'PTRACCAT', 'PTEDUCAT', 'PTDOBYY', and 'PHASE'

DEM<-read.csv('1-data/Covariate/Covariate_baseline.csv')

###########################################################################################
#############################     exposure：APOE      #####################################
###########################################################################################

### GENOTYPE
APOE<-read.csv('1-data/exposure/APOE.csv')

DEM_X<-merge(DEM,
             APOE[,c('RID','GENOTYPE')], 
             by = "RID", all = TRUE)

###########################################################################################
############################     outcome: question      ###################################
###########################################################################################

#################################### ADAS #################################################

#################################### ADNI 1
ADAS1<-read.csv('1-data/Outcome/ADAS/ADAS1.csv')
ADAS1$TOTAL13<-ADAS1$TOTALMOD

ADAS1 <- ADAS1 %>%
  mutate(EXAMDATE = as.Date(EXAMDATE))

# Eliminate missing
ADAS1$TOTAL13[ADAS1$TOTAL13<0]<-NA

ADAS1 <- ADAS1 %>%
  arrange(RID, desc(EXAMDATE)) %>% 
  group_by(RID) %>% 
  filter(!is.na(TOTAL13)) %>% 
  slice(1) %>% 
  ungroup()

#################################### ADNI GO 23
ADAS23GO<-read.csv('1-data/Outcome/ADAS/ADASGO23.csv')
ADAS23GO$EXAMDATE<-ADAS23GO$VISDATE

ADAS23GO <- ADAS23GO %>%
  mutate(EXAMDATE = as.Date(EXAMDATE))

ADAS23GO <- ADAS23GO %>%
  arrange(RID, desc(EXAMDATE)) %>% 
  group_by(RID) %>% 
  filter( !is.na(TOTAL13)) %>%
  slice(1) %>% 
  ungroup()

#################################### merge

ADAS<-merge(ADAS1[,c('RID','EXAMDATE','TOTAL13')],
            ADAS23GO[,c('RID','EXAMDATE','TOTAL13')], 
            by = "RID", all = TRUE)

ADAS <- ADAS %>%
  rowwise() %>%
  mutate(
    TOTAL13 = case_when(
      !is.na(TOTAL13.x) & !is.na(TOTAL13.y) ~ ifelse(
        as.Date(EXAMDATE.y) >= as.Date(EXAMDATE.x), 
        TOTAL13.y, 
        TOTAL13.x
      ),
      !is.na(TOTAL13.x) & is.na(TOTAL13.y) ~ TOTAL13.x,
      is.na(TOTAL13.x) & !is.na(TOTAL13.y) ~ TOTAL13.y,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

ADAS<-as.data.frame(ADAS)
ADAS$EXAMDATE<-as.Date(apply(cbind(ADAS$EXAMDATE.x,ADAS$EXAMDATE.y),1,function(x){max(x,na.rm = T)}))
ADAS$DATE_TOTAL13<-ADAS$EXAMDATE

DEM_XY<-merge(DEM_X,
              ADAS[,c("RID",'DATE_TOTAL13','TOTAL13')],
              by = "RID", all = TRUE)

################################## biomarker ##############################################

### TAU
tau<-read.csv('1-data/Outcome/biomarker/TAU.csv')

tau <- tau %>%
  mutate(EXAMDATE = as.Date(EXAMDATE))

tau <- tau %>%
  arrange(RID, desc(EXAMDATE)) %>% 
  group_by(RID) %>%
  filter(!is.na(TAU) | row_number() == 1) %>% 
  slice(1) %>% 
  ungroup()

tau$DATE_TAU<-tau$EXAMDATE

DEM_XY<-merge(DEM_XY,
              tau[,c("RID",'DATE_TAU','TAU')], 
              by = "RID", all = TRUE)


### AB42
AB42<-read.csv('1-data/Outcome/biomarker/AB42.csv')

AB42 <- AB42 %>%
  mutate(EXAMDATE = as.Date(EXAMDATE))

AB42 <- AB42 %>%
  arrange(RID, desc(EXAMDATE)) %>%
  group_by(RID) %>% 
  filter(!is.na(ABETA42) | row_number() == 1) %>% 
  slice(1) %>% 
  ungroup()

AB42$DATE_AB42<-AB42$EXAMDATE

DEM_XY<-merge(DEM_XY,
              AB42[,c("RID",'DATE_AB42','ABETA42')], 
              by = "RID", all = TRUE)

################################# brain volumn #############################################

white<-read.csv('1-data/Outcome/WMH/WMH.csv')
white <- white %>%
  mutate(EXAMDATE = as.Date(EXAMDATE))

white <- white %>%
  arrange(RID, desc(EXAMDATE)) %>%
  group_by(RID) %>%
  filter(!is.na(TOTAL_WMH) | row_number() == 1) %>% 
  slice(1) %>%
  ungroup()
white$DATE_WHM<-white$EXAMDATE

DEM_XY<-merge(DEM_XY,
              white[,c("RID",'DATE_WHM','TOTAL_WMH')], 
              by = "RID", all = TRUE)

###########################################################################################
###################################     lipid      ########################################
###########################################################################################

M<-read.csv('1-data/Lipid/lipids.csv')
M <- M %>% distinct(RID, .keep_all = TRUE)

intersect_id<-intersect(DEM_XY$RID,M[,'RID'])
position_id1<-match(intersect_id,DEM_XY$RID)
position_id2<-match(intersect_id,M[,'RID'])

data<-cbind(DEM_XY[position_id1,],M[position_id2,])

# Exposure classification
data$D<-0
data[data$GENOTYPE %in% c('2/4','3/4','4/4'),'D']<-1
data[data$GENOTYPE %in% c('4/4'),'D']<-2

# Inclusion and Exclusion
data<-data[which(!is.na(data$GENOTYPE)),]
data<-data[!is.na(data$RID) & !is.na(data$AGE),] 
data<-data[!is.na(data$PTEDUCAT),] 

### Take ADAS ending as an example
data<-data[which(!is.na(data$TOTAL13)),]
data<-data[data$EXAMDATE<=data$DATE_TOTAL13,]


##########################################################################################
###############################     TransHDM      ########################################
##########################################################################################
library(caret)
library(MASS)
library(glmnet)
library(qvalue)
library(HDMT)
library(doParallel)
library(foreach)
library(readxl)

D.names<-'D'
Y.names<-'TOTAL13'
X.names<-c('AGE','PTEDUCAT','PTGENDER')
M.names<-colnames(data)[match('SPH.D18.1.',colnames(data)):(match('SPH.D18.1.',colnames(data))+780)] 

data<-data[,c(Y.names,D.names,M.names,X.names,'PTRACCAT')]
colnames(data)<-c('Y','D',paste0('M',1:length(M.names)),paste0('X',1:length(X.names)),'PTRACCAT')

# Data segmentation
target_data<-data[data$PTRACCAT==4,c('Y','D',paste0('M',1:length(M.names)),paste0('X',1:length(X.names)))]
source_data<-data[data$PTRACCAT %in% c(5),c('Y','D',paste0('M',1:length(M.names)),paste0('X',1:length(X.names)))]

# Standardization
target_data[,c(paste0('M',1:length(M.names)))] <- target_data[,c(paste0('M',1:length(M.names)))] %>%
  mutate(across(where(is.numeric), ~ scale(.)[, 1])) 
source_data[,c(paste0('M',1:length(M.names)))] <- source_data[,c(paste0('M',1:length(M.names)))] %>%
  mutate(across(where(is.numeric), ~ scale(.)[, 1]))

# Identification of transferability
source_detection_result<-source_detection(target_data,list(source_data))
source_detection_result

# Modeling
### HDMA in target
black.fit<-TransHDMA(target_data,verbose = T,topN = 60,ncore = 5)

### TransHDMA
trans.fit<-TransHDMA(target_data,source_data,transfer=T,verbose = T,topN =150,ncore = 5,dblasso_SIS=T)




