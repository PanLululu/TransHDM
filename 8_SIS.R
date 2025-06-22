# SIS Function Documentation ---------------------------------------------------

#' SIS Function for High-Dimensional Mediation Analysis
#'
#' The SIS function performs dimension reduction in high-dimensional mediation analysis 
#' using Sure Independence Screening (SIS) to select mediators most strongly associated 
#' with both the exposure and outcome variables. It supports transfer learning to enhance 
#' screening robustness and accuracy by leveraging source data with target data.
#'
#' @param target_data A data frame containing the target dataset with columns: 
#' outcome variable (Y), exposure (D), mediators (M1, M2, ...), and covariates (X1, X2, ...). 
#' Required.
#' @param source_data An optional data frame (default: NULL) containing source dataset 
#' with the same structure as target_data. Used when transfer = TRUE.
#' @param topN An integer (default: NULL) specifying the number of mediators to retain. 
#' If NULL, the number is automatically determined based on sample size and outcome type.
#' @param transfer A logical value (default: FALSE) indicating whether to use transfer learning 
#' (incorporating source data) for mediator screening.
#' @param verbose A logical value (default: TRUE) controlling console output of screening progress.
#' @param ncore An integer (default: 1) specifying the number of CPU cores for parallel computation.
#' @param dblasso_method A logical value (default: FALSE) Select the function for transfer learning. 
#' The default is to use least squares regression for transfer learning. Another option is to use dblasso.
#'
#' @return A list containing:
#' \itemize{
#'   \item target_SIS: Subset of target_data containing Y, D, selected mediators, and covariates
#'   \item source_SIS: Subset of source_data (if provided), same structure as target_SIS
#'   \item M_ID_name_SIS: Names of selected mediators
#' }
#'
#' @examples
#' # Prepare target and source data
#' target_data <- data.frame(
#'   Y = rnorm(200), D = rnorm(200), 
#'   matrix(rnorm(200*50), 200, 50, dimnames = list(NULL, paste0("M", 1:50))), 
#'   X1 = rnorm(200)
#' )
#' source_data <- data.frame(
#'   Y = rnorm(300), D = rnorm(300),
#'   matrix(rnorm(300*50), 300, 50, dimnames = list(NULL, paste0("M", 1:50))),
#'   X1 = rnorm(300)
#' )
#'
#' # Execute SIS function
#' result <- SIS(
#'   target_data = target_data,
#'   source_data = source_data,
#'   transfer = TRUE,
#'   ncore = 1
#' )
#' selected_mediators <- result$M_ID_name_SIS
#' print(selected_mediators)
# ------------------------------------------------------------------------------

SIS <- function(target_data,
                source_data=NULL,
                topN=NULL,
                transfer=FALSE,
                verbose=T,
                ncore=1,
                dblasso_method=F
){
  if(verbose) message("Step 1: Sure Independent Screening ...", "  (", format(Sys.time(), "%X"), ")")
  
  # the number of top mediators that associated with exposure (X): d_0
  p_x<-length(grep("X", colnames(target_data), value = TRUE)) # dimension of X: p_m
  p_m<-length(grep("M", colnames(target_data), value = TRUE)) # dimension of M: p_m
  n_t<-nrow(target_data) # sample size of target_data
  
  if (is.null(topN)){
    d_0 <- ceiling(2 * n_t/log(n_t))
  }else{
    d_0 <- topN
  }
  d_0 <- min(p_m, d_0) # if d > p select all mediators
  
  ### 1. Estimate the regression coefficients beta (mediators --> outcome)
  beta_SIS <- matrix(0, 1, p_m)
  MXD_t<-target_data[,c(paste0('M',1:p_m),paste0('X',1:p_x),'D')]
  Y_t<-target_data[,'Y']
  
  if(transfer){
    
    if(!dblasso_method){
      MXD_s<-source_data[,c(paste0('M',1:p_m),paste0('X',1:p_x),'D')]
      MXD_ts<-rbind(MXD_t,MXD_s)
      Y_ts<-c(target_data[,'Y'],source_data[,'Y'])
      
      if(ncore==1){
        for (i in 1:p_m) {
          MXD_SIS_ts <- MXD_ts[, c(i, (p_m + 1):(p_m + p_x + 1))]
          MXD_SIS_t <- MXD_t[, c(i, (p_m + 1):(p_m + p_x + 1))]
          w <- glm(Y_ts~.,data=cbind(MXD_SIS_ts, Y_ts),family = "gaussian")$coefficients
          offset <- as.vector(as.matrix(MXD_SIS_t) %*% w[-1] + w[1])
          cv_glm<-cv.glmnet(x=as.matrix(cbind(1,MXD_SIS_t)),
                            y=Y_t,
                            alpha =0,
                            offset = offset,
                            family="gaussian")
          delta<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm$lambda.min)],
                   cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm$lambda.min)][-1])
          
          beta_SIS[i] <- (w+delta)[2]
        }
      }else{
        doParallel::registerDoParallel(ncore)
        system.time(beta_SIS<-foreach(i = 1:p_m, .combine = "c",.packages=c("glmnet","MASS","foreach","caret")) %dopar% {
          MXD_SIS_ts <- MXD_ts[, c(i, (p_m + 1):(p_m + p_x + 1))]
          MXD_SIS_t <- MXD_t[, c(i, (p_m + 1):(p_m + p_x + 1))]
          w <- glm(Y_ts~.,data=cbind(MXD_SIS_ts, Y_ts),family = "gaussian")$coefficients
          offset <- as.vector(as.matrix(MXD_SIS_t) %*% w[-1] + w[1])
          cv_glm<-cv.glmnet(x=as.matrix(cbind(1,MXD_SIS_t)), 
                            y=Y_t, 
                            alpha =0, 
                            offset = offset,
                            family="gaussian")
          delta<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm$lambda.min)],
                   cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm$lambda.min)][-1])
          
          (w+delta)[2]
        })
      }
      
    }else{
      MXD_s<-source_data[,c(paste0('M',1:p_m),paste0('X',1:p_x),'D')]
      Y_s<-source_data[,'Y']
      
      if(ncore==1){
        for (i in 1:p_m) {
          beta_SIS[i]<-dblasso(target = list(x=MXD_t[, c(i, (p_m + 1):(p_m + p_x + 1))],y=Y_t),
                               source = list(x=MXD_s[, c(i, (p_m + 1):(p_m + p_x + 1))],y=Y_s),
                               transfer=T)$dbcoef.hat[2]
        }
      }else{
        doParallel::registerDoParallel(ncore)
        system.time(beta_SIS<-foreach(i = 1:p_m, .combine = "c",.packages=c("glmnet","MASS","foreach","caret")) %dopar% {
          dblasso(target = list(x=MXD_t[, c(i, (p_m + 1):(p_m + p_x + 1))],y=Y_t),
                  source = list(x=MXD_s[, c(i, (p_m + 1):(p_m + p_x + 1))],y=Y_s),
                  transfer=T)$dbcoef.hat[2]
        })
      }
    }
   
  }else{
    
    if(!dblasso_method){
      for (i in 1:p_m) {
        MXD_SIS <- MXD_t[, c(i, (p_m + 1):(p_m + p_x + 1))]
        fit <- glm(Y_t~.,data=cbind(MXD_SIS, Y_t),family = "gaussian")
        beta_SIS[i] <- fit$coefficients[2]
      }
      
    }else{
      for (i in 1:p_m) {
        beta_SIS[i]<-dblasso(target = list(x=MXD_t[, c(i, (p_m + 1):(p_m + p_x + 1))],y=Y_t))$dbcoef.hat[2]
      }
    }
  }
  
  ### 2. Estimate the regression coefficients alpha (exposure --> mediators)
  alpha_SIS <- matrix(0, 1, p_m)
  DX_t<-target_data[,c('D',paste0('X',1:p_x))]
  M_t<-target_data[,paste0('M',1:p_m)]
  
  if(transfer){
    
    if(!dblasso_method){
      
      DX_s<-source_data[,c('D',paste0('X',1:p_x))]
      DX_ts<-rbind(DX_t,DX_s)
      M_ts<-rbind(target_data[,paste0('M',1:p_m)],source_data[,paste0('M',1:p_m)])
      
      if(ncore==1){
        for (i in 1:p_m) {
          fit <- lsfit(DX_ts, M_ts[, i], intercept = TRUE)
          w <- coef(fit)
          offset <- as.vector(as.matrix(DX_t) %*% w[-1] + w[1])
          cv_glm<-cv.glmnet(x=as.matrix(cbind(1,DX_t)),
                            y=M_t[, i],
                            alpha =0,
                            offset = offset,
                            family="gaussian")
          delta<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm$lambda.min)],
                   cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm$lambda.min)][-1])
          
          alpha_SIS[i] <- (w+delta)[2]
        }
      }else{
        system.time(alpha_SIS<-foreach(i = 1:p_m, .combine = "c",.packages=c("glmnet","MASS","foreach","caret")) %dopar% {
          fit <- lsfit(DX_ts, M_ts[, i], intercept = TRUE)
          w <- coef(fit)
          offset <- as.vector(as.matrix(DX_t) %*% w[-1] + w[1])
          cv_glm<-cv.glmnet(x=as.matrix(cbind(1,DX_t)), 
                            y=M_t[, i], 
                            alpha =0, 
                            offset = offset,
                            family="gaussian")
          delta<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm$lambda.min)],
                   cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm$lambda.min)][-1])
          
          (w+delta)[2]
        })
        closeAllConnections()
      }
      
    }else{
      DX_s<-source_data[,c('D',paste0('X',1:p_x))]
      M_s<-source_data[,paste0('M',1:p_m)]
      
      if(ncore==1){
        for (i in 1:p_m) {
          
          alpha_SIS[i]<-dblasso(target = list(x=DX_t,y=M_t[, i]),
                                source = list(x=DX_s,y=M_s[, i]),
                                transfer=T)$dbcoef.hat[2]
        }
      }else{
        system.time(alpha_SIS<-foreach(i = 1:p_m, .combine = "c",.packages=c("glmnet","MASS","foreach","caret")) %dopar% {
          dblasso(target = list(x=DX_t,y=M_t[, i]),
                  source = list(x=DX_s,y=M_s[, i]),
                  transfer=T)$dbcoef.hat[2]
        })
        closeAllConnections()
      }
    }
  
  }else{
    
    if(!dblasso_method){
      for (i in 1:p_m) {
        fit <- lsfit(DX_t, M_t[, i], intercept = TRUE)
        alpha_SIS[i] <- matrix(coef(fit))[2]
      }
    }else{
      for (i in 1:p_m) {
        alpha_SIS[i]<-dblasso(target = list(x=DX_t,y=M_t[, i]))$dbcoef.hat[2]
      }
    }
  }
  
  # 3. Select the d_0 number of mediators with top largest effect
  ab_SIS <- alpha_SIS * beta_SIS
  ID_SIS <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d_0])
  d <- length(ID_SIS)
  M_ID_name_SIS <- names(target_data)[grep('M',names(target_data))[ID_SIS]]
  target_SIS<-target_data[,c('Y','D',M_ID_name_SIS,paste0('X',1:p_x))]
  source_SIS<-source_data[,c('Y','D',M_ID_name_SIS,paste0('X',1:p_x))]
  
  if(verbose) message("Top ", d, " mediators are selected: ", paste0(M_ID_name_SIS, collapse = ", "), "  (", format(Sys.time(), "%X"), ")")
  
  return(list(target_SIS=target_SIS,source_SIS=source_SIS,M_ID_name_SIS=M_ID_name_SIS))
}

