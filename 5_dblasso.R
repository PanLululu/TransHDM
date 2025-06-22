
# dblasso Function Documentation -------------------------------------------------

#' Fit Debiased LASSO with Transfer Learning
#'
#' Fits a debiased LASSO regression model under transfer learning framework, supporting feature 
#' selection and coefficient estimation by combining target and source data.
#'
#' @param target A list containing two elements:
#'   \itemize{
#'     \item{x: Feature matrix of target data (numeric matrix)}
#'     \item{y: Response vector of target data (numeric vector)}
#'   } Required.
#' @param source A list (optional, default: NULL) containing two elements:
#'   \itemize{
#'     \item{x: Feature matrix of source data (numeric matrix)}
#'     \item{y: Response vector of source data (numeric vector)}
#'   } Used when transfer = TRUE.
#' @param transfer A logical value (default: FALSE) indicating whether to enable transfer learning 
#'   (combining source data with target data for estimation).
#' @param level A numeric value (default: 0.95) specifying confidence level for confidence intervals.
#' @param lambda A string specifying criterion for selecting regularization parameter:
#'   \itemize{
#'     \item{'lambda.min': Lambda value giving minimum cross-validation error}
#'     \item{'lambda.1se': Largest lambda within 1 standard error of minimum error}
#'   }
#'
#' @return A list containing:
#' \itemize{
#'   \item{dbcoef.hat: Debiased LASSO coefficient vector (including intercept)}
#'   \item{coef.hat: Original LASSO coefficient vector}
#'   \item{CI: Data frame with confidence intervals (lb = lower bound, ub = upper bound)}
#'   \item{var.est: Variance estimates for debiased coefficients}
#'   \item{se.est: Standard errors for debiased coefficients}
#'   \item{P.value: Vector of p-values for coefficients}
#' }
#'
#' @examples
#' # Prepare target and source data
#' target <- list(x = matrix(rnorm(100 * 20), 100, 20), y = rnorm(100))
#' source <- list(x = matrix(rnorm(200 * 20), 200, 20), y = rnorm(200))
#'
#' # Non-transfer mode
#' result_no_transfer <- dblasso(target = target, transfer = FALSE, 
#'                              level = 0.95, lambda = 'lambda.min')
#' print(result_no_transfer$dbcoef.hat)
#' print(result_no_transfer$CI)
#'
#' # Transfer learning mode
#' result_transfer <- dblasso(target = target, source = source, transfer = TRUE,
#'                           level = 0.95, lambda = 'lambda.min')
#' print(result_transfer$dbcoef.hat)
#' print(result_transfer$P.value)
# ------------------------------------------------------------------------------

dblasso<- function(target, source = NULL, transfer=FALSE, level = 0.95,lambda='lambda.1se') {
  
  gamma_est<-function(target,source=NULL,transfer=FALSE){
    p_x<-dim(target$x)[2]
    
    if(!transfer){
      target_data<-as.data.frame(cbind(target$y,target$x))
      colnames(target_data)<-c('Y',paste0('X',1:p_x))
      
      cv_glm<-cv.glmnet(x=as.matrix(cbind(target_data[,c(paste0('X',1:p_x))])), 
                        y=target_data[,'Y'],
                        alpha =1,
                        intercept=F,
                        family="gaussian")
      beta_TL<-cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm[[lambda]])]
      
    }else{
      
      target_data<-as.data.frame(cbind(target$y,target$x))
      source_data<-as.data.frame(cbind(source$y,source$x))
      colnames(target_data)<-c('Y',paste0('X',1:p_x))
      colnames(source_data)<-c('Y',paste0('X',1:p_x))
      
      data_all<-rbind(target_data,source_data)
      
      cv_glm<-cv.glmnet(x=as.matrix(cbind(data_all[,c(paste0('X',1:p_x))])), 
                        y=data_all[,'Y'],
                        alpha =1,
                        intercept=F,
                        family="gaussian")
      w<-cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm[[lambda]])]
      
      offset <- as.numeric(as.matrix(target_data[,c(paste0('X',1:p_x))]) %*% w)
      
      # Step2 de-bias rough estimator
      cv_glm<-cv.glmnet(x=as.matrix(cbind(target_data[,c(paste0('X',1:p_x))])), 
                        y=target_data[,'Y'],
                        alpha =1,
                        intercept=F,
                        offset = offset,
                        family="gaussian")
      
      delta<-cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm[[lambda]])]
      beta_TL<-w+delta
    }
    
    return(beta_TL)
  }
  
  # confidence interval: r.level
  r.level <- level + (1-level)/2
  j <- 0
  
  # lasso estimator
  coef.hat<-lasso(target=target, source=source, transfer=transfer)
  coef.name<-names(coef.hat)
  coef.hat <- as.vector(coef.hat)
  
  # data list: D
  
  D <- list(target = target, source = source)
  D$target$x <- cbind(1, as.matrix(D$target$x))
  if(transfer){D$source$x <- cbind(1, as.matrix(D$source$x))}
  
  # dimension of design matrix
  p <- ncol(D$target$x)
  
  # combination of design matrix
  if(transfer){
    X.comb <- rbind(D$source$x,D$target$x)
  }else{
    X.comb <- D$target$x
  }
  
  # Sigma.hat of design matrix
  Sigma.hat <- t(X.comb) %*% X.comb/nrow(X.comb)
  
  ### inverse fisher information matrix: Theta.hat
  L <- foreach(j = 1:p, .combine = "rbind") %do% {
    D1 <- D
    D1$target$y <- D1$target$x[, j]
    D1$target$x <- D1$target$x[, -j]
    
    if(transfer){
      D1$source$y <- D1$source$x[, j]
      D1$source$x <- D1$source$x[, -j]
      gamma<-gamma_est(target = D1$target, source = D1$source, transfer=T)
    }else{
      gamma<-gamma_est(target = D1$target)
    }
    
    tau2 <- Sigma.hat[j, j] - Sigma.hat[j, -j, drop = F] %*% gamma
    theta <- rep(1, p)
    theta[-j] <- -gamma
    c(theta, tau2)
  }
  
  Theta.hat <- solve(diag(L[,p+1])) %*% L[1:p, 1:p]
  
  ### debiased lasso estimator: dbcoef.hat
  residuals <- D$target$y - D$target$x %*% coef.hat
  dbcoef.hat <- as.vector(as.matrix(coef.hat) + Theta.hat %*% t(D$target$x) %*% residuals/length(D$target$y))
  
  ### variance estimator
  var.est <- diag(Theta.hat %*% Sigma.hat %*% t(Theta.hat))
  se.est <- sqrt(var.est/length(D$target$y))
  CI <- data.frame(lb = dbcoef.hat - qnorm(r.level)*se.est, 
                   ub = dbcoef.hat + qnorm(r.level)*se.est)
  Z<-(dbcoef.hat-0)/se.est
  P.value<-as.vector(2 * pnorm(abs(Z), lower.tail = FALSE))
  
  names(coef.hat)<-names(dbcoef.hat)<-rownames(CI)<-names(var.est)<-names(se.est)<-names(P.value)<-coef.name
  
  return(list(dbcoef.hat = dbcoef.hat, coef.hat = coef.hat, CI = CI, var.est = var.est, se.est=se.est,P.value=P.value))
}
