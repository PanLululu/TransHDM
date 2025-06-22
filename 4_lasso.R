
# LASSO Function Documentation -------------------------------------------------

#' Fit LASSO Regression with Transfer Learning
#'
#' Fits a LASSO (Least Absolute Shrinkage and Selection Operator) regression model 
#' under a transfer learning framework. Supports feature selection and coefficient 
#' estimation by combining target data and source data.
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
#' @param transfer A logical value (default: FALSE) indicating whether to enable 
#'   transfer learning mode (combine source data with target data).
#' @param lambda A string (default: 'lambda.min') specifying the criterion for 
#'   selecting regularization parameter:
#'   \itemize{
#'     \item{'lambda.min': Lambda value that gives minimum cross-validation error}
#'     \item{'lambda.1se': Largest lambda value within 1 standard error of the minimum error}
#'   }
#'
#' @return A numeric vector \code{coef} containing LASSO coefficient estimates (including intercept).
#'
#' @examples
#' # Prepare target and source data
#' target <- list(x = matrix(rnorm(100 * 20), 100, 20), y = rnorm(100))
#' source <- list(x = matrix(rnorm(200 * 20), 200, 20), y = rnorm(200))
#'
#' # Non-transfer mode
#' coef_no_transfer <- lasso(target = target, transfer = FALSE, lambda = 'lambda.min')
#' print(coef_no_transfer)
#'
#' # Transfer learning mode
#' coef_transfer <- lasso(target = target, source = source, transfer = TRUE, lambda = 'lambda.min')
#' print(coef_transfer)
# ------------------------------------------------------------------------------

lasso<-function(target, source = NULL, transfer=FALSE,lambda='lambda.min'){
  
  x_t<-target$x
  y_t<-target$y
  x_s<-source$x
  y_s<-source$y
  
  if(!transfer){
    cv_glm<-cv.glmnet(x=as.matrix(cbind(1,x_t)), 
                      y=y_t, 
                      alpha =1, 
                      family="gaussian")
    
    coef<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm[[lambda]])],
            cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm[[lambda]])][-1])
  }else{
    cv_glm<-cv.glmnet(x=as.matrix(cbind(1,rbind(x_t,x_s))), 
                      y=c(y_t,y_s), 
                      alpha =1, 
                      family="gaussian")
    w<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm[[lambda]])],
         cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm[[lambda]])][-1])
    
    offset <- as.numeric(as.matrix(cbind(x_t)) %*% w[-1] + w[1])
    
    cv_glm<-cv.glmnet(x=as.matrix(cbind(1,x_t)), 
                      y=y_t, 
                      alpha =1, 
                      offset = offset,
                      family="gaussian")
    
    delta<-c(cv_glm$glmnet.fit$a0[which(cv_glm$lambda ==cv_glm[[lambda]])],
             cv_glm$glmnet.fit$beta[, which(cv_glm$lambda ==cv_glm[[lambda]])][-1])
    
    coef<-w+delta
  }
  return(coef=coef)
}

# Internal function: null_estimation
# A function to estimate the proportions of the three component nulls
# This is from HDMT package version < 1.0.4

null_estimation <- function(input_pvalues, lambda = 0.5) {
  ## input_pvalues is a matrix with 2 columns of p-values, the first column is p-value for exposure-mediator association, the second column is p-value for mediator-outcome association adjusted for exposure
  ## lambda is the threshold for pi_{00} estimation, default 0.5
  # check input
  if (is.null(ncol(input_pvalues))) {
    stop("input_pvalues should be a matrix or data frame")
  }
  if (ncol(input_pvalues) != 2) {
    stop("inpute_pvalues should have 2 column")
  }
  input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
  if (sum(stats::complete.cases(input_pvalues)) < nrow(input_pvalues)) {
    warning("input_pvalues contains NAs to be removed from analysis")
  }
  input_pvalues <- input_pvalues[stats::complete.cases(input_pvalues), ]
  if (!is.null(nrow(input_pvalues)) && nrow(input_pvalues) < 1) {
    stop("input_pvalues doesn't have valid p-values")
  }
  
  pcut <- seq(0.1, 0.8, 0.1)
  frac1 <- rep(0, 8)
  frac2 <- rep(0, 8)
  frac12 <- rep(0, 8)
  for (i in 1:8) {
    frac1[i] <- mean(input_pvalues[, 1] >= pcut[i]) / (1 - pcut[i])
    frac2[i] <- mean(input_pvalues[, 2] >= pcut[i]) / (1 - pcut[i])
    frac12[i] <- mean(input_pvalues[, 2] >= pcut[i] & input_pvalues[, 1] >= pcut[i]) / (1 - pcut[i])^2
  }
  
  ## use the median estimates for pi00 ##
  
  alpha00 <- min(frac12[pcut == lambda], 1)
  
  ## alpha1 is the proportion of nulls for first p-value
  ## alpha2 is the proportion of nulls for second p-value
  
  if (stats::ks.test(input_pvalues[, 1], "punif", 0, 1, alternative = "greater")$p > 0.05) alpha1 <- 1 else alpha1 <- min(frac1[pcut == lambda], 1)
  if (stats::ks.test(input_pvalues[, 2], "punif", 0, 1, alternative = "greater")$p > 0.05) alpha2 <- 1 else alpha2 <- min(frac2[pcut == lambda], 1)
  
  
  if (alpha00 == 1) {
    alpha01 <- 0
    alpha10 <- 0
    alpha11 <- 0
  } else {
    if (alpha1 == 1 & alpha2 == 1) {
      alpha01 <- 0
      alpha10 <- 0
      alpha11 <- 0
      alpha00 <- 1
    }
    
    if (alpha1 == 1 & alpha2 != 1) {
      alpha10 <- 0
      alpha11 <- 0
      alpha01 <- alpha1 - alpha00
      alpha01 <- max(0, alpha01)
      alpha00 <- 1 - alpha01
    }
    
    if (alpha1 != 1 & alpha2 == 1) {
      alpha01 <- 0
      alpha11 <- 0
      alpha10 <- alpha2 - alpha00
      alpha10 <- max(0, alpha10)
      alpha00 <- 1 - alpha10
    }
    
    if (alpha1 != 1 & alpha2 != 1) {
      alpha10 <- alpha2 - alpha00
      alpha10 <- max(0, alpha10)
      alpha01 <- alpha1 - alpha00
      alpha01 <- max(0, alpha01)
      
      if ((1 - alpha00 - alpha01 - alpha10) < 0) {
        alpha11 <- 0
        alpha10 <- 1 - alpha1
        alpha01 <- 1 - alpha2
        alpha00 <- 1 - alpha10 - alpha01
      } else {
        alpha11 <- 1 - alpha00 - alpha01 - alpha10
      }
    }
  }
  alpha.null <- list(alpha10 = alpha10, alpha01 = alpha01, alpha00 = alpha00, alpha1 = alpha1, alpha2 = alpha2)
  return(alpha.null)
}
