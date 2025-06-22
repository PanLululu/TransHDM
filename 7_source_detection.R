# Source Detection Function Documentation ----------------------------------------

#' Detect Transferable Source Data via Cross-Validation
#' 
#' Determines whether external source datasets can be effectively transferred to the target data 
#' by comparing predictive performance using LASSO regression under a transfer learning framework.
#'
#' @param target_data A data frame containing the target dataset. Required columns: 
#'   \itemize{
#'     \item{Y: Outcome variable (first column)}
#'     \item{D: Exposure variable (second column)}
#'     \item{M1, M2, ...: Mediators (columns with 'M' prefix)}
#'     \item{X1, X2, ...: Covariates (columns with 'X' prefix)}
#'   }
#' @param source_data A list of data frames containing source datasets. Each element must have 
#'   the same structure as \code{target_data}.
#' @param kfold Integer (default: 5). Number of folds for cross-validation.
#' @param C0 Numeric (default: 1). Threshold constant for determining transferability. 
#'   Larger values make the criterion more lenient.
#'
#' @return A list containing:
#' \itemize{
#'   \item{transfer.source.id: Indices of transferable source datasets}
#'   \item{source.loss: Mean validation loss for each source dataset}
#'   \item{target.valid.loss: Mean validation loss using target-only model}
#'   \item{threshold: Calculated transferability threshold (target.valid.loss + C0*sd)}
#'   \item{loss.cv: Full k-fold cross-validation loss matrix}
#' }
#' 
#' @examples
#' # Generate synthetic data
#' target_data <- data.frame(
#'   Y = rnorm(200), D = rnorm(200),
#'   M1 = rnorm(200), M2 = rnorm(200),
#'   X1 = rnorm(200)
#' )
#' 
#' source1 <- data.frame(
#'   Y = rnorm(300), D = rnorm(300),
#'   M1 = rnorm(300), M2 = rnorm(300),
#'   X1 = rnorm(300)
#' )
#' source2 <- data.frame(
#'   Y = rnorm(250), D = rnorm(250),
#'   M1 = rnorm(250), M2 = rnorm(250),
#'   X1 = rnorm(250)
#' )
#' 
#' # Run source detection
#' result <- source_detection(
#'   target_data = target_data,
#'   source_data = list(source1, source2),
#'   kfold = 5,
#'   C0 = 1.5
#' )
#' 
#' # Display transferable sources
#' cat("Transferable source IDs:", result$transfer.source.id, "\n")
#' 
#' # Compare loss values
#' print(data.frame(
#'   Source = c(paste0("Source", 1:2), "Target"),
#'   Loss = c(result$source.loss, result$target.valid.loss),
#'   Threshold = result$threshold
#' ))
# ------------------------------------------------------------------------------

source_detection<-function(target_data,source_data,kfold=5,C0=1){
  
  target_data_split<-kfold_split(target_data,kfold=kfold)
  M_list<-colnames(target_data)[grepl('M',colnames(target_data))]
  
  loss<-function(data,coef_beta){
    n<-nrow(data)
    p_m<-length(grep('M',names(coef_beta)))
    p_x<-length(grep('X',names(coef_beta)))
    
    mu<-as.matrix(cbind(1,data[,c('D',paste0('M',1:p_m),paste0('X',1:p_x))]))%*%as.matrix(coef_beta)
    residuals<-data[,'Y']-mu
    Sigma <-(t(residuals)%*%residuals)/n
    l.Y<-n/2*log(2*pi)+n/2*log(Sigma)+n/2
    
    return(l.Y)
  }
  
  loss.cv<-t(sapply(1:kfold,function(k){
    
    train<-target_data_split$train_set[[k]]
    test<-target_data_split$test_set[[k]]
    
    source.loss<-sapply(1:length(source_data), function(j){
      #data_all<-list(x=rbind(train[,-1],source_data[[j]][,-1]),y=c(train[,1],source_data[[j]][,1]))
      train_list<-list(x=train[,-1],y=train[,1])
      source_data_list<-list(x=source_data[[j]][,-1],y=source_data[[j]][,1])
      coef_beta<-lasso(train_list,source_data_list,transfer = T)
      loss(test,coef_beta)
    })
    
    train_list<-list(x=train[,-1],y=train[,1])
    coef_beta<-lasso(train_list)
    target.loss<-loss(test,coef_beta)
    c(source.loss, target.loss)
  }))
  
  source.loss <- colMeans(loss.cv)[1:(ncol(loss.cv)-1)]
  target.valid.loss <- colMeans(loss.cv)[ncol(loss.cv)]
  target.valid.loss.sd <- sd(loss.cv[, ncol(loss.cv)])
  
  threshold <- target.valid.loss + C0*max(target.valid.loss.sd, 0.01)
  transfer.source.id <- which(source.loss <= threshold)
  
  obj <- list(transfer.source.id = transfer.source.id, source.loss = source.loss, target.valid.loss = target.valid.loss,
              threshold = threshold,loss.cv=loss.cv)
  obj
}



