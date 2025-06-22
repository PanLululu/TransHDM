# kfold_split Function Documentation -------------------------------------------

#' K-Fold Cross-Validation Data Splitting
#' 
#' Splits input data into k folds for cross-validation, generating training and test sets 
#' for each fold. Particularly useful for mediator selection stability assessment 
#' in high-dimensional mediation analysis.
#'
#' @param data A data frame or matrix containing the dataset to be split. 
#'            Rows represent observations, columns represent variables.
#' @param kfold Integer (default: 3). Number of folds for cross-validation. 
#'            Must be ≥2 and ≤nrow(data).
#' @param seed Integer (optional). Random seed for reproducible splitting. 
#'            Default NULL produces non-deterministic splits.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item{train_set: List of length kfold, each element is a training subset}
#'   \item{test_set: List of length kfold, each element is the corresponding test subset}
#' }
#' 
#' @examples
#' # Create sample mediation analysis data
#' mediation_data <- data.frame(
#'   Y = rnorm(300),
#'   D = rnorm(300),
#'   M1 = rnorm(300), M2 = rnorm(300),
#'   X1 = rnorm(300)
#' )
#'
#' # Create 5-fold split with reproducibility
#' splits <- kfold_split(mediation_data, kfold = 5, seed = 123)
# ------------------------------------------------------------------------------

kfold_split<-function(data,kfold=3,seed=NULL){ 
  set.seed(seed) 
  n=dim(data)[1] 
  split_indices <- createFolds(1:n, k = kfold, list = TRUE) 
  train_set<-list() 
  test_set<-list() 
  for(k in 1:kfold){ 
    train_set[[k]] <- data[-split_indices[[k]], ] 
    test_set[[k]] <- data[split_indices[[k]], ] 
    
  } 
  return(list(train_set=train_set,test_set=test_set)) 
}