#' TransHDM: High-Dimensional Mediation Analysis with Transfer Learning
#'
#' The \code{TransHDM} function performs high-dimensional mediation analysis under a transfer learning framework. 
#' It identifies and estimates indirect effects of mediators between exposure (\code{D}) and outcome (\code{Y}) 
#' by integrating target data and optional source data.
#'
#' @param target_data A data frame containing the target dataset with:
#'   \itemize{
#'     \item \code{Y}: Outcome variable (numeric vector)
#'     \item \code{D}: Exposure variable (numeric vector)
#'     \item \code{M1, M2, ...}: Mediator variables (numeric columns)
#'     \item \code{X1, X2, ...}: Covariates (numeric columns)
#'   }
#'   \emph{(Required)}
#' @param source_data A data frame (optional, default: \code{NULL}) containing source dataset with the same 
#'   structure as \code{target_data}. Used when \code{transfer = TRUE}.
#' @param transfer A logical value (default: \code{FALSE}) indicating whether to enable transfer learning mode.
#' @param verbose A logical value (default: \code{TRUE}) controlling progress message display in console.
#' @param ncore An integer (default: 1) specifying number of CPU cores for parallel computation.
#' @param topN An integer (default: \code{NULL}) specifying the number of mediators to retain after SIS screening.
#'   If \code{NULL}, this will be auto-calculated based on data dimensions.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{contributions}: Data frame of significant mediators with columns:
#'     \itemize{
#'       \item \code{alpha}: Exposure-mediator coefficients
#'       \item \code{beta}: Mediator-outcome coefficients
#'       \item \code{alpha_beta}: Indirect effects (αβ products)
#'       \item \code{ab_pv}: P-values for indirect effects
#'       \item \code{pa}: Proportion mediated
#'     }
#'   \item \code{contributions_select}: Data frame of all post-SIS mediators with additional FDR-adjusted values
#'   \item \code{effects}: Data frame summarizing:
#'     \itemize{
#'       \item Indirect effect (Total)
#'       \item Direct effect (\code{DE_est})
#'       \item Total effect
#'       \item Proportion mediated
#'     }
#'   \item \code{IDE_est}: Numeric vector of indirect effect estimates for all mediators (non-selected mediators set to 0)
#'   \item \code{DE_est}: Direct effect estimate of exposure on outcome
#' }
#'
#' @examples
#' # Prepare target and source data
#' target_data <- data.frame(
#'   Y = rnorm(200), 
#'   D = rnorm(200),
#'   matrix(rnorm(200*50), 200, 50, 
#'   dimnames = list(NULL, paste0("M", 1:50)),
#'   X1 = rnorm(200)
#' )
#' 
#' source_data <- data.frame(
#'   Y = rnorm(300),
#'   D = rnorm(300),
#'   matrix(rnorm(300*50), 300, 50,
#'   dimnames = list(NULL, paste0("M", 1:50)),
#'   X1 = rnorm(300)
#' )
#'
#' # Run analysis with transfer learning
#' result <- TransHDM(
#'   target_data = target_data,
#'   source_data = source_data,
#'   transfer = TRUE,
#'   ncore = 2,
#'   topN = 10
#' )
#' 
#' print(result$contributions)
#' print(result$effects)
#' 
#' @export


TransHDM<-function(
    target_data,
    source_data=NULL,
    transfer=FALSE,
    verbose=TRUE,
    ncore=1,
    topN=NULL
){
  p_m<-length(grep('M',colnames(target_data)))
 
  ### Step-1 Sure Independent Screening
  SIS_result<-SIS(target_data,source_data,transfer = transfer,verbose=verbose,ncore=ncore,topN=topN)
  target_SIS<-SIS_result$target_SIS
  source_SIS<-SIS_result$source_SIS
  
  target<-list(x=target_SIS[,-1],y=target_SIS[,1])
  source<-list(x=source_SIS[,-1],y=source_SIS[,1])
  
  ### Step-2  De-biased Lasso Estimates         
  if(verbose) message("Step 2: De-biased Lasso Estimates  ...","   (", format(Sys.time(), "%X"), ")")
  
  M_list<-colnames(target_SIS)[grepl('M',colnames(target_SIS))]
  p_m_SIS<-length(M_list)
  p_x<-length(colnames(target_SIS)[grepl('X',colnames(target_SIS))])
  
  # 2.1. Estimate beta 
  DBLASSO_fit<-dblasso(target,source,transfer)
  DE_est<-DBLASSO_fit$dbcoef.hat['D']
  beta_SIS_est <- DBLASSO_fit$dbcoef.hat[M_list]
  P_beta_SIS <- DBLASSO_fit$P.value[M_list]
  
  # if(length(which(P_beta_SIS<0.05))==0)
  #   stop("No mediatior is identified !")
  
  if(verbose)  message("Non-zero dblasso estimate(s) of mediator(s) found from Outcome model: ", 
                       paste0(names(which(P_beta_SIS<0.05)), collapse = ","))
  
  # 2.2 Estimate alpha 
  alpha_SIS_est <- rep(0, p_m_SIS)
  P_alpha_SIS <- rep(0,p_m_SIS)
  
  if(ncore==1){
    for (i in 1:p_m_SIS) {
      target<-list(x=target_SIS[,c('D',paste0('X',1:p_x))],y=target_SIS[,M_list[i]])
      source<-list(x=source_SIS[,c('D',paste0('X',1:p_x))],y=source_SIS[,M_list[i]])
      DBLASSO_fit<-dblasso(target,source,transfer)
      P_alpha_SIS[i] <- DBLASSO_fit$P.value['D'] ## the SIS for alpha
      alpha_SIS_est[i] <- DBLASSO_fit$dbcoef.hat['D']
    }
  }else{
    doParallel::registerDoParallel(ncore)
    system.time(alpha_SIS<-foreach(i = 1:p_m_SIS, .combine = "rbind",.packages=c("glmnet","MASS","foreach","caret","qvalue",'HDMT'),.export = c("dblasso","lasso")) %dopar% {
      target<-list(x=target_SIS[,c('D',paste0('X',1:p_x))],y=target_SIS[,M_list[i]])
      source<-list(x=source_SIS[,c('D',paste0('X',1:p_x))],y=source_SIS[,M_list[i]])
      DBLASSO_fit<-dblasso(target,source,transfer)
      c(DBLASSO_fit$P.value['D'],DBLASSO_fit$dbcoef.hat['D'])
    })
    P_alpha_SIS<-alpha_SIS[,1]
    alpha_SIS_est<-alpha_SIS[,2]
    closeAllConnections()
  }
  
  names(P_alpha_SIS)<-M_list
  names(alpha_SIS_est)<-M_list
  
  # if(length(which(P_alpha_SIS<0.05))==0)
  #   stop("No mediatior is identified !")
  
  if(verbose) message("Non-zero dblasso estimate(s) of mediator(s) found from Mediator model: ", 
                      paste0(names(which(P_alpha_SIS<0.05)), collapse = ","))
  
  ### STEP 3   The multiple-testing  procedure         
  if(verbose) message("Step 3: The multiple-testing  procedure ...", "   (", format(Sys.time(), "%X"), ")")
  
  PA <- cbind(c(P_alpha_SIS), c(P_beta_SIS))
  P_value <- apply(PA, 1, max) # The joint p-values for SIS variable
  N0 <- dim(PA)[1] * dim(PA)[2]
  input_pvalues <- PA + matrix(runif(N0, 0, 10^{-10}), dim(PA)[1], 2)
  
  # Estimate the proportions of the three component nulls
  nullprop <- null_estimation(input_pvalues)
  
  # Compute the estimated pointwise FDR for every observed p-max
  fdrcut <- HDMT::fdr_est(nullprop$alpha00,
                          nullprop$alpha01,
                          nullprop$alpha10,
                          nullprop$alpha1,
                          nullprop$alpha2,
                          input_pvalues,
                          exact = 0
  )
  FDRcut=0.05
  ID_fdr <- which(fdrcut <= FDRcut)
  # Following codes extract the estimates for mediators with fdrcut<=0.05
  beta_hat_est <- beta_SIS_est[ID_fdr]
  alpha_hat_est <- alpha_SIS_est[ID_fdr]
  P.value <- P_value[ID_fdr]
  
  # Indirect effect
  alpha_beta_hat_est <- beta_hat_est * alpha_hat_est # mediation(indirect) effect
  
  IDE_est<-rep(0,p_m)
  IDE_est[as.numeric(gsub('M','',names(alpha_beta_hat_est)))]<-alpha_beta_hat_est
  names(IDE_est)<-paste0('M',1:p_m)
  
  TE_est<-sum(IDE_est)+DE_est
  
  PE_est<-sum(IDE_est)/TE_est
  PE_est_perM<-alpha_beta_hat_est/TE_est
  
  effects<-data.frame(effect=c('indirect','direct','total','pe'),
                      estimate=c(sum(IDE_est),DE_est,TE_est,PE_est))
  
  contributions_select<-data.frame(
    mediator=names(alpha_SIS_est),
    alpha=alpha_SIS_est,
    alpha_pv=P_alpha_SIS,
    beta=beta_SIS_est,
    beta_pv=P_beta_SIS,
    alpha_beta=alpha_SIS_est*beta_SIS_est,
    pe=alpha_SIS_est*beta_SIS_est/TE_est,
    fdr=fdrcut
  )
  
  contributions<-data.frame(
    mediator=names(alpha_beta_hat_est),
    alpha=alpha_hat_est,
    alpha_pv=P_alpha_SIS[ID_fdr],
    beta=beta_hat_est,
    beta_pv=P_beta_SIS[ID_fdr],
    alpha_beta=alpha_beta_hat_est,
    ab_pv=P.value,
    pa=alpha_beta_hat_est/TE_est
  )
  
  
  if(verbose) message("Non-zero dblasso estimate(s) of mediator(s) found: ", paste0(names(alpha_beta_hat_est), collapse = ","))
  
  return(list(contributions=contributions,
              contributions_select=contributions_select,
              effects=effects,IDE_est=IDE_est,DE_est=DE_est))
}

