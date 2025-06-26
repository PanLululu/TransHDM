# --------------------------- Load All Custom Functions --------------------------- #
# These R scripts define core functions for simulation and estimation

source("1_generate_simulationData_homogeneous_design.R")
source("2_generate_simulationData_heterogeneous_design.R")
source("3_kfold_split.R")
source("4_lasso.R")
source("5_dblasso.R")
source("6_null_estimation.R")
source("7_source_detection.R")
source("8_SIS.R")
source("9_TransHDM.R")

# --------------------------- Load Required Libraries --------------------------- #
library(caret)
library(MASS)
library(glmnet)
library(parallel)
library(doParallel)
library(foreach)
library(qvalue)

# --------------------------- Initialize Parallel Backend --------------------------- #
cl <- makeCluster(50)             # Create a cluster with 50 cores
registerDoParallel(cl)            # Register the cluster for parallel processing

# --------------------------- Fixed Simulation Parameters --------------------------- #
p_x <- 5                          # Number of covariates
n_simu <- 200                     # Number of repetitions per condition
n_s <- 200                        # Sample size for each source dataset

# --------------------------- Simulation Loops Over Conditions --------------------------- #
for (p_m in c(1000, 2000)) {              # Number of mediators
  for (n in c(100)) {                     # Target sample size
    for (rho in c(0, 0.3, 0.6)) {         # Correlation level between mediators
      
      message('rho: ', rho)
      message('n: ', n)
      message('p_m: ', p_m)
      
      # Generate true effect coefficients based on structural model
      true_coef <- generate_simulationData(n = n, p_x = p_x, p_m = p_m, rho = rho)$coef
      true_effect <- true_coef$beta2 * true_coef$alpha1
      
      # ------------------- Simulation for HDM (without transfer learning) ------------------- #
      system.time(
        simu_HDMA <- foreach(s = 1:n_simu, .combine = "rbind",
                             .packages = c("glmnet", "MASS", "foreach", "caret", "qvalue", "HDMT")) %dopar% {
                               # Generate target data
                               target_data <- generate_simulationData(n = n, p_x = p_x, p_m = p_m, rho = rho, seed = s)$data
                               
                               # Fit HDM model
                               HDM.fit <- TransHDM(target_data, verbose = TRUE)
                               
                               # Return estimates of direct and indirect effects
                               c(HDM.fit$DE_est, HDM.fit$IDE_est)
                             }
      )
      
      # Save HDM results
      simu <- list(simu_HDMA = simu_HDMA, true_effect = true_effect)
      save(simu, file = paste0('simu_HDM_rho_', rho, '_n_', n, '_p_m_', p_m, '.RData'))
      message('simu_HDM end')
      
      # ---------------- Simulation for TransHDM (with transfer learning) ---------------- #
      for (s_n in c(1, 2, 3)) {  # Vary the number of source datasets
        
        system.time(
          simu_TransHDMA <- foreach(s = 1:n_simu, .combine = "rbind",
                                    .packages = c("glmnet", "MASS", "foreach", "caret", "qvalue", "HDMT")) %dopar% {
                                      # Generate target data
                                      target_data <- generate_simulationData(n = n, p_x = p_x, p_m = p_m, rho = rho, seed = s)$data
                                      
                                      # Generate s_n source datasets with transferable structure
                                      source_data <- lapply(1:3,function(k){
                                        if(k<=s_n){
                                          generate_simulationData(n = n_s, p_x = p_x, p_m = p_m, rho = rho,
                                                                  source = TRUE, transferable = TRUE, h = 2)$data
                                        }else{
                                          generate_simulationData(n = n_s, p_x = p_x, p_m = p_m, rho = rho,
                                                                  source = TRUE, transferable = FALSE, h = 2)$data
                                        }
                                      }) 
                                      
                                      # detect transferable source
                                      detect<-source_detection(target_data,source_data,kfold=5)
                                      
                                      source_data <- foreach(k=detect$transfer.source.id, .combine = "rbind") %do% {
                                        source_data[[k]]
                                      }
                                      
                                      # Fit TransHDM model with transfer learning enabled
                                      TransHDM.fit <- TransHDM(target_data = target_data,
                                                               source_data = source_data,
                                                               transferable = TRUE,
                                                               verbose = TRUE,
                                                               ncore = 1)
                                      
                                      # Return effect estimates
                                      c(TransHDM.fit$DE_est, TransHDM.fit$IDE_est)
                                    }
        )
        
        # Save TransHDM results
        simu <- list(simu_TransHDMA = simu_TransHDMA, true_effect = true_effect)
        save(simu, file = paste0('simu_TransHDM_rho_', rho, '_n_', n, '_p_m_', p_m, '_s_n_', s_n, '.RData'))
        
        message('s_n: ', s_n)
        message('simu_TransHDM end')
      }
    }
  }
}

# --------------------------- Cleanup --------------------------- #
stopCluster(cl)  # Shut down the parallel backend after simulation is complete
