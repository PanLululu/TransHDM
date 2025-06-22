# Initialize an empty matrix to store evaluation metrics across all settings
simu_index <- c()

# Loop over combinations of rho, n, p_m, and source sample size (s_n)
for (rho in c(0, 0.3, 0.6)) {
  for (n in c(100)) {
    for (p_m in c(1000, 2000)) {
      for (s_n in c(NA, 200, 400, 600)) {
        
        # Load result based on method type: HDM or TransHDM
        if (is.na(s_n)) {
          method <- 'simu_HDM'
          load(paste0('simu_HDM_rho_', rho, '_n_', n, '_p_m_', p_m, '.RData'))
        } else {
          method <- 'simu_TransHDM'
          load(paste0('simu_TransHDM_rho_', rho, '_n_', n, '_p_m_', p_m, '_n_s_', s_n, '.RData'))
        }
        
        result <- simu[[method]]
        n_simu <- dim(result)[1]
        true_effect <- simu$true_effect
        
        # True effects
        DE_true <- 1                              # True direct effect
        IDE_true <- sum(true_effect[1:4])         # True total indirect effect (mediators 1–4)
        PE_true <- IDE_true / (DE_true + IDE_true)# True proportion of indirect effect
        
        # Group-specific true IDE and PE for M1–4 (active), M5–6 (β≠0, α=0), M7–8 (α≠0, β=0), M9+ (irrelevant)
        IDE_perM_true <- c(mean(true_effect[1:4]), rep(0, 3))
        PE_perM_true <- IDE_perM_true / (DE_true + IDE_true)
        
        # ------- Evaluate Direct Effect (DE) -------
        DE_est <- result[,"D"]
        DE_bias <- abs(mean(DE_est) - DE_true) / DE_true
        DE_MSE <- sqrt(mean((DE_est - DE_true)^2))
        DE_sd <- sd(DE_est)
        
        # ------- Evaluate Indirect Effect (IDE) -------
        IDE_est <- rowSums(result[,-1])
        IDE_bias <- abs(mean(IDE_est) - IDE_true) / IDE_true
        IDE_MSE <- sqrt(mean((IDE_est - IDE_true)^2))
        IDE_SD <- sd(IDE_est)
        
        # ------- Evaluate Indirect Effect by Mediator Group -------
        IDE_perM_est <- cbind(
          rowMeans(result[, paste0('M', 1:4)]),       # α≠0, β≠0
          rowMeans(result[, paste0('M', 5:6)]),       # α=0, β≠0
          rowMeans(result[, paste0('M', 7:8)]),       # α≠0, β=0
          rowMeans(result[, paste0('M', 9:p_m)])      # α=0, β=0
        )
        IDE_perM_bias <- abs(colMeans(IDE_perM_est) - IDE_perM_true)
        IDE_perM_MSE <- sqrt(colMeans((t(t(IDE_perM_est) - IDE_perM_true)^2)))
        IDE_perM_SD <- apply(IDE_perM_est, 2, sd)
        
        # ------- Evaluate Proportion Effect (PE) -------
        PE_est <- IDE_est / (IDE_est + DE_est)
        PE_bias <- abs(mean(PE_est) - PE_true) / PE_true
        PE_MSE <- sqrt(mean((PE_est - PE_true)^2))
        PE_SD <- sd(PE_est)
        
        # ------- Proportion Effect by Mediator Group -------
        PE_perM_est <- cbind(
          rowMeans(result[, paste0('M', 1:4)]),
          rowMeans(result[, paste0('M', 5:6)]),
          rowMeans(result[, paste0('M', 7:8)]),
          rowMeans(result[, paste0('M', 9:p_m)])
        ) / (IDE_est + DE_est)
        PE_perM_bias <- abs(colMeans(PE_perM_est) - PE_perM_true)
        PE_perM_MSE <- sqrt(colMeans((t(t(PE_perM_est) - PE_perM_true)^2)))
        PE_perM_SD <- apply(PE_perM_est, 2, sd)
        
        # ------- Power (True Positive Rate): M1–M4 -------
        Power <- length(which(result[, paste0('M', 1:4)] > 0)) / length(result[, paste0('M', 1:4)])
        
        # ------- Type I Error (False Positives) by Group -------
        T1E_perM <- c(
          sum(rowSums(result[, paste0('M', 5:6)] > 0)) / n_simu / 2,
          sum(rowSums(result[, paste0('M', 7:8)] > 0)) / n_simu / 2,
          sum(rowSums(result[, paste0('M', 9:p_m)] > 0)) / n_simu / (p_m - 8)
        )
        
        # ------- Family-Wise Error Rate (FWER) -------
        FWER <- length(which(rowSums(result[, paste0('M', 5:p_m)]) != 0)) / n_simu
        
        # ------- Combine All Metrics into One Row -------
        index <- c(rho, n, p_m, s_n, method,
                   DE_bias, DE_MSE, DE_sd,
                   IDE_bias, IDE_MSE, IDE_SD,
                   PE_bias, PE_MSE, PE_SD,
                   Power, FWER,
                   IDE_perM_bias, IDE_perM_MSE, IDE_perM_SD,
                   PE_perM_bias, PE_perM_MSE, PE_perM_SD,
                   T1E_perM)
        
        names(index) <- c('rho', 'n', 'p_m', 's_n', 'method',
                          'DE_bias', 'DE_MSE', 'DE_SD',
                          'IDE_bias', 'IDE_MSE', 'IDE_SD',
                          'PE_bias', 'PE_MSE', 'PE_SD',
                          'Power', 'FWER',
                          paste0('IDE_perM_bias_', 1:4), paste0('IDE_perM_MSE_', 1:4), paste0('IDE_perM_SD_', 1:4),
                          paste0('PE_perM_bias_', 1:4), paste0('PE_perM_MSE_', 1:4), paste0('PE_perM_SD_', 1:4),
                          paste0('T1E_perM_', 2:4))
        
        # Append result to overall summary table
        simu_index <- rbind(simu_index, index)
      }
    }
  }
}
