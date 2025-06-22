TransHDM: Transfer Learning for High-Dimensional Mediation Analysis
1. Introduction
TransHDM is an R implementation of a transfer learning based high-dimensional mediation analysis framework that integrates transferable source data to improve the identification of mediators in small-sample target populations. This method is particularly useful when mediator-outcome mechanisms are shared across related domains. TransHDM supports both traditional HDM analysis and transfer learning-enhanced analysis, offering improved power and robustness in real data and simulation studies.

2. TransHDM Framework
The TransHDM pipeline consists of:

Generating simulated target and source datasets under both homogeneous and heterogeneous covariate distributions

Fitting HDM (non-transfer) or TransHDM (transfer) models

Comparing mediator selection results across settings

3. Installation
TransHDM depends on the following R packages. Please make sure they are installed:
install.packages(c("caret", "MASS", "glmnet", "parallel", "doParallel", "foreach", "qvalue"))

4. Usage
4.1. Simulated Data Generation
All necessary R functions are included in the directory. After loading the required functions, you can generate example datasets using:

For homogeneous covariate design:
source("code/1_generate_simulationData_homogeneous_design.R")

For heterogeneous covariate design (to simulate covariate shift):
source("code/2_generate_simulationData_heterogeneous_design.R")

target_data <- generate_simulationData(n = n, p_x = p_x, p_m = p_m, rho = rho, seed = s)$data
source_data <- generate_simulationData(n = n_s, p_x = p_x, p_m = p_m, rho = rho,
                                       source = TRUE, transferable = TRUE, h = h)$data

4.2. Fitting Models
Fit HDM model (no transfer learning):
HDM.fit <- TransHDM(target_data, verbose = TRUE)

Fit TransHDM model (with transfer learning):
TransHDM.fit <- TransHDM(target_data = target_data, source_data = source_data,
                         transferable = TRUE, verbose = TRUE, ncore = 1)
