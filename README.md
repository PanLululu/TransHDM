# **TransHDM**  
### *Transfer Learning for High-Dimensional Mediation Analysis*

---

## 1. Introduction

**TransHDM** is an **R implementation** of a **transfer learning-based high-dimensional mediation analysis** framework. It integrates **transferable source data** to improve the identification of mediators in **small-sample target populations**. This approach is particularly beneficial when **mediator–outcome mechanisms** are **shared across related domains**. Demonstrated improvements in **power** and **robustness** across real and simulated datasets.

TransHDM supports both:

- **Traditional high-dimensional mediation (HDM)** analysis  
- **Transfer learning-enhanced mediation** analysis  

---

##  2. Framework Overview

The TransHDM analysis pipeline consists of the following steps:

1. **Simulate** target and source datasets under **homogeneous** or **heterogeneous** covariate distributions  
2. **Fit** either the traditional HDM or the transfer-enhanced TransHDM models  
3. **Compare** mediator selection results across settings  

---

## 3. Installation

Please ensure the following R packages are installed:

```r
install.packages(c("caret", "MASS", "glmnet", 
                   "parallel", "doParallel", 
                   "foreach", "qvalue"))
```

## 4. Usage
### 4.1. Simulated Data Generation
All necessary R functions are included in the directory. After loading the required functions, you can generate example datasets using:

```r
# For homogeneous covariate design:
source("1_generate_simulationData_homogeneous_design.R")

# For heterogeneous covariate design (to simulate covariate shift):
source("2_generate_simulationData_heterogeneous_design.R")

target_data <- generate_simulationData(n = n, p_x = p_x, p_m = p_m, rho = rho, seed = s)$data
source_data <- generate_simulationData(n = n_s, p_x = p_x, p_m = p_m, rho = rho,
                                       source = TRUE, transferable = TRUE, h = h)$data
```

### 4.2. Fitting Models
```r
# Fit HDM model (no transfer learning):
HDM.fit <- TransHDM(target_data, verbose = TRUE)

# Fit TransHDM model (with transfer learning):
TransHDM.fit <- TransHDM(target_data = target_data, source_data = source_data,
                         transferable = TRUE, verbose = TRUE, ncore = 1)
```

## 5. High-Dimensional Mediation Analysis
We provide a target dataset (sim_target) and a source dataset (sim_source) that simulate the characteristics of real data, used to explore the association between 781 lipid mediators and APOE genotype with tau protein levels in CSF. These datasets are generated to reflect the underlying statistical features of the original data while ensuring privacy protection.

sim_target dataset represents the target group with 37 samples and sim_source dataset represents the source group with 1,059 samples. Both datasets includes the following variables:

- Exposure variable D
- Outcome variable Y
- Three covariates X1, X2, X3
- 781 mediator variables M1 to M781




