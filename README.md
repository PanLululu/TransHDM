# **TransHDM**  
### *Transfer Learning for High-Dimensional Mediation Analysis*

---

## Introduction

**TransHDM** is an **R implementation** of a **transfer learning-based high-dimensional mediation analysis** framework. It integrates **transferable source data** to improve the identification of mediators in **small-sample target populations**. This approach is particularly beneficial when **mediator–outcome mechanisms** are **shared across related domains**. Demonstrated improvements in **power** and **robustness** across real and simulated datasets.

TransHDM supports both:

- **Traditional high-dimensional mediation (HDM)** analysis  
- **Transfer learning-enhanced mediation** analysis  

---

##  Framework Overview

The TransHDM analysis pipeline consists of the following steps:

1. **Simulate** target and source datasets under **homogeneous** or **heterogeneous** covariate distributions  
2. **Fit** either the traditional HDM or the transfer-enhanced TransHDM models  
3. **Compare** mediator selection results across settings  

---

## Installation

Please ensure the following R packages are installed:

```r
install.packages(c("caret", "MASS", "glmnet", 
                   "parallel", "doParallel", 
                   "foreach", "qvalue"))
```

## Usage
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
