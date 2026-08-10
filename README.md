# Communication-efficient distributed hazard difference estimation for heterogeneous multi-site survival data

[![arXiv](https://img.shields.io/badge/arXiv-2601.14609-b31b1b.svg)](https://arxiv.org/abs/2601.14609)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Language](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)

This repository contains the official R implementation for **DiSAH**, a communication-efficient, server-independent federated learning framework designed for absolute risk assessment in privacy-restricted, multi-site time-to-event (survival) data.

---

## 📌 Overview

Traditional privacy-preserving multi-site training frameworks suffer from two fundamental bottlenecks:
1. **Architectural Constraints:** Iterative frameworks (such as FedAvg) require persistent external connections to a central coordinator, which are routinely blocked by rigid hospital firewalls.
2. **Clinical Interpretability:** Most survival models focus on relative effect measures (like hazard ratios from the Cox model), which lack direct clinical interpretability for absolute survival risk assessment.

**DiSAH** addresses these gaps by using an **additive risks model** tailored for right-censored survival data. It enables multi-site absolute risk difference estimation while keeping patient-level data secure and strictly localized.

### 🌟 Key Features
* **Minimal Communication:** Requires only **1 round** of summary statistics exchange for the stratified model and **3 rounds** for the unstratified model.
* **Server-Independent:** Eliminates the need for a persistent central server orchestration—making it highly compatible with isolated hospital networks.
* **Statistical Inference:** Unlike SGD-driven FL frameworks, DiSAH provides valid **confidence intervals** and **hypothesis testing**.
* **Theoretical Guarantees:** Proven to be asymptotically equivalent to pooled individual-level analysis (the "gold standard" where data is centralized).

---

## 📂 Repository Structure

```text
├── src/
│   ├── fedrd_stratified.R     # Core algorithm for the stratified FedRD estimator (1 round)
│   ├── fedrd_unstratified.R   # Core algorithm for the unstratified FedRD estimator (3 rounds)
│   └── utils.R                # Helper functions for data parsing and martingale calculations
├── simulation/
│   ├── generate_data.R        # Script to simulate heterogeneous multi-site survival data
│   └── run_simulations.R      # Main script to reproduce simulation results from the paper
├── data/
│   └── mock_survival_data.csv # Dummy dataset for demo/testing purposes
├── demo.R                     # Walkthrough execution script
└── README.md
```

---

## 🚀 Prerequisites & Installation

DiSAH is written entirely in **R**. Make sure you have R installed (version ≥ 4.0 recommended).

### Required R Packages
Open your R console and run the following command to install the required dependencies:

```R
install.packages(c("survival", "Matrix", "tidyverse"))
```

---

## 💻 Quick Start & Demo

We provide a complete walk-through demo in `demo.R` using simulated mock data to show how both stratified and unstratified models operate.

### Step 1: Clone the Repository
```bash
git clone [https://github.com/siqili0325/DiSAH.git](https://github.com/siqili0325/DiSAH.git)
cd FedRD
```

### Step 2: Run the Demo
You can run the demo directly via your terminal:
```bash
Rscript demo.R
```

### Step 3: Example Usage
Here is a brief look at how the estimators are invoked within an R script:

```R
source("src/fedrd_stratified.R")
source("src/fedrd_unstratified.R")

# 1. Load mock decentralized site data
# (In practice, this step represents summary statistics loaded from separate hospitals)
site_data_list <- readRDS("data/mock_site_data.rds")

# 2. Fit Stratified method (Only 1 round of summary statistics exchange)
stratified_results <- fedrd_stratified(site_data_list)
print(stratified_results$coefficients)
print(stratified_results$confidence_intervals)

# 3. Fit Unstratified method (3 rounds of communication)
unstratified_results <- fedrd_unstratified(site_data_list)
print(unstratified_results$coefficients)
```

---

## 📊 Reproducing Paper Results

To reproduce the simulation tables and figures presented in Section 4 of our paper:

1. Navigate to the `simulation/` directory.
2. Run the benchmarking pipeline to compare against local training and the centralized gold standard:
   ```bash
   Rscript simulation/run_simulations.R
   ```

---

## ✍️ Citation

If you find our framework or code useful for your research, please consider citing our work:

```bibtex
@article{wang2026fedrd,
  title={Communication-efficient distributed hazard difference estimation for heterogeneous multi-site survival data},
  author={Wang, Ziwen and Li, Siqi and Ong, Marcus Eng Hock and Liu, Nan},
  journal={arXiv preprint arXiv:2601.14609},
  year={2026}
}
```

---

## 👥 Contact

* **Ziwen Wang**  <wangziwen@tyut.edu.cn>
* **Siqi Li**  <siqili@u.duke.nus.edu>


For questions, issues, or collaborations regarding the methodology or codebase, please open a GitHub Issue or reach out via email.
