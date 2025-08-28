# Adolescent_EF_Mental_Health

This repository contains the code and data analysis scripts associated with the research paper:
"Normative Executive Function Development Reveals Age-Varying Mental Health Associations in Youth" 

## Table of Contents

- [About This Study](#about-this-study)
- [Key Findings](#key-findings)
- [Repository Structure](#repository-structure)
- [Getting Started](#getting-started)
- [Dependencies](#dependencies)
- [Authors](#authors)

## About This Study

This research establishes population-level normative developmental charts for executive function (EF) across adolescence (ages 11-18 years) and investigates how deviations from these norms relate to age-specific mental health symptoms. We utilized a large Chinese adolescent cohort (A-HELP study, N=33,650) and replicated findings in an independent U.S. cohort (Adolescent Brain Cognitive Development (ABCD) study, N=11,715). We assessed inhibitory control (Go/No-Go task) and working memory (1-back, 2-back tasks) and correlated EF deviations with various mental health dimensions measured by the Strengths and Difficulties Questionnaire (SDQ).

## Key Findings

* In a study of over 33,650 Chinese adolescents, all three EF tasks (Go/No-Go, 1-back, 2-back) showed age-related improvements with decreasing inter-individual variability.
* Higher EF deviation scores (indicating better-than-expected performance relative to age norms) were significantly associated with fewer peer problems, conduct problems, lower hyperactivity/inattention, and greater prosocial behavior.
* Age-specific analyses revealed that these associations varied across development, with stronger effect sizes observed in early adolescence (approximately ages 11-13) that declined by late adolescence.
* Findings were replicated in a large U.S. adolescent sample from the ABCD study, showing consistent developmental patterns and EF-mental health associations.

## Repository Structure

This repository is organized into several steps, reflecting the analytical pipeline of the study:

* `Step1_clean_data/`: Contains scripts for cleaning and preprocessing raw executive function data (Go/No-Go, 1-back, 2-back).
  * `GNG_clean_yf.mlx`: MATLAB scripts for Go/No-Go data cleaning.
  * `oneback_clean_yf.mlx`, `twoback_clean_yf.mlx`: MATLAB scripts for 1-back and 2-back data cleaning.
  * `check_demo.R`: R script for checking demographic data.
* `Step2_construct_normative_model/`: Scripts for constructing normative developmental models using GAMLSS.
  * `Step2_1_select_parameter.R`: R script for selecting model parameters.
  * `Step2_2_bootstrap.R`: R script for bootstrap analysis.
  * `Step2_2_1_read_derivative.R`: R script for reading derivatives.
  * `Step2_3_construct_model.R`: R script for constructing the normative model.
* `Step3_correlation/`: Scripts for analyzing the correlation between EF deviations and mental health symptoms.
  * `Run_step3.sh`: Shell script to execute Step 3 analyses.
  * `Step3_corr_slope.R`, `Step3_corr_slope10000.R`: R scripts for correlation and slope analysis.
* `Step4_interaction/`: Scripts for analyzing age-varying interaction effects.
  * `Run_step4.sh`: Shell script to execute Step 4 analyses.
  * `Step4_EF_psy_int_post_slope_compare.R`: R script for comparing interaction effects.
* `Step5_ABCD/`: Scripts specifically for the replication analysis using the ABCD dataset.
  * `EF_psy_int_flanker_compare.R`: R script for comparing EF and mental health interactions with Flanker task.
  * `Step3_EF_psychiatry_corr_flanker.R`: R script for correlation analysis with Flanker task.
* `functions/`: Contains custom R functions used across different analysis steps.
  * `gamm_factor_interaction_deviation.R`: R function for GAMM factor interaction deviation.
  * `gamm_varyingcoefficients_new.R`: R function for GAMM with varying coefficients.

## Getting Started

To run the analyses presented in this study, you will need to have R and MATLAB installed. The scripts are designed to be run in a sequential manner, following the `Step` folders.

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/CuiLabCIBR/Adolescent_EF_Mental_Health.git](https://github.com/CuiLabCIBR/Adolescent_EF_Mental_Health.git)
    cd Adolescent_EF_Mental_Health
    ```
2.  **Data Acquisition:** The raw data for the A-HELP study is not publicly available due to data privacy policies. The ABCD data can be accessed via the NDA ([https://nda.nih.gov/abcd](https://nda.nih.gov/abcd)). You will need to process your own data according to the methods described in the manuscript and the cleaning scripts in `Step1_clean_data/`.
3.  **Run Scripts:** Follow the numerical order of the `Step` folders. Most R scripts can be executed directly. Shell scripts (`.sh`) can be run from the terminal.
    ```bash
    # Example for Step 3
    cd Step3_correlation
    ./Run_step3.sh
    ```

## Dependencies

* **R (v4.2.2 or later):** With packages such as `gamlss`, `mgcv`, etc. (specific packages are loaded within individual R scripts).
* **MATLAB:** Required for running the `.mlx` cleaning scripts.

## Authors

* Yang Li
* Lirou Tan
* Yinan Duan
* Xiaoyu Xu
* Haoshu Xu
* Mei Yu
* Luxia Jia
* Zhilin Li
* Chenguang Zhao
* Qunlin Chen
* Bart Larsen
* Adam Pines
* Tengfei Wang
* Runsen Chen (Correspondence)
* Zaixu Cui (Correspondence)
