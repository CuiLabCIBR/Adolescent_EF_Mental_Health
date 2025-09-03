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

This research establishes population-level normative developmental charts for executive function (EF) across adolescence (ages 11-18 years) and investigates how deviations from these norms relate to age-specific mental health symptoms. We utilized a large Chinese adolescent cohort (A-HELP study, N=33,622) and replicated findings in an independent U.S. cohort (Adolescent Brain Cognitive Development (ABCD) study, N=11,549). We assessed inhibitory control (Go/No-Go task) and working memory (1-back, 2-back tasks) and correlated EF deviations with various mental health dimensions measured by the Strengths and Difficulties Questionnaire (SDQ).

## Key Findings

* In a study of over 33,622 Chinese adolescents, all three EF tasks (Go/No-Go, 1-back, 2-back) showed age-related improvements with decreasing inter-individual variability.
* Higher EF deviation scores (indicating better-than-expected performance relative to age norms) were significantly associated with fewer peer problems, conduct problems, lower hyperactivity/inattention, and greater prosocial behavior.
* Age-specific analyses revealed that these associations varied across development, with stronger effect sizes observed in early adolescence (approximately ages 11-13) that declined by late adolescence.
* Findings were replicated in a large U.S. adolescent sample from the ABCD study, showing consistent developmental patterns and EF-mental health associations.

## Repository Structure

This repository is organized into several steps, reflecting the analytical pipeline of the study:

+ `Step2_construct_normative_model/`: Scripts for constructing normative developmental models using GAMLSS.
  + `Step2_1_select_parameter.R`: R script for selecting model parameters (distribution family and degrees of freedom).
  + `Step2_2_bootstrap.R`: R script for bootstrap analysis to determine reliability and stability of the GAMLSS fitted trajectories.
  + `Step2_3_model_evaluate_GNG.Rmd, Step2_3_model_evaluate_1back.Rmd, Step2_3_model_evaluate_2back.Rmd`: R Markdown files for model evaluation, sensitivity analyses, and visualization of normative trajectories for each task.
  + `Step2_3_read_derivative.R`: R script for reading and summarizing the derivatives from the bootstrap analysis.
  + `Step2_4_construct_model.R`: R script for constructing the final normative model and calculating deviation scores.
  + `Step2_5_plot_Normativemodel.R`: R script for plotting the GAMLSS normative curve of A-HELP datasets.

+ `Step3_correlation/`: Scripts for analyzing the correlation between EF deviations and mental health symptoms.
  + `Step3_corr_slope.R`: R script for correlation and slope analysis between EF deviation scores and mental health outcomes.

+ `Step4_interaction/`: Scripts for analyzing age-varying interaction effects.
  + `Step4_EF_psy_int_post_slope_compare.R`: R script for fitting varying-coefficient models to analyze age-specific associations.

+ `Step5_correlation/`: Scripts for analyzing the correlation between EF deviations and mental health symptoms.
  + `Step5_2_1_1_select_parameter_Flanker.R`: R script for selecting model parameters(distribution family) of ABCD dataset.
  + `Step5_2_1_2_select_parameter_Flanker.R`: R script for selecting model parameters(degrees of freedom) of ABCD dataset.
  + `Step5_2_2_bootstrap_Flanker.R`: R script for bootstrap analysis to determine reliability and stability of the GAMLSS fitted trajectories of ABCD dataset.
  + `Step5_2_3_model_evaluate_Flanker.Rmd`: R Markdown files for model evaluation, sensitivity analyses, and visualization of normative trajectories for ABCD flanker task.
  + `Step5_2_4_construct_model_Flanker.R`: R script for constructing the final normative model and calculating deviation scores for flanker task.
  + `Step5_2_5_plot_Normativemodel_Flanker.R`: R script for plotting the GAMLSS normative curve of ABCD datasets.
  + `Step5_3_EF_psychiatry_corr_flanker.R`: R script for correlation and slope analysis between EF deviation scores(flanker task) and mental health outcomes.
  + `Step5_4_EF_psy_int_flanker_compare.R`: R script for fitting varying-coefficient models to analyze age-specific associations with random effects considered.

+ `functions/`: Contains custom R functions used across different analysis steps.
  + `Compare_distributions_gamlss.R`: R functions to compare GAMLSS with different distributions and degrees of freedom for model selection.
  + `Construct_gamlss_set.R`: An R function to construct and return GAMLSS model objects and their performance metrics.
  + `gam_varyingcoefficients.R`: An R function for fitting a GAM with varying coefficients to model interaction effects.
  + `gamcog_withsmoothvar_deviation.R`: An R function to fit and evaluate the linear effects of a continuous variable while controlling for a smooth term.
  + `gamm_varyingcoefficients.R`: An R function for fitting a GAMM with varying coefficients to model interaction effects.
  + `gammcog_withsmoothvar_deviation.R`: An R function to fit and evaluate the linear effects and random effects of a continuous variable while controlling for a smooth term.

## Getting Started

To run the analyses presented in this study, you will need to have R and MATLAB installed. The scripts are designed to be run in a sequential manner, following the `Step` folders.

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/CuiLabCIBR/Adolescent_EF_Mental_Health.git](https://github.com/CuiLabCIBR/Adolescent_EF_Mental_Health.git)
    cd Adolescent_EF_Mental_Health
    ```
2.  **Data Acquisition:** The raw data for the A-HELP study is not publicly available due to data privacy policies. The ABCD data can be accessed via the NDA ([https://nda.nih.gov/abcd](https://nda.nih.gov/abcd)). You will need to process your own data according to the methods described in the manuscript.
3.  **Run Scripts:** Follow the numerical order of the `Step` folders. Most R scripts can be executed directly. Shell scripts (`.sh`) can be run from the terminal.
    ```bash
    # Example for Step 3
    cd Step3_correlation
    Rscript Step3_corr_slope.R
    ```

## Dependencies

* **R (v4.2.2 or later):** With packages such as `gamlss`, `mgcv`, etc. (specific packages are loaded within individual R scripts).

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
