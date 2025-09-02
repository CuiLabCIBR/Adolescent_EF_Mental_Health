rm(list = ls())
library(tidyverse)
library(ggplot2)
wd <- getwd()
if (str_detect(wd,"cuizaixu_lab")) {
  interfileFolder <- " "
} else {
  interfileFolder <- " "   
}

## ========= Select TASK (uncomment one at a time) =========
# # ---- Go/No-Go ----
# TASK <- list(
#   out_subdir = "GNGd_prime",  # Subfolder containing bootstrap .rds files
#   label      = "Go/No-Go",    # Used for figure titles
#   file_tag   = "GNGd"         # Used for file names
# )

# # ---- 1-back ----
# TASK <- list(
#   out_subdir = "1-back",
#   label      = "1-back",
#   file_tag   = "back1"
# )

# # ---- 2-back ----
# TASK <- list(
#   out_subdir = "2-back",
#   label      = "2-back",
#   file_tag   = "back2"
# )

# ---- Flanker ----
TASK <- list(
  out_subdir = "Flanker",
  label      = "Flanker",
  file_tag   = "Flanker"
)

## ========= Directories & file list =========
bootstrap_dir <- file.path(interfileFolder, "bootstrap", TASK$out_subdir)
output_dir    <- file.path(interfileFolder, "bootstrap")   # Save all outputs here
if (!dir.exists(bootstrap_dir)) stop("Directory not found: ", bootstrap_dir)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

file_list <- list.files(path = bootstrap_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(file_list) == 0) stop("No .rds files found in: ", bootstrap_dir)

## ========= Read all bootstrap results =========
mu_derivatives_all    <- vector("list", length(file_list))  # For P50 derivatives
sigma_derivatives_all <- vector("list", length(file_list))  # For sigma derivatives
sigma_all             <- vector("list", length(file_list))  # For sigma predictions
age_points            <- NULL

for (i in seq_along(file_list)) {
  x <- readRDS(file_list[i])
  
  # Extract age grid once
  if (is.null(age_points)) age_points <- x$Age_points
  
  # Extract derivatives and predictions
  mu_derivatives_all[[i]]    <- x$P50_derivative
  sigma_derivatives_all[[i]] <- if (!is.null(x$sigma_derivative_all)) x$sigma_derivative_all else x$sigma_derivative
  sigma_all[[i]]             <- x$sigma_pred
}

# Convert lists to matrices
mu_derivatives_matrix    <- do.call(cbind, mu_derivatives_all)
sigma_derivatives_matrix <- do.call(cbind, sigma_derivatives_all)
sigma_matrix             <- do.call(cbind, sigma_all)

## ========= Compute mean & 95% CI across bootstraps =========
P50_mean <- rowMeans(mu_derivatives_matrix, na.rm = TRUE)
P50_ci   <- apply(mu_derivatives_matrix, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

sigma_deriv_mean <- rowMeans(sigma_derivatives_matrix, na.rm = TRUE)
sigma_deriv_ci   <- apply(sigma_derivatives_matrix, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

sigma_pred_mean <- rowMeans(sigma_matrix, na.rm = TRUE)
sigma_pred_ci   <- apply(sigma_matrix, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

# Combine all results into one dataframe
derivative_summary <- data.frame(
  Age               = age_points,
  P50_mean          = P50_mean,
  P50_lower         = P50_ci[1, ],
  P50_upper         = P50_ci[2, ],
  sigma_deriv_mean  = sigma_deriv_mean,
  sigma_deriv_lower = sigma_deriv_ci[1, ],
  sigma_deriv_upper = sigma_deriv_ci[2, ],
  sigma_pred_mean   = sigma_pred_mean,
  sigma_pred_lower  = sigma_pred_ci[1, ],
  sigma_pred_upper  = sigma_pred_ci[2, ]
)

## ========= Identify significance where 95% CI excludes 0 =========
derivative_summary <- derivative_summary %>%
  mutate(
    P50_sig         = (P50_lower > 0) | (P50_upper < 0),
    sigma_deriv_sig = (sigma_deriv_lower > 0) | (sigma_deriv_upper < 0)
  )

## ========= Plot results =========
p1 <- ggplot(derivative_summary, aes(x = Age, y = P50_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = P50_lower, ymax = P50_upper), alpha = 0.2) +
  labs(title = paste0(TASK$label, ": d(μ)/dAge (Mean & 95% CI)"),
       x = "Age", y = "d(μ)/dAge") +
  theme_minimal()

p2 <- ggplot(derivative_summary, aes(x = Age, y = sigma_deriv_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = sigma_deriv_lower, ymax = sigma_deriv_upper), alpha = 0.2) +
  labs(title = paste0(TASK$label, ": d(σ)/dAge (Mean & 95% CI)"),
       x = "Age", y = "d(σ)/dAge") +
  theme_minimal()

p3 <- ggplot(derivative_summary, aes(x = Age, y = sigma_pred_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = sigma_pred_lower, ymax = sigma_pred_upper), alpha = 0.2) +
  labs(title = paste0(TASK$label, ": σ (Mean & 95% CI)"),
       x = "Age", y = "σ") +
  theme_minimal()

## ========= Save outputs =========
# Save summary table
write.csv(
  derivative_summary,
  file = file.path(output_dir, paste0("derivative_summary_", TASK$file_tag, ".csv")),
  row.names = FALSE
)

# Save plots
ggsave(file.path(output_dir, paste0("mu_derivative_age_", TASK$file_tag, ".png")),
       plot = p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, paste0("sigma_derivative_age_", TASK$file_tag, ".png")),
       plot = p2, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, paste0("sigma_pred_age_", TASK$file_tag, ".png")),
       plot = p3, width = 8, height = 6, dpi = 300)

cat("Done for task:", TASK$label, "\n",
    "Input dir:", bootstrap_dir, "\n",
    "Output dir:", output_dir, "\n")