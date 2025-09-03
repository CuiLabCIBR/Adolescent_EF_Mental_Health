rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(mgcv)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)
library(psych)
library(reshape)
library(cowplot)
library(showtext)
wd <- getwd()
interfileFolder <- ' '
FigureFolder  <- ' '

# --- Load Go/No-Go Data ---
sumGNGd_prime_deviation <- readRDS(paste0(interfileFolder,"/Gonogo/GNGd_prime.deviations.rds"))
modGNGd_prime.set1.sum <- readRDS(paste0(interfileFolder, "/Gonogo/GAMLSS.GNGd_primeset1.sum.rds"))
modGNGd_prime.set2.sum <- readRDS(paste0(interfileFolder, "/Gonogo/GAMLSS.GNGd_primeset2.sum.rds"))
modGNGd_prime.all.sum <- readRDS(paste0(interfileFolder, "/Gonogo/GAMLSS.GNGd_prime.all.sum.rds"))
GNGd_prime_data1 <- read_csv(paste0(interfileFolder, "/Gonogo/GNGd_prime.data1.csv"))
GNGd_prime_data2 <- read_csv(paste0(interfileFolder, "/Gonogo/GNGd_prime.data2.csv"))
derivative_GNGd <- read_csv(paste0(interfileFolder, "/Gonogo/derivative_summary_GNGd.csv"))

# --- 1. Plot Age and Sex Distribution ---
fillcolor <- c("#E5E9F2", "#A4C5DF")
gng_plot <- ggplot(data = sumGNGd_prime_deviation, aes(x = Age_year, y = after_stat(count), fill = Sex)) +
  geom_histogram(binwidth = 1, color = "black", boundary = 11, position = "stack", linewidth = 0.5) + 
  labs(x = "Age (year)", y = NULL, title = paste0("Go/No-go, N=", nrow(sumGNGd_prime_deviation))) +
  scale_fill_manual(values = fillcolor, name = "Sex") +
  scale_x_continuous(breaks = seq(11, max(sumGNGd_prime_deviation$Age_year, na.rm = TRUE), by = 1)) + 
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    plot.title = element_text(color = "black", size = 9, hjust = 0.5),
    axis.title = element_text(color = "black", size = 9),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(color = "black", size = 9)
  )
print(gng_plot)
ggsave(filename = paste0(FigureFolder, "/gngnumber.pdf"), plot = gng_plot, width = 9, height=7, units="cm")


# --- 2. Model Performance and Distribution Analysis ---
modGNGd_prime.set1 <- modGNGd_prime.set1.sum$performance.tb
modGNGd_prime.set2 <- modGNGd_prime.set2.sum$performance.tb
print(modGNGd_prime.set1)
print(modGNGd_prime.set2)

# Extract and print performance metrics
partial_Rsq_set1 <- modGNGd_prime.set1$partialRsq
partial_Rsq_set2 <- modGNGd_prime.set2$partialRsq
print(paste("R² for Set1: ", partial_Rsq_set1))
print(paste("R² for Set2: ", partial_Rsq_set2))

# Calculate skewness
skewness_value <- skewness(sumGNGd_prime_deviation$d_prime_deviationZ)
print(paste("Skewness: ", skewness_value))

# Calculate kurtosis
kurtosis_value <- kurtosis(sumGNGd_prime_deviation$d_prime_deviationZ)
print(paste("Kurtosis: ", kurtosis_value))
# Perform Shapiro-Wilk test on a random sample
set.seed(123)
sample_size <- 5000
sample_data <- sample(sumGNGd_prime_deviation$d_prime_deviationZ, size = sample_size, replace = FALSE)
shapiro_test_result <- shapiro.test(sample_data)
print(shapiro_test_result)

# --- 3. Plot Normative Trajectory and Derivative ---
# Calculate centiles from the full model
n_points <- 1000
centiles.all.tmp <- modGNGd_prime.all.sum$centiles_strat  # full model centiles
Centiles <- (centiles.all.tmp[[1]] + centiles.all.tmp[[2]]) / 2
X  <- seq(min(sumGNGd_prime_deviation$Age_year), max(sumGNGd_prime_deviation$Age_year), length.out=1000)
sumGNGd_prime_deviation$Sex <- factor(sumGNGd_prime_deviation$Sex, levels = c("F", "M"))
# Plotting
p_main <- ggplot() +
  geom_hex(data = sumGNGd_prime_deviation, aes(x = Age_year, y = d_prime), bins = 20) +
  scale_fill_gradientn(colors = c("#E5E9F2", "#A4C5DF", "#4980B5")) +
  geom_line(aes(x=X, y=Centiles[3,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[4,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[5,]), linetype="solid", color = "black", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[6,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[7,]), linetype="dashed", color = "gray10", size = 0.5) +
  scale_y_continuous(name="d'", limits = c(-1, 5.5),breaks = seq(-1, 5, by = 2))+
  scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 2))+
  labs(NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(size =9, hjust = 0.5, color = "black"),
    axis.text = element_text(colour = "black",size=9), 
    axis.title = element_text(size =9,face = "plain"),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.ticks = element_line(colour = "black", size = 0.25), 
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.y = element_text(colour = "black",size =9),  
    axis.line.y = element_line(colour = "black", size = 0.25),
    plot.background=element_rect(fill="white",color = NA),
    legend.position = "none" ,
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
x_length <- 7  
y_length <- 6
print(p_main)
ggsave(
  filename = paste0(paste0(FigureFolder, "/NormativeDevCurve_GNGd_prime_blue.pdf")),
  plot = p_main,
  width = x_length,  
  height = y_length, 
  units = "cm"
)

# Bar plot for derivative significance
derivative_GNGd <- derivative_GNGd %>% drop_na(P50_lower, P50_upper)
derivative_GNGd <- derivative_GNGd %>%
  mutate(significance = ifelse(P50_lower > 0, "Increasing",
                               ifelse(P50_upper < 0, "Decreasing", "Non-significant")))
derivative_GNGd$color_group <- ifelse(derivative_GNGd$significance == "Non-significant", NA, derivative_GNGd$P50_mean)
p_bar <- ggplot(derivative_GNGd) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="#A4C5DF", low="white", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="#A4C5DF", low="white", midpoint=0, na.value="white", labels=NULL) +
  annotate("rect", xmin=min(derivative_GNGd$Age)-0.025, xmax=max(derivative_GNGd$Age)+0.025, ymin=-0.025, ymax=0.525, 
           color="black", fill=NA, size=0.25) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab("Age (years)") +
  scale_x_continuous(breaks=NULL) +
  labs(fill="mu Mean") +
  theme_void() +
  theme(
    axis.text=element_text(size=9, color='black'),
    legend.position="none",
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    plot.margin=unit(c(0, 0, 0, 0), "cm")
  )
# Combine and save normative curve plot
x_length <- 7
y_length_main <- 6
y_length_bar <- 0.3 
final_plot <- p_main / p_bar + plot_layout(heights = c(y_length_main, y_length_bar))
print(final_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_GNGd_prime_Addbar_blue.pdf"), 
  plot = final_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)

# --- 4. Plot Variance Trajectory and Derivative ---
lower_CI <- derivative_GNGd$sigma_pred_lower
upper_CI <- derivative_GNGd$sigma_pred_upper
p_sigma <- ggplot(derivative_GNGd, aes(x = Age, y = sigma_pred_mean)) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "gray80") +
  #geom_hline(yintercept = seq(0.8, 1, by = 0.05), linetype = "solid", color = "gray90", size = 0.5) + 
  geom_line(size = 0.5, color = "black") +
  labs(NULL) +
  theme_minimal() +
  scale_y_continuous(name="Variance", limits = c(0.65, 1.05), breaks = seq(0.7, 1, by = 0.1)) +
  scale_x_continuous(name="", limits = c(11, 18), breaks = seq(11, 18, by = 2)) +
  theme(
    plot.title = element_text(size =9, hjust = 0.5),
    axis.text = element_text(colour = "black",size =9),
    axis.title = element_text(size =9, hjust = 0.5),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.ticks = element_line(colour = "black", size = 0.25), 
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.y = element_text(colour = "black",size =9),  
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
x_length <- 7
y_length <- 6
p_sigma
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_GNGd_prime.pdf"),
  plot = p_sigma,
  width = x_length, 
  height = y_length,  
  units = "cm"
)
# Bar plot for variance derivative significance
derivative_GNGd <- derivative_GNGd %>% drop_na(sigma_deriv_lower, sigma_deriv_upper)
derivative_GNGd <- derivative_GNGd %>%
  mutate(significance = ifelse(sigma_pred_lower > 0, "Increasing",
                               ifelse(sigma_pred_upper < 0, "Decreasing", "Non-significant")))
derivative_GNGd$color_group <- ifelse(derivative_GNGd$significance == "Non-significant", NA, derivative_GNGd$sigma_deriv_mean)
p_sigmabar <- ggplot(derivative_GNGd) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="white", low="#A4C5DF", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="white", low="#A4C5DF", midpoint=0, na.value="white", labels=NULL) +
  annotate("rect", xmin=min(derivative_GNGd$Age)-0.025, xmax=max(derivative_GNGd$Age)+0.025, ymin=-0.025, ymax=0.525, 
           color="black", fill=NA, size=0.25) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab("Age (years)") +
  scale_x_continuous(breaks=NULL) +
  labs(fill="mu Mean") +
  theme_void() +
  theme(
    axis.text=element_text(size=9, color='black'),
    legend.position="none",
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    plot.margin=unit(c(0, 0, 0, 0), "cm")
  )
x_length <- 7 
y_length_main <- 6
y_length_bar <- 0.3
sigma_plot <- p_sigma / p_sigmabar + plot_layout(heights = c(y_length_main, y_length_bar))
print(sigma_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_GNGd_prime_Addbar.pdf"), 
  plot = sigma_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)


# --- Load 1-back Data ---
sum1backAcc_deviation <- readRDS(paste0(interfileFolder,"/1-back/back1Acc.deviations.rds"))
mod1backAcc.set1.sum <- readRDS(paste0(interfileFolder, "/1-back/GAMLSS.back1Accset1.sum.rds"))
mod1backAcc.set2.sum <- readRDS(paste0(interfileFolder, "/1-back/GAMLSS.back1Accset2.sum.rds"))
mod1backAcc.all.sum <- readRDS(paste0(interfileFolder, "/1-back/GAMLSS.back1Acc.all.sum.rds"))
back1Acc_data1 <- read_csv(paste0(interfileFolder, "/1-back/back1Acc.data1.csv"))
back1Acc_data2 <- read_csv(paste0(interfileFolder, "/1-back/back1Acc.data2.csv"))
derivative_back1 <- read_csv(paste0(interfileFolder, "/1-back/derivative_summary_back1.csv"))

# --- 1. Plot Age and Sex Distribution ---
fillcolor <- c("#E5E9F2", "#A4C5DF")
back1_plot <- ggplot(data = sum1backAcc_deviation, aes(x = Age_year, y = after_stat(count), fill = Sex)) +
  geom_histogram(binwidth = 1, color = "black", boundary = 11, position = "stack", linewidth = 0.5) + 
  labs(x = "Age (year)", y = NULL, title = paste0("1-back, N=", nrow(sum1backAcc_deviation))) +
  scale_fill_manual(values = fillcolor, name = "Sex") +
  scale_x_continuous(breaks = seq(11, max(sum1backAcc_deviation$Age_year, na.rm = TRUE), by = 1)) + 
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    plot.title = element_text(color = "black", size = 9, hjust = 0.5),
    axis.title = element_text(color = "black", size = 9),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(color = "black", size = 9),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )
print(back1_plot)
ggsave(filename = paste0(FigureFolder, "/onebacknumber.pdf"), plot = back1_plot, width = 9, height=7, units="cm")

# --- 2. Model Performance and Distribution Analysis ---
mod1backAcc.set1 <- mod1backAcc.set1.sum$performance.tb
mod1backAcc.set2 <- mod1backAcc.set2.sum$performance.tb
print(mod1backAcc.set1)
print(mod1backAcc.set2)

# Extract and print performance metrics
partial_Rsq_set1 <- mod1backAcc.set1$partialRsq
partial_Rsq_set2 <- mod1backAcc.set2$partialRsq
print(paste("R² for Set1: ", partial_Rsq_set1))
print(paste("R² for Set2: ", partial_Rsq_set2))

# Calculate skewness 
skewness_value <- skewness(sum1backAcc_deviation$Oneback_acc_deviationZ)
print(paste("Skewness: ", skewness_value))

# Calculate kurtosis
kurtosis_value <- kurtosis(sum1backAcc_deviation$Oneback_acc_deviationZ)
print(paste("Kurtosis: ", kurtosis_value))

# Perform Shapiro-Wilk test on a random sample
set.seed(123)
sample_size <- 5000
sample_data <- sample(sum1backAcc_deviation$Oneback_acc_deviationZ, size = sample_size, replace = FALSE)
shapiro_test_result <- shapiro.test(sample_data)
print(shapiro_test_result)

# --- 3. Plot Normative Trajectory and Derivative ---
# Calculate centiles from the full model
n_points <- 1000
centiles.all.tmp <- mod1backAcc.all.sum$centiles_strat  # full model centiles

Centiles <- (centiles.all.tmp[[1]] + centiles.all.tmp[[2]]) / 2
X  <- seq(min(sum1backAcc_deviation$Age_year), max(sum1backAcc_deviation$Age_year), length.out=1000)

# Main plot for the normative curve
p_main <- ggplot() +
  geom_hex(data = sum1backAcc_deviation, aes(x = Age_year, y = Oneback_acc), bins = 20) +
  scale_fill_gradientn(colors = c("#E5E9F2", "#A4C5DF", "#4980B5")) +
  geom_line(aes(x=X, y=Centiles[3,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[4,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[5,]), linetype="solid", color = "black", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[6,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[7,]), linetype="dashed", color = "gray10", size = 0.5) +
  scale_y_continuous(name="Accuracy",limits = c(0.1, 1),breaks = seq(0.2, 1, by = 0.2))+
  scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 2))+
  labs(NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5, color = "black"),
    axis.text = element_text(colour = "black",size=9), 
    axis.title = element_text(size=9,face = "plain"),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.ticks = element_line(colour = "black", size = 0.25), 
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.y = element_text(colour = "black",size = 9),
    plot.background=element_rect(fill="white",color = NA),
    legend.position = "none" ,
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
x_length <- 7 
y_length <- 6
print(p_main)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_1backAcc_blue.pdf"),
  plot = p_main,
  width = x_length,  
  height = y_length, 
  units = "cm"
)
# Bar plot for derivative significance
derivative_back1 <- derivative_back1 %>% drop_na(P50_lower, P50_upper)
derivative_back1 <- derivative_back1 %>%
  mutate(significance = ifelse(P50_lower > 0, "Increasing",
                               ifelse(P50_upper < 0, "Decreasing", "Non-significant")))
derivative_back1$color_group <- ifelse(derivative_back1$significance == "Non-significant", NA, derivative_back1$P50_mean)
p_bar <- ggplot(derivative_back1) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="#A4C5DF", low="white", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="#A4C5DF", low="white", midpoint=0, na.value="white", labels=NULL) +
  annotate("rect", xmin=min(derivative_back1$Age)-0.025, xmax=max(derivative_back1$Age)+0.025, ymin=-0.025, ymax=0.525, 
           color="black", fill=NA, size=0.25) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab("Age (years)") +
  scale_x_continuous(breaks=NULL) +
  labs(fill="mu Mean") +
  theme_void() +
  theme(
    axis.text=element_text(size=9, color='black'),
    legend.position="none",
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    plot.margin=unit(c(0, 0, 0, 0), "cm")
  )
# Combine and save normative curve plot
x_length <- 7
y_length_main <- 6
y_length_bar <- 0.3 
final_plot <- p_main / p_bar + plot_layout(heights = c(y_length_main, y_length_bar))
print(final_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_1backAcc_Addbar_blue.pdf"), 
  plot = final_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)

# --- 4. Plot Variance Trajectory and Derivative ---
# Main plot for variance
lower_CI <- derivative_back1$sigma_pred_lower
upper_CI <- derivative_back1$sigma_pred_upper
p_sigma <- ggplot(derivative_back1, aes(x = Age, y = sigma_pred_mean)) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "gray80") +
  geom_line(size = 0.5, color = "black") +
  labs(NULL) +
  theme_minimal() +
  scale_y_continuous(name="Variance", limits = c(0.1, 0.28), breaks = seq(0.1, 0.27, by = 0.05)) +
  scale_x_continuous(name="", limits = c(11, 18), breaks = seq(11, 18, by = 2)) +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5),
    axis.text = element_text(colour = "black",size = 9),
    axis.title = element_text(size = 9, hjust = 0.5),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.ticks = element_line(colour = "black", size = 0.25), 
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.y = element_text(colour = "black",size = 9),   
    panel.grid = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"))
x_length <- 7
y_length <- 6
p_sigma
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_back1Acc.pdf"),
  plot = p_sigma,
  width = x_length, 
  height = y_length,  
  units = "cm"
)
# Bar plot for variance derivative significance
derivative_back1 <- derivative_back1 %>% drop_na(sigma_deriv_lower, sigma_deriv_upper)
derivative_back1 <- derivative_back1 %>%
  mutate(significance = ifelse(sigma_deriv_lower > 0, "Increasing",
                               ifelse(sigma_deriv_upper < 0, "Decreasing", "Non-significant")))

derivative_back1$color_group <- ifelse(derivative_back1$significance == "Non-significant", NA, derivative_back1$sigma_deriv_mean)

p_sigmabar <- ggplot(derivative_back1) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="white", low="#A4C5DF", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="white", low="#A4C5DF", midpoint=0, na.value="white", labels=NULL) +
  annotate("rect", xmin=min(derivative_back1$Age)-0.025, xmax=max(derivative_back1$Age)+0.025, ymin=-0.025, ymax=0.525, 
           color="black", fill=NA, size=0.25) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks=NULL) +
  labs(fill="P50 Mean") +
  theme_void() +
  theme(
    axis.text=element_text(size=9, color='black'),
    legend.position="none",
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    plot.margin=unit(c(0,0,0,0), "cm")
  )
# Combine and save variance plot
x_length <- 7
y_length_main <- 6
y_length_bar <- 0.3
sigma_plot <- p_sigma / p_sigmabar + plot_layout(heights = c(y_length_main, y_length_bar))
print(sigma_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_1backAcc_Addbar.pdf"), 
  plot = sigma_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)




# --- Load 2-back Data ---
sum2backAcc_deviation <- readRDS(paste0(interfileFolder,"/2-back/back2Acc.deviations.rds"))
mod2backAcc.set1.sum <- readRDS(paste0(interfileFolder, "/2-back/GAMLSS.back2Accset1.sum.rds"))
mod2backAcc.set2.sum <- readRDS(paste0(interfileFolder, "/2-back/GAMLSS.back2Accset2.sum.rds"))
mod2backAcc.all.sum <- readRDS(paste0(interfileFolder, "/2-back/GAMLSS.back2Acc.all.sum.rds"))
back2Acc_data1 <- read_csv(paste0(interfileFolder, "/2-back/back2Acc.data1.csv"))
back2Acc_data2 <- read_csv(paste0(interfileFolder, "/2-back/back2Acc.data2.csv"))
derivative_back2 <- read_csv(paste0(interfileFolder, "/2-back/derivative_summary_back2.csv"))

# --- 1. Plot Age and Sex Distribution ---
# Create a histogram to visualize the sample demographics
fillcolor <- c("#E5E9F2", "#A4C5DF")
back2_plot <- ggplot(data = sum2backAcc_deviation, aes(x = Age_year, y = after_stat(count), fill = Sex)) +
  geom_histogram(binwidth = 1, color = "black", boundary = 11, position = "stack", linewidth = 0.5) + 
  labs(x = "Age (year)", y = NULL, title = paste0("2-back, N=", nrow(sum2backAcc_deviation))) +
  scale_fill_manual(values = fillcolor, name = "Sex") +
  scale_x_continuous(breaks = seq(11, max(sum2backAcc_deviation$Age_year, na.rm = TRUE), by = 1)) + 
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    plot.title = element_text(color = "black", size = 8.5, hjust = 0.5),
    axis.title = element_text(color = "black", size = 8.5),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(color = "black", size = 8.5),
    legend.text = element_text(size = 8.5),
    legend.title = element_text(size = 8.5)
  )
print(back2_plot)
ggsave(filename = paste0(FigureFolder, "/twobacknumber.pdf"), plot = back2_plot, width = 9, height=7, units="cm")


# --- 2. Model Performance and Distribution Analysis ---
mod2backAcc.set1 <- mod2backAcc.set1.sum$performance.tb
mod2backAcc.set2 <- mod2backAcc.set2.sum$performance.tb
print(mod2backAcc.set1)
print(mod2backAcc.set2)

# Extract and print performance metrics
partial_Rsq_set1 <- mod2backAcc.set1$partialRsq
partial_Rsq_set2 <- mod2backAcc.set2$partialRsq
print(paste("R² for Set1: ", partial_Rsq_set1))
print(paste("R² for Set2: ", partial_Rsq_set2))

# Calculate skewness
skewness_value <- skewness(sum2backAcc_deviation$Twoback_acc_deviationZ)
print(paste("Skewness: ", skewness_value))

# Calculate kurtosis
kurtosis_value <- kurtosis(sum2backAcc_deviation$Twoback_acc_deviationZ)
print(paste("Kurtosis: ", kurtosis_value))
# Perform Shapiro-Wilk test on a random sample
set.seed(123)
sample_size <- 5000
sample_data <- sample(sum2backAcc_deviation$Twoback_acc_deviationZ, size = sample_size, replace = FALSE)
shapiro_test_result <- shapiro.test(sample_data)
print(shapiro_test_result)


# --- 3. Plot Normative Trajectory and Derivative ---
# Calculate centiles from the full model
n_points <- 1000
centiles.all.tmp <- mod2backAcc.all.sum$centiles_strat  # full model centiles

Centiles <- (centiles.all.tmp[[1]] + centiles.all.tmp[[2]]) / 2
X  <- seq(min(sum2backAcc_deviation$Age_year), max(sum2backAcc_deviation$Age_year), length.out=1000)
# Main plot for the normative curve
p_main <- ggplot() +
  geom_hex(data = sum2backAcc_deviation, aes(x = Age_year, y = Twoback_acc), bins = 20) +
  scale_fill_gradientn(colors = c("#E5E9F2", "#A4C5DF", "#4980B5")) +
  geom_line(aes(x=X, y=Centiles[3,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[4,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[5,]), linetype="solid", color = "black", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[6,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[7,]), linetype="dashed", color = "gray10", size = 0.5) +
  scale_y_continuous(name="Accuracy",limits = c(0.1, 1),breaks = seq(0.2, 1, by = 0.2))+
  scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 2))+
  labs(NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5, color = "black"),
    axis.text = element_text(colour = "black",size=9), 
    axis.title = element_text(size=9,face = "plain"),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.ticks = element_line(colour = "black", size = 0.25), 
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.y = element_text(colour = "black",size = 9),
    plot.background=element_rect(fill="white",color = NA),
    legend.position = "none" ,
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
x_length <- 7  
y_length <- 6
print(p_main)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_2backAcc_blue.pdf"),
  plot = p_main,
  width = x_length,
  height = y_length, 
  units = "cm"
)

# Bar plot for derivative significance
derivative_back2 <- derivative_back2 %>% drop_na(P50_lower, P50_upper)
derivative_back2 <- derivative_back2 %>%
  mutate(significance = ifelse(P50_lower > 0, "Increasing",
                               ifelse(P50_upper < 0, "Decreasing", "Non-significant")))
derivative_back2$color_group <- ifelse(derivative_back2$significance == "Non-significant", NA, derivative_back2$P50_mean)

p_bar <- ggplot(derivative_back2) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="#A4C5DF", low="white", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="#A4C5DF", low="white", midpoint=0, na.value="white", labels=NULL) +
  # 使用 annotate() 绘制矩形边框
  annotate("rect", xmin=min(derivative_back2$Age)-0.025, xmax=max(derivative_back2$Age)+0.025, ymin=-0.025, ymax=0.525, 
           color="black", fill=NA, size=0.25) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks=NULL) +
  labs(fill="mu Mean") +
  theme_void() +
  theme(
    axis.text=element_text(size=9, color='black'),
    legend.position="none",
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    plot.margin=unit(c(0, 0, 0, 0), "cm")
  )
# Combine and save normative curve plot
x_length <- 7
y_length_main <- 6
y_length_bar <- 0.3
final_plot <- p_main / p_bar + plot_layout(heights = c(y_length_main, y_length_bar))
print(final_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_2backAcc_Addbar_blue.pdf"), 
  plot = final_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)

# --- 4. Plot Variance Trajectory and Derivative ---
# Main plot for variance
lower_CI <- derivative_back2$sigma_pred_lower
upper_CI <- derivative_back2$sigma_pred_upper
p_sigma <- ggplot(derivative_back2, aes(x = Age, y = sigma_pred_mean)) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "gray80") +
  geom_line(size = 0.5, color = "black") +
  labs(NULL) +
  theme_minimal() +
  scale_y_continuous(name="Variance", limits = c(0.13, 0.2), breaks = seq(0.12, 0.2, by = 0.02)) +
  scale_x_continuous(name="", limits = c(11, 18), breaks = seq(11, 18, by = 2)) +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5),
    axis.text = element_text(colour = "black",size = 9),
    axis.title = element_text(size = 9, hjust = 0.5),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.ticks = element_line(colour = "black", size = 0.25), 
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.y = element_text(colour = "black",size = 9),   
    panel.grid = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"))
x_length <- 7
y_length <- 6
p_sigma
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_back2Acc.pdf"),
  plot = p_sigma,
  width = x_length, 
  height = y_length,  
  units = "cm"
)
# Bar plot for variance derivative significance
derivative_back2 <- derivative_back2 %>% drop_na(sigma_deriv_lower, sigma_deriv_upper)
derivative_back2 <- derivative_back2 %>%
  mutate(significance = ifelse(sigma_deriv_lower > 0, "Increasing",
                               ifelse(sigma_deriv_upper < 0, "Decreasing", "Non-significant")))

derivative_back2$color_group <- ifelse(derivative_back2$significance == "Non-significant", NA, derivative_back2$sigma_deriv_mean)

p_sigmabar <- ggplot(derivative_back2) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="white", low="#A4C5DF", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="white", low="#A4C5DF", midpoint=0, na.value="white", labels=NULL) +
  annotate("rect", xmin=min(derivative_back2$Age)-0.025, xmax=max(derivative_back2$Age)+0.025, ymin=-0.025, ymax=0.525, 
           color="black", fill=NA, size=0.25) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks=NULL) +
  labs(fill="P50 Mean") +
  theme_void() +
  theme(
    axis.text=element_text(size=9, color='black'),
    legend.position="none",
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    plot.margin=unit(c(0,0,0,0), "cm")
  )
# Combine and save variance plot
x_length <- 7
y_length_main <- 6
y_length_bar <- 0.3
sigma_plot <- p_sigma / p_sigmabar + plot_layout(heights = c(y_length_main, y_length_bar))
print(sigma_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_2backAcc_Addbar.pdf"), 
  plot = sigma_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)
