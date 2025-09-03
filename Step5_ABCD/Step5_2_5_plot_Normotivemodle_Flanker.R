rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(mgcv)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)
library(psych)
library(reshape)
library(cowplot)
library(writexl)
library(extrafont)
library(patchwork)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- ' '
  demopath <- ' '
  interfileFolder <-  ' '
  functionFolder <- ' '
  resultFolder <- ' '
}else{
  FigureFolder <- ' '
  interfileFolder <- ' '
  functionFolder <- ' '
}
sumFlanker_deviation <- readRDS(paste0(interfileFolder,"/Flanker.deviations.rds"))
modFlanker.set1.sum <- readRDS(paste0(interfileFolder, "/GAMLSS.Flankerset1.sum.rds"))
modFlanker.set2.sum <- readRDS(paste0(interfileFolder, "/GAMLSS.Flankerset2.sum.rds"))
modFlanker.all.sum <- readRDS("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/data_2507/GAMLSS.Flanker.all.sum.rds")
Flanker_data1 <- read_csv(paste0(interfileFolder, "/Flanker_data.set1.csv"))
Flanker_data2 <- read_csv(paste0(interfileFolder, "/Flanker_data.set2.csv"))
derivative_Flanker <- read_csv(paste0(interfileFolder, "/derivative_summary_Flanker.csv"))


# --- 1. Model Performance and Distribution Analysis ---
modFlanker.set1 <- modFlanker.set1.sum$performance.tb
modFlanker.set2 <- modFlanker.set2.sum$performance.tb
print(modFlanker.set1)
print(modFlanker.set2)

# Extract and print performance metrics
partial_Rsq_set1 <- modFlanker.set1$partialRsq
partial_Rsq_set2 <- modFlanker.set2$partialRsq
print(paste("R² for Set1: ", partial_Rsq_set1))
print(paste("R² for Set2: ", partial_Rsq_set2))

# Calculate skewness
skewness_value <- skewness(sumFlanker_deviation$nihtbx_flanker_uncorrected_deviationZ)
print(paste("Skewness: ", skewness_value))

# Calculate kurtosis
kurtosis_value <- kurtosis(sumFlanker_deviation$nihtbx_flanker_uncorrected_deviationZ)
print(paste("Kurtosis: ", kurtosis_value))
# Perform Shapiro-Wilk test on a random sample
set.seed(123)
sample_size <- 5000
sample_data <- sample(sumFlanker_deviation$nihtbx_flanker_uncorrected_deviationZ, size = sample_size, replace = FALSE)
shapiro_test_result <- shapiro.test(sample_data)
print(shapiro_test_result)


# --- 2. Plot Normative Trajectory and Derivative ---
# Calculate centiles from the full model
n_points <- 1000
centiles.all.tmp <- modFlanker.all.sum$centiles_strat  # full model centiles

Centiles <- (centiles.all.tmp[[1]] + centiles.all.tmp[[2]]) / 2

X  <- seq(min(sumFlanker_deviation$Age_year), max(sumFlanker_deviation$Age_year), length.out=1000)
sumFlanker_deviation$Sex <- factor(sumFlanker_deviation$Sex, levels = c("F", "M"))
# Plotting
p_main <- ggplot() +
  geom_hex(data = sumFlanker_deviation, aes(x = Age_year, y = nihtbx_flanker_uncorrected), bins = 40) +
  scale_fill_gradientn(colors = c("#E5E9F2", "#A4C5DF", "#4980B5")) + 
  geom_line(aes(x=X, y=Centiles[3,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[4,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[5,]), linetype="solid", color = "black", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[6,]), linetype="dashed", color = "gray10", size = 0.5) +
  geom_line(aes(x=X, y=Centiles[7,]), linetype="dashed", color = "gray10", size = 0.5) +
  scale_y_continuous(name="Score",limits = c(60, 120),breaks = seq(60, 120, by = 15))+
  scale_x_continuous(name="", limits = c(8.8, 16),breaks = seq(9, 16, by = 2))+
  labs(NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5, color = "black"),
    axis.text = element_text(colour = "black",size=9), 
    axis.title = element_text(size=8.5,face = "plain"),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.ticks = element_line(colour = "black", size = 0.25), 
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.y = element_text(colour = "black",size = 9),
    plot.background=element_rect(fill="white",color = NA),
    legend.position = "none" ,
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
x_length <- 7  
y_length <- 6
print(p_main)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_Flanker_blue.pdf"),
  plot = p_main,
  width = x_length, 
  height = y_length, 
  units = "cm"
)

# Bar plot for derivative significance
derivative_Flanker <- derivative_Flanker %>% drop_na(P50_lower, P50_upper)
derivative_Flanker <- derivative_Flanker %>%
  mutate(significance = ifelse(P50_lower > 0, "Increasing",
                               ifelse(P50_upper < 0, "Decreasing", "Non-significant")))
derivative_Flanker$color_group <- ifelse(derivative_Flanker$significance == "Non-significant", NA, derivative_Flanker$P50_mean)

p_bar <- ggplot(derivative_Flanker) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="#A4C5DF", low="white", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="#A4C5DF", low="white", midpoint=0, na.value="white", labels=NULL) +
  annotate("rect", xmin=min(derivative_Flanker$Age)-0.025, xmax=max(derivative_Flanker$Age)+0.025, ymin=-0.025, ymax=0.525, 
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
  filename = paste0(FigureFolder, "/NormativeDevCurve_Flanker_Addbar.pdf"), 
  plot = final_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)

# --- 3. Plot Variance Trajectory and Derivative ---
lower_CI <- derivative_Flanker$sigma_pred_lower
upper_CI <- derivative_Flanker$sigma_pred_upper

p_sigma <- ggplot(derivative_Flanker, aes(x = Age, y = sigma_pred_mean))+
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "gray80") +
  geom_line(size = 0.5, color = "black") +
  labs(NULL) +
  theme_minimal() +
  scale_y_continuous(name="Variance", limits = c(0.05, 0.08), breaks = seq(0.05, 0.08, by = 0.01)) +
  scale_x_continuous(name="", limits = c(8.8, 16), breaks = seq(9, 16, by = 2)) +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5),
    axis.text = element_text(colour = "black",size = 9),
    axis.title = element_text(size = 9, hjust = 0.5),
    axis.line = element_line(colour = "black", size = 0.25),
    axis.ticks = element_line(colour = "black", size = 0.25), 
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.y = element_text(colour = "black",size = 9),   
    legend.position = "none" ,
    legend.title=element_text(size=9),
    legend.text=element_text(size=9),
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
x_length <- 7
y_length <- 6
p_sigma
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_Flanker.pdf"),
  plot = p_sigma,
  width = x_length, 
  height = y_length,  
  units = "cm"
)
# Bar plot for variance derivative significance
derivative_Flanker <- derivative_Flanker %>% drop_na(sigma_deriv_lower, sigma_deriv_upper)
derivative_Flanker <- derivative_Flanker %>%
  mutate(significance = ifelse(sigma_deriv_lower > 0, "Increasing",
                               ifelse(sigma_deriv_upper < 0, "Decreasing", "Non-significant")))
derivative_Flanker$color_group <- ifelse(derivative_Flanker$significance == "Non-significant", NA, derivative_Flanker$sigma_deriv_mean)

p_sigmabar <- ggplot(derivative_Flanker) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="white", low="#A4C5DF", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="#A4C5DF", low="#A4C5DF", midpoint=0, na.value="white", labels=NULL) +
  annotate("rect", xmin=min(derivative_Flanker$Age)-0.025, xmax=max(derivative_Flanker$Age)+0.025, ymin=-0.025, ymax=0.525, 
           color="black", fill=NA, size=0.25) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab("Age (years)") +
  scale_x_continuous(breaks=NULL) +
  labs(fill="P50 mean") +
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
sigma_plot <- p_sigma / p_sigmabar + plot_layout(heights = c(y_length_main, y_length_bar))
print(sigma_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_Flanker_Addbar.pdf"), 
  plot = sigma_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)


