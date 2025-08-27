# -------------------------------------------------------------------------
# 1. load library
# -------------------------------------------------------------------------
library(readxl)
library(dplyr)

# -------------------------------------------------------------------------
# 2. read daya
# -------------------------------------------------------------------------
GNGd_path <- "D:/datasets/yunfu/raw_data/cleaned_table/GNGd_table.xlsx"
oneback_path <- "D:/datasets/yunfu/raw_data/cleaned_table/oneback_table.xlsx"
twoback_path <- "D:/datasets/yunfu/raw_data/cleaned_table/twoback_table.xlsx"
demo_path <- "D:/datasets/yunfu/raw_data/demo/sum_desensitization_data.xlsx"

GNGd_table <- read_excel(GNGd_path)
oneback_table <- read_excel(oneback_path)
twoback_table <- read_excel(twoback_path)
demo_table <- read_excel(demo_path)

# -------------------------------------------------------------------------
# 3. data preprocess
# -------------------------------------------------------------------------
ethnicity_info <- demo_table %>%
  select(ON, Ethnicity = demo_ethnic_y)

GNGd_table <- inner_join(GNGd_table, ethnicity_info, by = "ON")
oneback_table <- inner_join(oneback_table, ethnicity_info, by = "ON")
twoback_table <- inner_join(twoback_table, ethnicity_info, by = "ON")


# -------------------------------------------------------------------------
# 4. function for data analysis
# -------------------------------------------------------------------------
calculate_stats <- function(df, ef_col_name) {
  
  # sum of samples
  N <- nrow(df)
  
  age_stat <- sprintf("%.2f (%.2f)", mean(df$Age_year, na.rm = TRUE), sd(df$Age_year, na.rm = TRUE))
  
  male_n <- sum(df$Sex == "M", na.rm = TRUE)
  male_pct <- (male_n / N) * 100
  sex_stat <- sprintf("%d (%.2f)", male_n, male_pct)
  
  right_handed_n <- sum(df$Hand == "右", na.rm = TRUE)
  right_handed_pct <- (right_handed_n / N) * 100
  right_handed_stat <- sprintf("%d (%.2f)", right_handed_n, right_handed_pct)
  
  left_handed_n <- sum(df$Hand == "左", na.rm = TRUE)
  left_handed_pct <- (left_handed_n / N) * 100
  left_handed_stat <- sprintf("%d (%.2f)", left_handed_n, left_handed_pct)
  
  han_n <- sum(df$Ethnicity == 1, na.rm = TRUE)
  han_pct <- (han_n / N) * 100
  han_stat <- sprintf("%d (%.2f)", han_n, han_pct)
  
  others_n <- N - han_n
  others_pct <- (others_n / N) * 100
  others_stat <- sprintf("%d (%.2f)", others_n, others_pct)
  
  ef_perf_stat <- sprintf("%.2f (%.2f)", mean(df[[ef_col_name]], na.rm = TRUE), sd(df[[ef_col_name]], na.rm = TRUE))
  
  emotional_problems_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_ES_sum, na.rm = TRUE), sd(df$SDQ_ES_sum, na.rm = TRUE))
  peer_problems_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_PP_sum, na.rm = TRUE), sd(df$SDQ_PP_sum, na.rm = TRUE))
  conduct_problems_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_CP_sum, na.rm = TRUE), sd(df$SDQ_CP_sum, na.rm = TRUE))
  hyperactivity_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_H_sum, na.rm = TRUE), sd(df$SDQ_H_sum, na.rm = TRUE))
  prosocial_behavior_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_PB_sum, na.rm = TRUE), sd(df$SDQ_PB_sum, na.rm = TRUE))
  
  return(c(
    N,
    age_stat,
    sex_stat,
    "", 
    right_handed_stat,
    left_handed_stat,
    "", 
    han_stat,
    others_stat,
    ef_perf_stat,
    "", 
    emotional_problems_stat,
    peer_problems_stat,
    conduct_problems_stat,
    hyperactivity_stat,
    prosocial_behavior_stat
  ))
}


# -------------------------------------------------------------------------
# 5. call function
# -------------------------------------------------------------------------
gng_stats <- calculate_stats(GNGd_table, "d_prime")
oneback_stats <- calculate_stats(oneback_table, "Oneback_acc")
twoback_stats <- calculate_stats(twoback_table, "Twoback_acc")

# -------------------------------------------------------------------------
# 6. merge table
# -------------------------------------------------------------------------
row_labels <- c(
  "N",
  "Age(years) (mean (SD))",
  "Sex = Male (%)",
  "Handedness (%)",
  "  Right-handed",
  "  Left-handed",
  "Ethnicity (%)",
  "  Han nationality",
  "  Others",
  "EF performance (mean (SD))",
  "Mental health (mean (SD))",
  "  Emotional problems",
  "  Peer problems",
  "  Conduct problems",
  "  Hyperactivity",
  "  Prosocial behavior"
)

summary_table <- data.frame(
  Characteristics = row_labels,
  `Participants with Go/No-Go task` = gng_stats,
  `Participants with 1-back task` = oneback_stats,
  `Participants with 2-back task` = twoback_stats,
  check.names = FALSE 
)

# -------------------------------------------------------------------------
# 7. print table
# -------------------------------------------------------------------------
print(summary_table)