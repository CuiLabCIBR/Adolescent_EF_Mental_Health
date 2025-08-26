# -------------------------------------------------------------------------
# 1. 加载所需的库
# -------------------------------------------------------------------------
# 如果没有安装，请先运行: install.packages(c("readxl", "dplyr"))
library(readxl)
library(dplyr)

# -------------------------------------------------------------------------
# 2. 读取所有数据文件
# -------------------------------------------------------------------------
GNGd_path <- "D:/datasets/yunfu/raw_data/cleaned_table/GNGd_table.xlsx"
oneback_path <- "D:/datasets/yunfu/raw_data/cleaned_table/oneback_table.xlsx"
twoback_path <- "D:/datasets/yunfu/raw_data/cleaned_table/twoback_table.xlsx"
demo_path <- "D:/datasets/yunfu/raw_data/demo/执行功能结果变量计分_总_脱敏.xlsx"

GNGd_table <- read_excel(GNGd_path)
oneback_table <- read_excel(oneback_path)
twoback_table <- read_excel(twoback_path)
demo_table <- read_excel(demo_path)

# -------------------------------------------------------------------------
# 3. 数据预处理：合并民族信息
# -------------------------------------------------------------------------
# 从 GNG 人口学数据中选择 ON 和 民族 列，并重命名
ethnicity_info <- demo_table %>%
  select(ON, Ethnicity = demo_ethnic_y)

# 使用内连接 (inner_join) 将民族信息合并到主数据表中
# 这会保留两个表中都存在的 ON 值的行
GNGd_table <- inner_join(GNGd_table, ethnicity_info, by = "ON")
oneback_table <- inner_join(oneback_table, ethnicity_info, by = "ON")
twoback_table <- inner_join(twoback_table, ethnicity_info, by = "ON")


# -------------------------------------------------------------------------
# 4. 定义一个函数来计算每个数据集的统计数据
# -------------------------------------------------------------------------
calculate_stats <- function(df, ef_col_name) {
  
  # 总样本量
  N <- nrow(df)
  
  # 年龄 (均值, 标准差)
  age_stat <- sprintf("%.2f (%.2f)", mean(df$Age_year, na.rm = TRUE), sd(df$Age_year, na.rm = TRUE))
  
  # 性别 = 男性 (%)
  # 假设: 性别列为字符型, 男性记为 "M"
  male_n <- sum(df$Sex == "M", na.rm = TRUE)
  male_pct <- (male_n / N) * 100
  sex_stat <- sprintf("%d (%.2f)", male_n, male_pct)
  
  # 利手 (%)
  # 假设: 利手列为字符型, 右利手记为 "右", 左利手记为 "左"
  right_handed_n <- sum(df$Hand == "右", na.rm = TRUE)
  right_handed_pct <- (right_handed_n / N) * 100
  right_handed_stat <- sprintf("%d (%.2f)", right_handed_n, right_handed_pct)
  
  left_handed_n <- sum(df$Hand == "左", na.rm = TRUE)
  left_handed_pct <- (left_handed_n / N) * 100
  left_handed_stat <- sprintf("%d (%.2f)", left_handed_n, left_handed_pct)
  
  # 民族 (%)
  # 使用合并后的 'Ethnicity' 列, 汉族编码为 1
  han_n <- sum(df$Ethnicity == 1, na.rm = TRUE)
  han_pct <- (han_n / N) * 100
  han_stat <- sprintf("%d (%.2f)", han_n, han_pct)
  
  others_n <- N - han_n
  others_pct <- (others_n / N) * 100
  others_stat <- sprintf("%d (%.2f)", others_n, others_pct)
  
  # EF 表现 (均值, 标准差)
  ef_perf_stat <- sprintf("%.2f (%.2f)", mean(df[[ef_col_name]], na.rm = TRUE), sd(df[[ef_col_name]], na.rm = TRUE))
  
  # 心理健康 (均值, 标准差)
  emotional_problems_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_ES_sum, na.rm = TRUE), sd(df$SDQ_ES_sum, na.rm = TRUE))
  peer_problems_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_PP_sum, na.rm = TRUE), sd(df$SDQ_PP_sum, na.rm = TRUE))
  conduct_problems_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_CP_sum, na.rm = TRUE), sd(df$SDQ_CP_sum, na.rm = TRUE))
  hyperactivity_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_H_sum, na.rm = TRUE), sd(df$SDQ_H_sum, na.rm = TRUE))
  prosocial_behavior_stat <- sprintf("%.2f (%.2f)", mean(df$SDQ_PB_sum, na.rm = TRUE), sd(df$SDQ_PB_sum, na.rm = TRUE))
  
  # 将所有结果作为一个向量返回
  return(c(
    N,
    age_stat,
    sex_stat,
    "", # 为 Handedness (%) 留空
    right_handed_stat,
    left_handed_stat,
    "", # 为 Ethnicity (%) 留空
    han_stat,
    others_stat,
    ef_perf_stat,
    "", # 为 Mental health (mean(SD)) 留空
    emotional_problems_stat,
    peer_problems_stat,
    conduct_problems_stat,
    hyperactivity_stat,
    prosocial_behavior_stat
  ))
}


# -------------------------------------------------------------------------
# 5. 为每个数据集计算统计数据
# -------------------------------------------------------------------------
gng_stats <- calculate_stats(GNGd_table, "d_prime")
oneback_stats <- calculate_stats(oneback_table, "Oneback_acc")
twoback_stats <- calculate_stats(twoback_table, "Twoback_acc")

# -------------------------------------------------------------------------
# 6. 创建最终的汇总表格
# -------------------------------------------------------------------------
# 定义行名
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

# 合并成数据框
summary_table <- data.frame(
  Characteristics = row_labels,
  `Participants with Go/No-Go task` = gng_stats,
  `Participants with 1-back task` = oneback_stats,
  `Participants with 2-back task` = twoback_stats,
  check.names = FALSE # 防止 R 自动修改列名
)

# -------------------------------------------------------------------------
# 7. 打印结果
# -------------------------------------------------------------------------
print(summary_table)