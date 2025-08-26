library(readxl)
library(dplyr)
library(stringr)

# 设置路径
path <- "D:/datasets/yunfu/raw_data/twoback_data/raw_data"

# 获取所有xlsx文件
files <- list.files(path = path, pattern = "\\.xlsx$", full.names = TRUE)

# 存储符合条件的被试信息
result_list <- list()

# 遍历每个文件
for (file in files) {
  # 读取数据
  data <- readxl::read_excel(file, col_types = "text")
  
  # 按照 id 分组处理
  ids <- unique(data$id)
  
  for (id in ids) {
    subj_data <- data[data$id == id, ]
    
    # 筛选正式试次和测试试次
    formal_data <- subj_data[grepl("正式", subj_data$trial), ]
    test_data <- subj_data[grepl("测试", subj_data$trial), ]
    
    # 判断是否正式试次为空但测试试次不为空
    if (nrow(formal_data) == 0 && nrow(test_data) > 0) {
      result_list[[length(result_list) + 1]] <- data.frame(
        File = basename(file),
        ID = id,
        stringsAsFactors = FALSE
      )
    }
  }
}

# 转换为数据框并输出结果
if (length(result_list) > 0) {
  result_df <- do.call(rbind, result_list)
  print(result_df)
} else {
  cat("没有找到符合条件的被试。\n")
}