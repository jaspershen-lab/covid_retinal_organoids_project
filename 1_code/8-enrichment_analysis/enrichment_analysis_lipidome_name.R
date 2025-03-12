library(r4projects)
library(stringr)
library(dplyr)

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# 加载数据
load("3_data_analysis/2-data_cleaning/3-lipidome/lipidome_data.RData")
cluster_1 <- read.csv("3_data_analysis/7-heatmap_analysis/lipidome_heatmap_analysis/lipidome_data/clustering_results/cluster_1.csv", 
                 header = TRUE, 
                 stringsAsFactors = FALSE,
                 fileEncoding = "UTF-8")

# 创建输出目录
dir.create("3_data_analysis/8-enrichment_analysis/lipidome_data/lipid_minion_input", 
           recursive = TRUE, 
           showWarnings = FALSE)

# 定义路径常量
CLUSTER_PATH <- "3_data_analysis/7-heatmap_analysis/lipidome_heatmap_analysis/lipidome_data/clustering_results"
OUTPUT_PATH <- "3_data_analysis/8-enrichment_analysis/lipidome_data/lipid_minion_input"

# 改进的脂质名称标准化函数
standardize_lipid_name <- function(name) {
  # 初始检查
  if (is.null(name) || is.na(name) || name == "") {
    return(NA_character_)
  }
  
  # 去除前后空格
  name <- trimws(name)
  
  # 更新需要跳过的模式
  skip_patterns <- c(
    "^RIKEN|",
    "unknown|",
    "interference|",
    "\\[standard\\s*confirm\\]|",
    "\\[need\\s*to\\s*validate\\]|",
    "^SFTSV|",
    "^maybe|",
    "\\(interference\\)|",
    "^\\d{8}-\\w+-unknown-\\d+|",
    "\\[standardconfirm\\]|",
    "interference$|",
    "\\[needtovalidate\\]|",
    "\\+\\[M\\+Na\\]\\+|",  # 质谱信息
    "\\+\\[M\\+H\\]\\+|",   # 质谱信息
    "\\+\\[M\\+NH4\\]\\+|", # 氨加合物
    "^Lyso-PAF|",           # 非标准脂质命名
    "or|OR|",               # 含多种可能性
    "^Phosphocholine$|",    # 基团名称
    "and|",                 # 混合物
    "\\+Na|",              # 质谱加合物信息
    "eor"                  # 处理 eor 情况
  )
  
  if (str_detect(name, paste0(skip_patterns, collapse = ""))) {
    return(NA_character_)
  }
  
  # 处理特殊后缀
  suffix <- ""
  if (str_detect(name, "_(A|B|C)$")) {
    suffix <- str_extract(name, "_(A|B|C)$")
    name <- str_remove(name, "_(A|B|C)$")
  }
  
  # 修复数字和括号之间的问题
  name <- str_replace(name, "(\\d+)\\(:", "($1:")
  
  # 处理双键位置信息
  name <- str_remove(name, "\\([5-9]Z[^\\)]*\\)")
  name <- str_remove(name, "\\d+Z,\\d+Z[^\\)]*")
  
  # 标准化复杂脂质类名称的格式
  complex_lipids <- c("HBMP", "DGCC", "DGTS", "DGGA", "SQDG", "MGDG")
  for (lipid in complex_lipids) {
    if (str_detect(name, paste0("^", lipid))) {
      content <- str_remove(name, paste0("^", lipid, "\\s*"))
      content <- str_remove_all(content, "[\\(\\)]")
      content <- str_replace_all(content, "_", "/")
      name <- paste0(lipid, "(", content, ")")
      break
    }
  }
  
  # 处理e结尾到O-的转换
  if (str_detect(name, "e\\)$")) {
    type <- str_extract(name, "^[A-Z]+")
    content <- str_remove(name, paste0("^", type, "\\("))
    content <- str_remove(content, "e\\)$")
    name <- paste0(type, "(O-", content, ")")
  }
  
  # TAG统一为TG
  name <- str_replace(name, "^TAG", "TG")
  
  # 处理一般形式
  if (!str_detect(name, "\\(")) {
    type <- str_extract(name, "^[A-Z][A-Za-z0-9-]*")
    content <- str_remove(name, paste0("^", type, "\\s*"))
    content <- str_replace_all(content, "_", "/")
    name <- paste0(type, "(", content, ")")
  }
  
  # 处理不完整的脂质表示
  name <- str_replace(name, "/0\\)$", "/0:0)")
  name <- str_replace(name, "/0/0\\)$", "/0:0/0:0)")
  
  # 清理空格
  name <- str_replace_all(name, "\\s+", "")
  
  # 添加回特殊标记
  if (suffix != "") {
    name <- paste0(name, suffix)
  }
  
  return(name)
}

# 改进的文件写入函数
write_output_files <- function(lipids, base_path) {
  # 确保没有空值和重复
  lipids <- unique(lipids[!is.na(lipids)])
  
  # 构建文件路径
  txt_path <- paste0(base_path, ".txt")
  csv_path <- paste0(base_path, ".csv")
  
  # 写入TXT文件
  txt_content <- c('"Lipid"', sprintf('"%s"', lipids))
  writeLines(txt_content, txt_path)
  
  # 写入CSV文件
  csv_df <- data.frame(Lipid = lipids)
  write.csv(csv_df, csv_path, row.names = FALSE, quote = TRUE)
  
  # 返回文件路径
  return(list(txt = txt_path, csv = csv_path))
}

# 处理universe文件
process_universe <- function() {
  cat("Processing universe file...\n")
  
  # 从lipidome_data提取化合物名称
  compound_names <- lipidome_data@variable_info$Compound.name
  
  # 转换名称
  standardized_names <- sapply(compound_names, standardize_lipid_name)
  
  # 保存结果
  base_path <- file.path(OUTPUT_PATH, "lipids_universe")
  files <- write_output_files(standardized_names, base_path)
  
  # 打印信息
  final_count <- length(unique(standardized_names[!is.na(standardized_names)]))
  cat("Universe file created with", final_count, "lipids\n")
  cat("Files created:", files$txt, "and", files$csv, "\n")
}

# 处理cluster文件
process_cluster_file <- function(cluster_number) {
  # 构建文件路径
  input_file <- file.path(CLUSTER_PATH, paste0("cluster_", cluster_number, ".csv"))
  base_path <- file.path(OUTPUT_PATH, paste0("cluster", cluster_number, "_query"))
  
  if (!file.exists(input_file)) {
    warning(paste("File not found:", input_file))
    return(NULL)
  }
  
  # 读取数据
  data <- read.csv(input_file, stringsAsFactors = FALSE)
  
  # 标准化脂质名称
  standardized_names <- sapply(data$lipid_name, standardize_lipid_name)
  
  # 保存结果
  files <- write_output_files(standardized_names, base_path)
  
  # 打印处理信息
  final_count <- length(unique(standardized_names[!is.na(standardized_names)]))
  cat(sprintf("\nProcessed Cluster %d:\n", cluster_number))
  cat(sprintf("Original lipids: %d\n", nrow(data)))
  cat(sprintf("Converted unique lipids: %d\n", final_count))
  cat(sprintf("Files created: %s and %s\n", files$txt, files$csv))
  
  return(standardized_names[!is.na(standardized_names)])
}

# 批量处理所有cluster文件
process_all_clusters <- function(cluster_count) {
  for(i in 1:cluster_count) {
    result <- process_cluster_file(i)
  }
}

# 主程序
main <- function() {
  # 创建输出目录
  dir.create(OUTPUT_PATH, recursive = TRUE, showWarnings = FALSE)
  
  # 处理universe文件
  process_universe()
  
  # 处理cluster文件
  process_all_clusters(6)
}
# 运行主程序
main()