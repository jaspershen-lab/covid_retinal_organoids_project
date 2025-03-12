library(r4projects)
library(tidymass)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

load("3_data_analysis/1-data_preparation/2-metabolome/metabolome_data.rda")

# 创建目录并设置工作路径
dir.create("3_data_analysis/2-data_cleaning/2-metabolome",
           recursive = TRUE,
           showWarnings = FALSE)
setwd("3_data_analysis/2-data_cleaning/2-metabolome/")

# 读取数据

variable_info <- extract_variable_info(metabolome_data)

# 定义代谢物的KEGG ID映射
kegg_mapping <- c(
  # NEG模式
  "guanosine_NEG" = "C00387",
  "Uracil_NEG" = "C00106",
  "D-glucose-6-phosphate_NEG" = "C00092",
  "mannose_NEG" = "C00159",
  "glucose_NEG" = "C00031",
  "D-Fructose-6-phosphate_NEG" = "C05345",
  "Uridine_NEG" = "C00299",
  "2-deoxyguanosine_NEG" = "C00330",
  "Thymidine_NEG" = "C00214",
  "Guanine_NEG" = "C00242",
  "Aspartic acid_NEG" = "C00049",
  "Cystine_NEG" = "C00491",
  "O-phospho-L-threonine_NEG" = "C12147",
  "2-deoxycytidine_NEG" = "C00881",
  "Cholesterol sulphate_NEG" = "C18043",
  "Phospho-choline_NEG" = "C00588",
  
  # POS模式
  "Glutamine area_POS" = "C00064",      
  "Alanine area_POS" = "C00041",        
  "Glutamic acid area_POS" = "C00025",  
  "Glycine area_POS" = "C00037",        
  "Lysine area_POS" = "C00047",         
  "Methionine area_POS" = "C00073",     
  "Isoleucine area_POS" = "C00407",     
  "leucine area_POS" = "C00123",        
  "Phenylalanine area_POS" = "C00079",  
  "Proline area_POS" = "C00148",        
  "Valine area_POS" = "C00183",         
  "urea area_POS" = "C00086",           
  "L-carnitine area_POS" = "C00487",    
  "2-deoxyadenosine area_POS" = "C00559", 
  "Cytosine area_POS" = "C00380",       
  "Tryptophan area_POS" = "C00078"      
)

# 重复ID以匹配数据中的重复项
KEGG_ID <- rep(kegg_mapping, each = 2)


# 首先只保留基础列
variable_info <- variable_info %>%
  dplyr::select(variable_id, Compound.name) %>%
  dplyr::mutate(KEGG_ID = KEGG_ID[match(Compound.name, names(KEGG_ID))])

# 转换KEGG到HMDB ID
HMDB_ID <- lapply(variable_info$KEGG_ID, function(x) {
  if(!is.na(x)) {
    trans_ID(x, "KEGG", "Human Metabolome Database", server = "cts.fiehnlab")
  } else {
    data.frame(KEGG = NA, `Human Metabolome Database` = NA)
  }
}) %>%
  do.call(rbind, .) %>%
  as.data.frame()

# 添加HMDB ID到variable_info
variable_info <- variable_info %>%
  dplyr::left_join(HMDB_ID, by = c("KEGG_ID" = "KEGG")) %>%
  dplyr::rename("HMDB_ID" = "Human Metabolome Database")

# 更新metabolome_data并进行标准化
metabolome_data@variable_info <- variable_info
metabolome_data <- metabolome_data %>% 
  normalize_data(method = "median")

# closeAllConnections()

# 保存数据
save(metabolome_data, file = "metabolome_data.rda")

