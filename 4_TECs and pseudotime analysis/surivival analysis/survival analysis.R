library(tidyverse)
library(dplyr)
library(survival)
library(survminer)

# 
load("tcga_data.RData")
expr = log2(expr1+1)
tumor = gsub("-",".",samplesNT)
normal = gsub("-",".",samplesNT)

# 
sur = read.delim("TCGA-LIHC.survival.tsv") %>% 
  select(c(sample,OS,OS.time))
head(sur)
names(sur) = c("sample","status","time")
sur$sample = gsub("-",".",sur$sample)
head(sur)
head(expr)

# 
genes = c("PDGFRB","PGF","JUN","NR4A1")
expr$gene = rownames(expr)
select_expr = subset(expr, gene %in% genes) %>% 
  t() %>% 
  as.data.frame()
select_expr$sample = rownames(select_expr)
head(select_expr)

# 
data = merge(select_expr, sur, by = "sample")
head(data)
dim(data)


# 
dim(data) # 415   7
for (i in colnames(data)[2:5]){
  name = paste0(i,"_","group")
  # value = data[,i]
  data[, name] = ifelse(data[,i] > median(data[,i]), "high", "low")
}
head(data)

table(data$JUN_group)

# 
data_value = data.frame()
for (i in colnames(data)[2:5]){
  name = paste0(i,"_","group")
  formula_str <- paste("Surv(time, status) ~", name)
  cox_model1 <- coxph(as.formula(formula_str), data = data)
  data_value[i,"pvalue"] = summary(cox_model1)$coefficients[,5]
}
head(data_value)
dim(data_value) # 
data_value_diff = filter(data_value, pvalue < 0.05)
dim(data_value_diff) # 
head(data_value_diff)


for (i in colnames(data)[2:5]){
  name = paste0(i,"_","group")
  formula_str <- paste("Surv(time, status) ~", name)
  cox_model2 <- survfit(as.formula(formula_str), data = data)
  title = paste0(i," Survival Curve")
  p = ggsurvplot(cox_model2, # 
                 data = data,  # 
                 conf.int = TRUE, # 
                 pval = TRUE, # 
                 risk.table = TRUE, # 
                 surv.median.line = "hv", # 
                 # add.all = TRUE, # 
                 palette = "hue") +  # 
    ggtitle(title)
  print(p)
}
