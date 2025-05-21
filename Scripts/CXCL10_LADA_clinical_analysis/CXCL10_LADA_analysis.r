# 加载必要包
library(pROC)    # 用于AUC计算
library(caret)   # 用于数据分割和Bootstrap
library(tidyverse)
library(boot)
library(performanceEstimation)
library(precrec)

LADA_data <- read.csv('CXCL10_LADA_data.csv', stringsAsFactors = F, header = T)
### check the number of diabetes cases ###
table(LADA_data$group)

LADA_data$BMI <- as.numeric(LADA_data$BMI)
cor.test(LADA_data$CXCL10, LADA_data$BMI, method = 'spearman')

#### select healthy control and LADA as a cohort ####
#### the association between CXCL10 and LADA ###
data <- LADA_data[(LADA_data$group==0)|(LADA_data$group==3), ]
data$group <- ifelse(data$group==3, 1, 0)
data$group <- as.factor(data$group)
data$BMI <- as.numeric(data$BMI)
data <- data[!is.na(data$BMI), ]
# data$CXCL10 <- scale(data$CXCL10)
### CXCL10: 1.04 [1.03, 1.06]
model1 <- glm(group~Age+Sex+BMI+CXCL10, family = binomial(link = 'logit'), 
              data = data)
summary(model1)
### CXCL10: 1.05 [1.02, 1.10]
model2 <- glm(group~Age+Sex+BMI+FBG+HbA1c+CXCL10, family = binomial(link = 'logit'), 
              data = data)
summary(model2)

# 单变量逻辑回归模型（仅用CXCL10）
model <- glm(group ~ Age+Sex+BMI+CXCL10, data = data, family = binomial)

# 预测概率
pred_prob <- predict(model, type = "response")

# 计算AUC（需指定真实标签）
auc_result <- roc(data$group, pred_prob)
plot(auc_result, main = "ROC Curve for CXCL10")  # 可视化
auc_value <- auc_result$auc
cat("Baseline AUC:", round(auc_value, 3), "")

set.seed(123)
# 定义Bootstrap函数
bootstrap_auc <- function(data, indices) {
  # 有放回抽样（样本量与原数据相同）
  sampled_data <- data[indices, ]
  
  # 拟合模型（仅用CXCL10）
  model <- glm(group ~ Age+Sex+BMI+CXCL10, data = sampled_data, family = binomial)
  
  # 预测概率
  pred_prob <- predict(model, newdata = data, type = "response")
  
  # 计算AUC（使用原始数据作为测试集）
  roc_obj <- roc(data$group, pred_prob)
  return(roc_obj$auc)
}

# 执行Bootstrap（10000次重采样）
n_bootstraps <- 10000
boot_results <- boot(data = data, statistic = bootstrap_auc, R = n_bootstraps)

# 提取结果
boot_auc_mean <- mean(boot_results$t)
boot_ci <- boot.ci(boot_results, type = "bca")$bca  # BCa置信区间

cat("Bootstrap Results:
")
cat("Mean AUC:", round(boot_auc_mean, 3), "
")
cat("95% CI:", round(boot_ci[4], 3), "-", round(boot_ci[5], 3), "
")


#### 计算 PRROC

# 计算Precision-Recall AUC
pr_auc <- precrec::evalmod(labels = data$group, scores = pred_prob)
precrec::auc(pr_auc)$aucs[2]

# 定义PR曲线AUC的Bootstrap函数
bootstrap_pr_auc <- function(data, indices) {
  sampled_data <- data[indices, ]
  
  # 在新样本中拟合模型
  model <- glm(group ~ Age+Sex+BMI+CXCL10, data = sampled_data, family = binomial)
  
  # 使用原始数据进行预测（保持与ROC计算一致）
  pred_prob <- predict(model, newdata = data, type = "response")
  
  # 计算Precision-Recall AUC
  pr_auc <- precrec::evalmod(labels = data$group, scores = pred_prob)
  
  pr_auc <- precrec::auc(pr_auc)$aucs[2]
  
  return(pr_auc)
}

# 设置随机种子保证可重复性
set.seed(123)

# 执行Bootstrap重采样（10,000次）
n_bootstraps <- 10000
boot_results_pr <- boot(data = data, statistic = bootstrap_pr_auc, R = n_bootstraps)

# 提取结果
boot_pr_mean <- mean(boot_results_pr$t)
boot_ci_pr <- boot.ci(boot_results_pr, type = "bca")$bca  # BCa置信区间

# 输出结果
cat("Precision-Recall AUC Results:
")
cat("Mean AUC:", round(boot_pr_mean, 3), "
")
cat("95% CI:", round(boot_ci_pr[4], 3), "-", round(boot_ci_pr[5], 3), "
")



#### select healthy T2D and LADA as a cohort ####
#### the association between CXCL10 and T2D ###
data <- LADA_data[(LADA_data$group==2)|(LADA_data$group==3), ]
data$group <- ifelse(data$group==3, 1, 0)
data$BMI <- as.numeric(data$BMI)
data <- data[!is.na(data$BMI), ]
# data <- data[!is.na(data$FCP), ]

# 单变量逻辑回归模型（仅用CXCL10）
model <- glm(group ~ Age+Sex+BMI+CXCL10, data = data, family = binomial)

# 预测概率
pred_prob <- predict(model, type = "response")

# 计算AUC（需指定真实标签）
auc_result <- roc(data$group, pred_prob)
plot(auc_result, main = "ROC Curve for CXCL10")  # 可视化
auc_value <- auc_result$auc
cat("Baseline AUC:", round(auc_value, 3), "")

set.seed(123)
# 定义Bootstrap函数
bootstrap_auc <- function(data, indices) {
  # 有放回抽样（样本量与原数据相同）
  sampled_data <- data[indices, ]
  
  # 拟合模型（仅用CXCL10）
  model <- glm(group ~ Age+Sex+BMI+CXCL10, data = sampled_data, family = binomial)
  
  # 预测概率
  pred_prob <- predict(model, newdata = data, type = "response")
  
  # 计算AUC（使用原始数据作为测试集）
  roc_obj <- roc(data$group, pred_prob)
  return(roc_obj$auc)
}

# 执行Bootstrap（10000次重采样）
n_bootstraps <- 10000
boot_results <- boot(data = data, statistic = bootstrap_auc, R = n_bootstraps)

# 提取结果
boot_auc_mean <- mean(boot_results$t)
boot_ci <- boot.ci(boot_results, type = "bca")$bca  # BCa置信区间

cat("Bootstrap Results:
")
cat("Mean AUC:", round(boot_auc_mean, 3), "
")
cat("95% CI:", round(boot_ci[4], 3), "-", round(boot_ci[5], 3), "
")

#### 计算 PRROC

# 计算Precision-Recall AUC
pr_auc <- precrec::evalmod(labels = data$group, scores = pred_prob)
precrec::auc(pr_auc)$aucs[2]

# 定义PR曲线AUC的Bootstrap函数
bootstrap_pr_auc <- function(data, indices) {
  sampled_data <- data[indices, ]
  
  # 在新样本中拟合模型
  model <- glm(group ~ Age+Sex+BMI+CXCL10, data = sampled_data, family = binomial)
  
  # 使用原始数据进行预测（保持与ROC计算一致）
  pred_prob <- predict(model, newdata = data, type = "response")
  
  # 计算Precision-Recall AUC
  pr_auc <- precrec::evalmod(labels = data$group, scores = pred_prob)
  
  pr_auc <- precrec::auc(pr_auc)$aucs[2]
  
  return(pr_auc)
}

# 设置随机种子保证可重复性
set.seed(123)

# 执行Bootstrap重采样（10,000次）
n_bootstraps <- 10000
boot_results_pr <- boot(data = data, statistic = bootstrap_pr_auc, R = n_bootstraps)

# 提取结果
boot_pr_mean <- mean(boot_results_pr$t)
boot_ci_pr <- boot.ci(boot_results_pr, type = "bca")$bca  # BCa置信区间

# 输出结果
cat("Precision-Recall AUC Results:
")
cat("Mean AUC:", round(boot_pr_mean, 3), "
")
cat("95% CI:", round(boot_ci_pr[4], 3), "-", round(boot_ci_pr[5], 3), "
")


# cross-validation
# 检查模型在不同子样本中的表现（分层抽样）
data$group <- factor(data$group, levels=c(1, 0), labels = c('LADA', 'T2D'))
train_control <- trainControl(
  method = "boot",
  number = 1000,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  sampling = "smote"  # 处理不平衡数据（可选）
)

model_full <- train(
  group ~ Age + Sex + BMI + CXCL10,
  data = data,
  method = "glm",
  family = "binomial",
  trControl = train_control,
  metric = "ROC"
)

cat("Model 2 (Full) AUC:", model_full$results$ROC, "")
