# GOLDBioAge 中文使用说明

[English README](README.md) · [模拟 Gompertz 完整教程](docs/simulated-gompertz.md) · [NHANES 完整教程](docs/nhanes.md) · [Cox 专题教程](docs/cox-bioage.md) · [已发表系数](Coefficients.md)

GOLDBioAge 可以用生存随访数据训练 GOLD BioAge 模型，也可以将已经固定的系数用于新的生物标志物数据。

## 先确定要做哪一种计算

| 目的 | 所需数据 | 方法 |
|---|---|---|
| 在自己的队列重新训练模型 | 年龄、生物标志物、随访时间、事件状态 | `gold_bioage()` |
| 用训练好的模型计算新人群 | 年龄和入选生物标志物 | `predict_gold_bioage(fit, newdata)` |
| 直接使用已发表的 NHANES/UKB 系数 | 对应系数向量、年龄和相同单位的指标 | `predict_gold_bioage(coefficients, newdata)` |
| 使用 Cox 回归拟合年龄分数 | 年龄、生物标志物、随访时间、事件状态 | `cox_bioage()` |
| 将固定 Cox 模型应用到新数据 | 年龄和拟合时的指标 | `predict_cox_bioage(fit, newdata)` |

调用 `gold_bioage()` 会重新估计模型，不会自动使用已发表系数。调用 `predict_gold_bioage()` 使用固定系数，不需要随访结局，也不会在新样本中重新校准截距。

## 安装

```r
install.packages("remotes")
remotes::install_github("Jerryhaom/GOLDBioAge")
library(GOLDBioAge)

# 只有使用 LASSO / Elastic Net 筛选时才需要：
install.packages(c("glmnet", "ggplot2"))
```

## 案例一：模拟 Gompertz 生存数据

下面的年龄是入组时年龄，`time` 是从入组开始经过的随访时间，两者都以年为单位。Gompertz 随访风险为 `h(t) = rate × exp(shape × t)`，这里 `shape = 0.08`，所以风险随随访时间增加。

```r
library(GOLDBioAge)
set.seed(20261106)
n <- 2000L
age <- runif(n, 40, 80)
marker_risk <- rnorm(n)
marker_protective <- rnorm(n)
individual_rate <- exp(-9 + 0.08 * age +
                       0.4 * marker_risk - 0.3 * marker_protective)
event_time <- flexsurv::rgompertz(n, shape = 0.08, rate = individual_rate)
censor_time <- runif(n, 5, 20)
simulated <- data.frame(
  time = pmin(event_time, censor_time),
  status = as.integer(event_time <= censor_time),
  age, marker_risk, marker_protective
)

# 先按参与者划分训练集和测试集。
train_id <- sample.int(n, floor(0.8 * n))
variables <- c("age", "marker_risk", "marker_protective")
fit_sim <- gold_bioage(simulated[train_id, ], variables)
fit_sim$coefficients

# 预测只提供生物标志物，保留训练集的系数和截距。
test <- simulated[-train_id, ]
test$GOLDbioage <- predict_gold_bioage(fit_sim, test[, variables])
test$age_difference <- test$GOLDbioage - test$age
summary(test$age_difference)
```

这里的指标是人为构造的标准化数值。模拟案例用于展示计算流程，不能说明模型已经具备临床有效性。完整教程和 [可运行脚本](examples/01_simulated_gompertz.R) 还包含无效噪声指标、LASSO 和 Elastic Net。

## 案例二：使用包内 NHANES 数据重新拟合

```r
library(GOLDBioAge)
data("NHANES4", package = "GOLDBioAge")
variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")

# 检查样本量、事件编码和缺失值。
dim(NHANES4)
table(NHANES4$status)
colSums(is.na(NHANES4[, c("time", "status", variables)]))

fit_nhanes <- gold_bioage(NHANES4, variables)
round(fit_nhanes$coefficients, 8)
head(fit_nhanes$GOLDBioage[, c("age", "GOLDbioage")])

# 保持原始行顺序；缺失预测指标的行返回 NA。
NHANES4$GOLDbioage_refit <- predict_gold_bioage(fit_nhanes, NHANES4)
length(fit_nhanes$rows_used)
length(fit_nhanes$rows_omitted)
```

包内数据有 **54,945 行、30 列**。按上述指标进行完整病例拟合，使用 **39,348 人、4,716 例死亡**，因模型变量缺失排除 15,597 行。拟合不会自动插补，也没有加入 NHANES 调查权重、分层或聚类设计。这是分发文件的分析示例，不是从 NHANES 原始下载开始的完整复现流程。

## 案例三：使用 Cox 回归拟合

```r
library(GOLDBioAge)
data("NHANES4", package = "GOLDBioAge")
cox_variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
                   "lymph", "mcv", "rdw", "wbc", "ggt")
fit_cox <- cox_bioage(NHANES4, cox_variables)
fit_cox$coefficients       # 保留原有系数矩阵格式，单位为年。
fit_cox$intercept          # 训练集截距，预测时不能省略。
fit_cox$score_coefficients # 包含 beta0 的完整命名系数向量。

NHANES4$Coxbioage <- predict_cox_bioage(fit_cox, NHANES4)
head(fit_cox$CoxBioage[, c("age", "Coxbioage")])
```

Cox 方法先用仅含年龄的模型估计年龄效应，再将这个效应作为 offset 固定到年龄加生物标志物的模型中。将指标系数除以参考年龄系数，即得到以年为单位的权重。Cox 的基线风险函数不指定为 Gompertz 形式。

默认截距在训练集上计算，使训练集平均 Cox BioAge 等于平均实际年龄；预测新数据时保留训练截距。**0.2.1 修正了旧版 Cox 预测值中心化与年龄 offset 不一致的问题，因此默认分数水平会改变，拟合权重不变。** 如需复现旧版结果，使用 `cox_bioage(..., calibration = "legacy")`。已经保存的旧版拟合对象可直接传给 `predict_cox_bioage()`，仍使用其原截距。

Cox 分数不是绝对死亡概率，训练集内与 Gompertz 分数相关也不能证明两种方法的外部预测性能等价。

两套完整案例脚本均已加入相同训练样本上的 Cox 拟合、固定系数评分和 Gompertz/Cox 比较。详见 [Cox 专题教程](docs/cox-bioage.md) 和[在线公式推导](https://jerryhaom.github.io/GOLDBioAge/gompertz_cox_formulas.html)。最终 Cox 拟合与 GOLD 流程中的 Cox-LASSO 筛选是两个步骤，不能混为一谈。

## 案例四：直接使用 NHANES 已发表系数

如果只想算 GOLD BioAge，无需在自己的样本上重新拟合生存模型。先核对列名和单位，然后直接应用系数：

```r
published_nhanes <- c(
  age = 1, albumin = -4.73575224, alp = 0.02604091, creat = 5.26909180,
  glucose_mmol = 0.57969768, lymph = -0.20318306, mcv = 0.33888497,
  rdw = 2.64452031, wbc = 0.44588830, ggt = 0.06075152, beta0 = -53.62873
)
person <- data.frame(age = 50, albumin = 4, alp = 80, creat = 1,
                     glucose_mmol = 5, lymph = 30, mcv = 90, rdw = 13, wbc = 7, ggt = 25)
predict_gold_bioage(published_nhanes, person)
# 51.10203967 年；网页保留一位小数时显示 51.1 年。
```

| 指标 | 单位 |
|---|---|
| age | 年 |
| albumin | g/dL |
| creat | mg/dL |
| glucose_mmol | mmol/L |
| alp、ggt | U/L |
| mcv | fL |
| lymph、rdw | 百分数，30% 写成 30 |
| wbc | 1000 cells/uL |

NHANES GOLD 中，**LYM% 的系数为负，GGT 的系数为正**。程序不会自动转换单位或进行对数转换。已发表截距这里采用 `-53.62873`；重新拟合时可得到更多小数位，因此可能有极小精度差异。

## 可选：变量筛选

`var` 是指标列名向量：

```r
sim_variables <- c("age", "marker_risk", "marker_protective")
fit_lasso <- gold_bioage(
  simulated[train_id, ], var = sim_variables,
  feature_selection = TRUE, selection_method = "lasso", nfolds = 5, seed = 123
)
fit_lasso$selected_vars
print(fit_lasso$cv_plot)
```

筛选阶段使用仅包含候选生物标志物的 Cox 模型，在 `lambda.min` 取非零系数，最终 Gompertz 模型保留年龄。Elastic Net 可设置 `selection_method = "elasticnet", alpha = 0.5`。筛选需要至少两个候选生物标志物、严格大于 0 的随访时间和不少于 `nfolds` 的事件。变量筛选必须在训练集内完成。

## 理解输出与结果

- `fit$coefficients`：命名系数向量，包含截距 `beta0`。
- `fit$GOLDBioage$GOLDbioage`：拟合完整病例的分数，注意字段大小写。
- `fit$rows_used` / `fit$rows_omitted`：用于拟合或被排除的原始行位置。
- `predict_gold_bioage()`：按新数据原始行顺序返回数值向量，缺失指标返回 `NA`。

训练集的平均 GOLD BioAge 与平均实际年龄相等，是截距校准的结果。将模型转移到新数据时保留这一截距，不要求新队列的平均差值也为零。`GOLDbioage - age` 是年龄差值，不等同于对年龄回归后的残差，也不能单独证明模型的临床效用。对开发数据重新拟合或回代评分不属于独立验证。

## 运行完整案例

在克隆仓库的目录下，安装包后运行：

```r
source("examples/01_simulated_gompertz.R")
source("examples/02_nhanes.R")
```

脚本会保留拟合对象和评分结果；保存 CSV 或 RDS 的代码默认注释，按需要启用。详细输入与返回值可查看 `?gold_bioage`、`?predict_gold_bioage`、`?cox_bioage`、`?predict_cox_bioage`。版本变化见 [NEWS.md](NEWS.md)。

论文：[Hao et al., Advanced Science, 2025](https://doi.org/10.1002/advs.202501765)。
