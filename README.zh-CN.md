# GOLDBioAge：从公式到可运行的 R 代码

本包可以用 Gompertz 回归拟合 **GOLD BioAge**，用 Cox 回归拟合另一种年龄评分，也可以把已有系数应用到新参与者的数据。

[English](README.md) · [数学公式与代码](docs/mathematics-and-code.md) · [各数据库系数](Coefficients.md) · [更新与版本迁移](docs/package-updates.md) · [在线 calculator](https://jerryhaom.github.io/GOLDBioAgeWeb.io/)

## 先选择你要做的事

| 你的目标 | 使用方法 | 需要的数据 |
|---|---|---|
| 用原有 NHANES 或 UKB GOLD 模型计算年龄 | 固定系数 + `predict_gold_bioage()` | 实际年龄及 9 项指标 |
| 用原有 NHANES、UKB 或 CHARLS Light 模型 | 固定系数 + 正确的 CRP 转换 | 实际年龄、肌酐、血糖、CRP |
| 用自己的队列建立 GOLD 或 Cox 模型 | `gold_bioage()` 或 `cox_bioage()` | 年龄、指标、随访时间、结局 |
| 把已经拟合的模型用于新数据 | 对应的 `predict_*()` 函数 | 模型所需指标；沿用训练截距 |

所有这些年龄评分，最终都通过一个加权和计算：

$$\mathrm{BioAge}=A+\sum_jw_jX_j+b_0,\qquad
\mathrm{BioAgeDiff}=\mathrm{BioAge}-A.$$

$A$ 是实际年龄（岁），$X_j$ 是指定单位和转换下的指标，$w_j$ 是换算为年龄单位的权重，$b_0$ 是截距。年龄的系数为 1。BioAgeDiff 是年龄之差，不是回归残差，也不能直接解释为剩余寿命。

## 1. 安装

```r
install.packages("remotes")
remotes::install_github("Jerryhaom/GOLDBioAge")
library(GOLDBioAge)
packageVersion("GOLDBioAge")
```

只有需要 LASSO/弹性网筛选时，才需额外安装：

```r
install.packages(c("glmnet", "ggplot2"))
```

## 2. 先算一个人：NHANES GOLD

下面是完整代码，直接使用原有 NHANES 系数，不重新拟合：

```r
library(GOLDBioAge)
beta <- c(age = 1, creat = 5.26909180, glucose_mmol = 0.57969768,
          mcv = 0.33888497, rdw = 2.64452031, albumin = -4.73575224,
          alp = 0.02604091, lymph = -0.20318306, wbc = 0.44588830,
          ggt = 0.06075152, beta0 = -53.62873)
person <- data.frame(age = 50, creat = 1, glucose_mmol = 5,
                     mcv = 90, rdw = 13, albumin = 4,
                     alp = 80, lymph = 30, wbc = 7, ggt = 25)
score <- predict_gold_bioage(beta, person)
c(BioAge = score, BioAgeDiff = score - person$age)
# BioAge = 51.10203967 岁，BioAgeDiff = 1.10203967 岁。
```

| 指标 | R 列名 | 必须使用的单位 |
|---|---|---|
| 实际年龄 | `age` | 岁 |
| 肌酐 | `creat` | mg/dL |
| 血糖 | `glucose_mmol` | mmol/L |
| 平均红细胞体积 | `mcv` | fL |
| 红细胞分布宽度 | `rdw` | % |
| 白蛋白 | `albumin` | g/dL |
| 碱性磷酸酶 | `alp` | U/L |
| 淋巴细胞百分比 | `lymph` | %，30% 输入 30 |
| 白细胞计数 | `wbc` | 1000 cells/uL |
| γ-谷氨酰转移酶 | `ggt` | U/L |

NHANES GOLD 的 **LYM% 权重为负，GGT 权重为正**。截距只加一次。批量计算时，每行一个参与者，使用同样列名和单位即可。预测保留原行序；缺少所需指标的行返回 `NA`，不会自动转换单位或填补缺失。

## 3. 原来各数据库的系数在哪里

[Coefficients.md](Coefficients.md) 保留原有五套系数，并给出公式、单位、计算案例：NHANES/UKB 的 GOLD，以及 NHANES/UKB/CHARLS 的 Light。另列一张表展示包内 NHANES 的 Cox 重新拟合结果。

0.2.2 起，五套原有系数也随包安装，可直接读取。例如切换到 UKB GOLD：

```r
reference <- read.csv(system.file("extdata", "published_coefficients.csv",
                                  package = "GOLDBioAge"))
selected <- subset(reference, dataset == "ukb" & model == "gold")
beta_ukb <- setNames(selected$coefficient, selected$term)
predict_gold_bioage(beta_ukb, person)
# 50.7472 岁。
```

Light 的公式是：

$$\mathrm{Light}=A+w_C\mathrm{CREA}+w_G\mathrm{GLU}
++w_R\ln(1+\mathrm{CRP})+b_0.$$

先换算原始浓度，再取自然对数。NHANES 要求 CRP 为 **mg/dL**，UKB/CHARLS 要求 **mg/L**；`log1p(0.3)` 与 `log1p(3)` 不同。在相同浓度 CRP = 0.3 mg/dL = 3 mg/L 下：

```r
nhanes_light <- 50 + 8.3313 * 1 + 0.8270 * 5 + 5.7305 * log1p(0.3) - 13.5298
ukb_light <- 50 + 9.0873 * 1 + 1.7311 * 5 + 4.1689 * log1p(3) - 20.1395
c(NHANES = nhanes_light, UKB = ukb_light)
# 50.43997842 岁与 53.38262256 岁。
```

在线 calculator 对两个数据库均要求输入 mg/dL，并在 UKB Light 内部乘以 10。R 预测函数要求你事先完成转换。不同数据库模型有不同权重和截距，结果不同是正常的。

## 4. 重新拟合 NHANES：Gompertz 与 Cox

```r
data("NHANES4", package = "GOLDBioAge")
variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")
gold_fit <- gold_bioage(NHANES4, variables)
cox_fit <- cox_bioage(NHANES4, variables)

gold_fit$coefficients
cox_fit$score_coefficients
NHANES4$GOLD_refit <- predict_gold_bioage(gold_fit, NHANES4)
NHANES4$Cox_refit <- predict_cox_bioage(cox_fit, NHANES4)
```

包内数据为 54,945 行；上述指标组合的完整样本为 39,348 人、4,716 例死亡。示例进行不加调查权重的完整案例拟合。重新拟合产生的系数与保留有限小数的原有固定模型会有细微差别。

两种方法都先拟合仅含年龄的模型，得到年龄效应 $\beta_A$；再把完整模型的年龄效应固定为它，将指标系数 $\theta_j$ 换算为年，并在训练样本计算截距：

$$w_j=\theta_j/\beta_A,\qquad
b_0=\operatorname{mean}_{\mathrm{train}}\left(A-\left[A+\sum_jw_jX_j\right]\right).$$

所以训练样本平均 BioAge 与平均实际年龄相同。新数据必须沿用这个截距，不能为了得到零均差再重新居中。Gompertz 的随访形状参数与年龄效应是两个参数；完整模型固定的是年龄效应。GOLD 选择基线瞬时风险（$t=0$），此时两个模型的形状因子都等于 1，所以计算公式不含随访形状，也不需要固定共同形状。形状参数仍参与拟合，固定后重新回归可能改变其他系数。Cox 通过 offset 固定年龄效应，基线风险函数不指定。

自己的数据需要数值型指标、以岁表示的 `age`、单位明确的非负随访 `time`，以及 0 = 删失、1 = 事件的 `status`。`rows_used` 和 `rows_omitted` 可用于追踪拟合时保留和排除的行。训练、筛选和预处理仅在开发样本内进行。

## 5. 建议阅读顺序

1. [公式与代码对照](docs/mathematics-and-code.md)：每个符号对应什么、加权和怎样算、权重和截距怎样得到。
2. [模拟 Gompertz 教程](docs/simulated-gompertz.md)：生成事件和删失，分训练/测试，拟合 GOLD/Cox 后固定系数预测；[完整脚本](examples/01_simulated_gompertz.R)。
3. [NHANES 教程](docs/nhanes.md)：检查数据、复现系数、保留行序、固定模型计算；[完整脚本](examples/02_nhanes.R)。
4. [原有数据库系数与案例](Coefficients.md)：五套模型的单位、转换、截距；[完整脚本](examples/03_existing_coefficients.R)。
5. [Cox 专门教程](docs/cox-bioage.md)：推导、校准、保存模型、旧版本兼容；[在线公式页](https://jerryhaom.github.io/GOLDBioAge/gompertz_cox_formulas.html)。

## 6. 更新说明

**0.2.2** 增加随包安装的系数文件，重组中英文说明，并将公式与可运行代码放在一起；保留原有基线算法和拟合结果。**0.2.1** 增加 Cox 固定模型预测，修正默认截距校准；使用 `calibration = "legacy"` 可复现较早版本的分数。**0.2.0** 改进输入校验、变量顺序、GOLD 固定预测和可选筛选，修复 Cox 筛选遗漏首个指标的问题。

筛选结果和旧 Cox 分数的整体水平可能因更新而变化；无筛选的 GOLD 结果保持到数值精度。详见[更新与迁移](docs/package-updates.md)、[NEWS](NEWS.md) 和[一致性检查记录](docs/consistency-audit.md)。

引用：Hao M. et al. *Advanced Science* (2025), e01765，[doi:10.1002/advs.202501765](https://doi.org/10.1002/advs.202501765)。联系：haombio@gmail.com。
