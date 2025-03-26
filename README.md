# GOLD BioAge <br>
This package provides biological aging clocks measurements based on Gompertz or Cox regression. <br>
The example data was downloaded from the National Health and Nutrition Examination Survey (NHANES). 

# INSTALL <br>
install.packages("devtools") <br>
devtools::install_github("Jerryhaom/GOLDBioAge") <br>

# Examples
library(GOLDBioAge) <br>
head(NHANES4) <br>
var <- c("age", "albumin", "alp", "creat","glucose_mmol","lymph","mcv", "rdw", "wbc", "ggt") <br>
bioage <- gold_bioage(NHANES4, var) <br>

# Citation <br>
Meng Hao et al. Gompertz Law-Based Biological Age (GOLD BioAge): A Simple and Practical Measurement of Biological Aging to Capture Morbidity and Mortality Risks. doi: 10.1101/2024.11.14.24317305

