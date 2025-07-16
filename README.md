# GOLD BioAge <br>
This package provides measurements of biological aging clocks based on Gompertz or Cox regression. <br>
The example data was downloaded from the National Health and Nutrition Examination Survey (NHANES). 

# INSTALL <br>
install.packages("devtools") <br>
devtools::install_github("Jerryhaom/GOLDBioAge") <br>

# Examples
library(GOLDBioAge) <br>
head(NHANES4) <br>
var <- c("age", "albumin", "alp", "creat","glucose_mmol","lymph","mcv", "rdw", "wbc", "ggt") <br>
#calculate bioage based on Gompertz models <br>
bioage <- gold_bioage(NHANES4, var) <br>
#calculate bioage based on cox regressions <br>
bioage <- cox_bioage(NHANES4, var) <br>

# More information
https://github.com/Jerryhaom/GOLDBioAgeWeb.io <br>

Contact Information: haombio@gmail.com <br>
Citation: M. Hao, H. Zhang, J. Wu, Y. Huang, X. Li, M. Wang, S. Wang, J. Wang, J. Chen, Z. jun Bao, L. Jin, X. Wang, Z. Hu, S. Jiang, Y. Li, Gompertz Law-Based Biological Age (GOLD BioAge): A Simple and Practical Measurement of Biological Ageing to Capture Morbidity and Mortality Risks. Adv. Sci. 2025, e01765. https://doi.org/10.1002/advs.202501765 <br>
