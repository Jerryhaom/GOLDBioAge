# GOLD BioAge <br>
This package measures biological aging based on Gompertz law. <br>
The example data was downloaded from the National Health and Nutrition Examination Survey (NHANES). 

# INSTALL <br>
install.packages("devtools") <br>
devtools::install_github("Jerryhaom/GOLDBioAge") <br>

# Examples
head(NHANES4) <br>
var <- c("age", "albumin", "alp", "creat","glucose_mmol","lymph","mcv", "rdw", "wbc", "ggt") <br>
bioage <- gold_bioage(NHANES4, var) <br>

# Citation <br>
Meng Hao et al. Gompertz law based biological age (GOLD BioAge): simple and practical biological aging clocks to capture morbidity and mortality risks. In preparation.

