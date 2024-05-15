# Reproducibility for Lin, Sung and Chen (2024) Category Tree Gaussian Process for Computer Experiments with Many-Category Qualitative Factors and Application to Cooling System Design
This folder consists of the R code for the manuscript Category Tree Gaussian Process for Computer Experiments with Many-Category Qualitative Factors and Application to Cooling System Design by Lin, Sung and Chen (2024).

* The R script files reproduce the results in the Sections 5 and 6 of the manuscript. 
  * The required R packages include `tgp`, `DiceKriging`, `SLHD`, `LVGP`, `lhs`, `AdequacyModel`, `globpso`, `RcppArmadillo`, `plgp`, `foreach`, `doParallel`, `mvtnorm`, `parallel`, `dplyr`, `ggplot2`, `caret`, `GPareto` and `MASS`.
  * The `globpso` package can be installed via https://github.com/PingYangChen/globpso.
  * `TBGP.R` is the master file implementing the model introduced in sections 3 and 4.
  * `generatedSim.R` reproduces the results in section 5.1.
  * `5L.R` reproduces the results in section 5.2.
  * `corrExampleToy.R` reproduces the results in section 5.3.
  * `corrExampleBorehole3q.R` and `corrExampleBorehole3qUnion.R` reproduces the results in section 5.4.
  * While `realExample.R` reproduces the results in section 6, the real data was collected through a commercial engineering software which was under license of ANSYS Icepak, and therefore the data will not be shared publicly.
