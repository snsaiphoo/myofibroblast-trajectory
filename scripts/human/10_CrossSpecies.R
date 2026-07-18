library(Hmisc)
library(dplyr)

comparison <- readRDS("../../data/CrossSpecies_PathwayMatrix.rds")

corr <- rcorr(
  as.matrix(comparison),
  type = "spearman"
)

mouse <- c(
  "Mouse_Repair_Fibroblasts",
  "Mouse_Repair_Activated_Stromal"
)

human <- c(
  "ADAM12hi fibroblasts",
  "NR4A1hi fibroblasts",
  "ABCA10hi fibroblasts",
  "FBLN1hi fibroblasts"
)

results <- expand.grid(
  Mouse = mouse,
  Human = human,
  stringsAsFactors = FALSE
)

results$Spearman_r <- mapply(function(x,y)
  corr$r[x,y],
  results$Mouse,
  results$Human)

results$P_value <- mapply(function(x,y)
  corr$P[x,y],
  results$Mouse,
  results$Human)

results <- results %>%
  mutate(
    Mouse = recode(Mouse,
                   Mouse_Repair_Fibroblasts = "Repair Fibroblasts",
                   Mouse_Repair_Activated_Stromal = "Activated Stromal"
    ),
    Human = recode(Human,
                   `ADAM12hi fibroblasts` = "ADAM12hi",
                   `NR4A1hi fibroblasts` = "NR4A1hi",
                   `ABCA10hi fibroblasts` = "ABCA10hi",
                   `FBLN1hi fibroblasts` = "FBLN1hi"
    ),
    Spearman_r = round(Spearman_r,2),
    P_value = signif(P_value,3)
  )

results

write.csv(
  results,
  "../../results/07_snRNA/CrossSpecies_SpearmanCorrelation.csv",
  row.names = FALSE
)

