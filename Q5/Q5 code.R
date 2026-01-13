
library(readxl)

setwd("~/Desktop")

geo <- read_excel(
  "GSE39645_GEO2R.xlsx",
  skip = 1
)
head(geo)

geo_clean <- geo[, c(
  "Gene.symbol",
  "logFC",
  "P.Value",
  "adj.P.Val",
  "t"
)]

colnames(geo_clean) <- c(
  "gene_symbol",
  "logFC",
  "p_value",
  "adj_p_value",
  "t_stat"
)

geo_clean$geo_id <- "GSE39645"

geo_clean$regulation <- ifelse(
  geo_clean$logFC > 0,
  "Up",
  "Down"
)

geo_sorted <- geo_clean[order(geo_clean$adj_p_value), ]
geo_top250 <- geo_sorted[1:250, ]

nrow(geo_top250)
head(geo_top250)
tail(geo_top250)



write.table(
  geo_top250,
  file = "Top250.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  na = ""
)

