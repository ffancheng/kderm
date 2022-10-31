rm(list = ls())
library(tidyverse)
library(dimRed)
library(dtplyr)
library(data.table)
library(intRinsic)
Jmisc::sourceAll("R/sources")
scen <- as.numeric(commandArgs()[[6]])
method <- c("mle", "linfit", "bayes")[scen]

load("data/half_count_ratio_3639id336tow.rda")
nid <- 3639
ntow <- 336

if(nid == 1) {
  train <- spdemand %>%
    lazy_dt() %>%
    filter(tow <= ntow,
           # id <= sort(unique(spdemand[,id]))[nid]
           id == 1003
    ) %>%
    dplyr::select(-id, -tow) %>%
    as.data.table()
} else if(nid == 3639) {
  train <- spdemand %>%
    lazy_dt() %>%
    # mutate(id_tow = paste0(id, "_", tow)) %>%
    as.data.table() %>%
    column_to_rownames("id")
}
dim(train)
# # saveRDS(train, file = paste0("data/spdemand_", nid, "id", ntow, "tow_train.rds"))

N <- nrow(train)
d <- intRinsic::twonn(train, method = method) # 6 for 1ID, 7 for all spdemand with 3639*336 rows, 91 for 3639 rows
print(d)
