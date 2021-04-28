# To transfer data/spdemand_3639id336tow.rda -> data/spdemand_3639id336tow.feather with rownames as ID_TOW
library(tidyverse)
library(dtplyr)
load("data/spdemand_3639id336tow.rda") # data.table
nid <- 3639
ntow <- 336
maxid <- sort(unique(spdemand[,id]))[nid]
spdemand <- 
  spdemand[id <= maxid & tow <= ntow,] %>%
  # spdemand %>% 
  lazy_dt() %>% 
  # filter(id == 1003) %>%
  mutate(index = paste0(id, "_", tow)) %>% 
  # as_tibble() %>% 
  column_to_rownames(var = "index") %>% 
  select(-id, -tow)
# save(spdemand, file = "data/spdemand_3639id336tow_rowname.rda")
feather::write_feather(spdemand, paste0("data/spdemand_", nid, "id", ntow, "tow.feather"))


## Below are the running time testing for feather
# # Install feather
# devtools::install_github("wesm/feather/R")

library(feather)
path <- "~/git/kderm/data/spdemand_3id10tow.feather"
write_feather(spdemand, path)
tictoc::tic()
# df <- read_feather(path) # 3.578s
df <- load("data/spdemand_3639id336tow.rda") # 7.098s
tictoc::toc()

# .feather is faster but the file size is larger than .rda
library(feather)

x <- runif(1e7)
x[sample(1e7, 1e6)] <- NA # 10% NAs
df <- as.data.frame(replicate(10, x))
write_feather(df, 'test.feather')

system.time(read_feather('test.feather'))
#>   user  system elapsed
#>  0.731   0.287   1.020

# $ pip install feather-format

# import feather
# path = 'your_file_path'
# datafile = feather.read_dataframe(path)

# import feather
# import pandas as pd
# import numpy as np
# arr = np.random.randn(10000000) # 10% nulls
# arr[::10] = np.nan
# df = pd.DataFrame({'column_{0}'.format(i): arr for i in range(10)})
# feather.write_dataframe(df, 'test.feather')
