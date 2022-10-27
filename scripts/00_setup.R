
library(tidyverse)
library(lubridate)
library(stars)
library(furrr)
library(units)

options(future.fork.enable = T)


dir_bucket_cmip5 <- "/mnt/bucket_cmip5"
dir_bucket_mine <- "/mnt/bucket_mine"
dir_pers_disk <- "/mnt/pers_disk"
