source("scripts/00_setup.R")

plan(multicore)

dom <- "EUR"

dir_raw_data <- str_glue("{dir_pers_disk}/raw_data")


dir_raw_data %>% 
  list.files() %>% 
  str_subset("RegCM") %>% 
  str_subset("MPI") %>% #.[1] %>% {str_glue("{dir_raw_data}/{.}")} %>% read_ncdf(proxy = F) -> s
  .[1:3] %>% 
  {str_glue("{dir_raw_data}/{.}")} %>% 
  str_flatten(" ") -> ff


str_glue("cdo cat {ff} test.nc") %>% system()

"test.nc" %>% read_ncdf(ncsub = cbind(start = c(1,1,1), count = c(2,2,NA))) -> s
s %>% st_get_dimension_values("time") -> time
time %>% as_date() -> time

time[1]
time[length(time)]

#unlink("test.nc")



set_units(20, degC) %>% 
  set_units(K) %>% 
  drop_units() -> lim_k

str_glue("cdo gec,{lim_k} test.nc test_gec.nc") %>% system()

"test_gec.nc" %>% read_ncdf(ncsub = cbind(start = c(100,100,1), count = c(150,150,NA))) -> s





str_glue("cdo yearsum test_gec.nc test_ysum.nc") %>% system()

"test_ysum.nc" %>% read_ncdf(ncsub = cbind(start = c(100,100,1), count = c(150,150,NA))) -> s






str_glue("cdo -yearsum -gec,{lim_k} -cat {ff} test_chain.nc") %>% system()

"test_chain.nc" %>% read_ncdf() %>% slice(time, 2) %>% mapview::mapview()
"test_ysum.nc" %>% read_ncdf() %>% slice(time, 2) %>% mapview::mapview()
