
source("scripts/00_setup.R")
source("scripts/functions.R")

plan(multicore)


dom <- "EUR"


# DOWNLOAD RAW DATA -------------------------------------------------------------------------------

dir_raw_data <- str_glue("{dir_pers_disk}/raw_data")

{
  # dir.create(dir_raw_data)
  # 
  # # Download REMO data
  # 
  # dir_tmp <- str_glue("{dir_bucket_cmip5}/RCM_regridded_data/REMO2015/{dom}/daily/maximum_temperature")
  # 
  # dir_tmp %>% 
  #   list.files() %>% 
  #   .[str_length(.) > 80] %>% 
  #   
  #   future_walk(function(f){
  #     
  #     
  #     str_glue("{dir_tmp}/{f}") %>% 
  #       {system(str_glue("gsutil cp {.} {dir_raw_data}"), ignore.stdout = T, ignore.stderr = T)}
  #     
  #     
  #   })
  # 
  # 
  # # Download RegCM4 data
  # 
  # dir_tmp <- str_glue("{dir_bucket_cmip5}/RCM_regridded_data/CORDEX_22/{dom}/daily/maximum_temperature")
  # 
  # dir_tmp %>% 
  #   list.files() %>% 
  #   .[str_length(.) > 80] %>% 
  #   str_subset("EC-EARTH", negate = T) %>% 
  #   str_subset("CERFACS", negate = T) %>% 
  #   
  #   future_walk(function(f){
  #     
  #     
  #     str_glue("{dir_tmp}/{f}") %>% 
  #       {system(str_glue("gsutil cp {.} {dir_raw_data}"), ignore.stdout = T, ignore.stderr = T)}
  #     
  #     
  #   })
}





# TILE DATA ---------------------------------------------------------------------------------------

f <-
  str_glue("{dir_pers_disk}/raw_data") %>%
  list.files(full.names = T) %>%
  str_subset("REMO") %>% 
  .[1]

source("scripts/tiling.R")




# THRESHOLD TABLE ---------------------------------------------------------------------------------

str_glue("{dir_bucket_mine}/misc_data/CMIP5_model_temp_thresholds.csv") %>% 
  read_delim() %>% #filter(str_detect(Model, "MPI"))
  select(1:6) %>% 
  pivot_longer(-Model, names_to = "wl") %>% 
  mutate(wl = str_sub(wl, 3)) -> thresholds

thresholds %>% 
  mutate(wl = ifelse(str_length(wl) == 1, str_glue("{wl}.0"), wl)) -> thresholds

thresholds %>% 
  mutate(Model = case_when(str_detect(Model, "HadGEM") ~ str_glue("MOHC-{Model}"),
                           str_detect(Model, "MPI") ~ str_glue("MPI-M-{Model}"),
                           str_detect(Model, "NorESM") ~ str_glue("NCC-{Model}"),
                           TRUE ~ Model)) -> thresholds




# MODELS TABLE ------------------------------------------------------------------------------------

dir_raw_data %>% 
  list.files() %>% 
  str_split("_", simplify = T) %>% 
  .[,c(6,3)] %>% 
  as_tibble(.name_repair = ~c("rcm", "gcm")) %>%
  {unique(.[c("rcm", "gcm")])} -> tb_models


tb_models %>% 
  mutate(calendar = case_when(str_detect(gcm, "Had") ~ "360",
                              str_detect(rcm, "REMO") & str_detect(gcm, "Had", negate = T) ~ "gregorian",
                              str_detect(rcm, "RegCM") & str_detect(gcm, "MPI") ~ "gregorian",
                              str_detect(rcm, "RegCM") & str_detect(gcm, "Nor") ~ "noleap")) -> tb_models





# PRE-PROCESS FOR BOOTSTRAP -----------------------------------------------------------------------

set_units(32, degC) %>% 
  set_units(K) %>% 
  drop_units() -> lim_k


dir_tiles <- str_glue("{dir_pers_disk}/tiles")
dir.create(dir_tiles)


# Function to obtain statistics
fn_statistics <- function(ts, index, j = F){
  
  if(any(is.na(ts))){
    
    c(mean = NA,
      p05 = NA, 
      p50 = NA,
      p95 = NA)
    
  } else {
    
    sub_ts <- ts[index]
    
    if(isTRUE(j)){
      sub_ts <- jitter(sub_ts)
    }
    
    c(mean = mean(sub_ts),
      quantile(sub_ts, c(0.05, 0.5, 0.95)) %>% setNames(c("p05", "p50", "p95")))
    
  }
}



# LOOP THROUGH TILES

for(i in seq_len(nrow(chunks_ind))){
  
  print(str_glue(" "))
  print(str_glue("*******************"))
  print(str_glue("PROCESSING TILE {i}"))
  
  
  r <- chunks_ind$r[i]
  lon_ch <- chunks_ind$lon_ch[i]
  lat_ch <- chunks_ind$lat_ch[i]
  
  
  
  # READ DATA
  
  cbind(start = c(lon_chunks[[lon_ch]][1], lat_chunks[[lat_ch]][1], 1),
        count = c(lon_chunks[[lon_ch]][2] - lon_chunks[[lon_ch]][1]+1,
                  lat_chunks[[lat_ch]][2] - lat_chunks[[lat_ch]][1]+1,
                  NA)) -> ncs
  
  
  tb_models %>% 
    pmap(function(rcm, gcm, calendar){
      
      print(str_glue("Importing {rcm} - {gcm}"))
      
      dir_raw_data %>% 
        list.files() %>% 
        str_subset(rcm) %>% 
        str_subset(gcm) %>% 
        
        future_map(function(f){
          
          str_c(dir_raw_data, "/", f) %>% 
            read_ncdf(proxy = T) %>% 
            suppressMessages() %>% 
            st_get_dimension_values("time") %>% 
            as.character() -> time_dim
          
          fn_dates(time_dim, calendar) -> fixed_time
          
          str_c(dir_raw_data, "/", f) %>% 
            read_ncdf(ncsub = ncs) %>% 
            suppressMessages() %>% 
            st_set_dimensions("time", values = fixed_time)
          
        }) %>% 
        do.call(c, .)
      
    }) -> l_s
  
  
  
  # LOOP THROUGH WL 
  
  walk(c("1.0", "3.0"), function(wl_){
    
    print(str_glue(" "))
    print(str_glue("Processing WL {wl_}"))
    
    
    # Slice WL years and aggregate by years
    
    print(str_glue("Slicing + Aggregating"))
    
    future_map2(l_s, seq(1,6), function(s, i){
      
      rcm <- tb_models$rcm[i]
      gcm <- tb_models$gcm[i]
      
      thresholds %>% 
        filter(Model == gcm) %>% 
        filter(wl == wl_) %>% 
        pull(value) -> thres_val
      
      s %>% 
        filter(year(time) >= thres_val - 10,
               year(time) <= thres_val + 10) -> s_wl
      
      s_wl %>% 
        aggregate(FUN = function(x) sum(x >= lim_k), 
                  by = "1 year") %>% 
        aperm(c(2,3,1))
      
    }) -> l_s_wl
    
   
    
    
    # CALCULATE STATS -----------------------------------------------------------------------------
    
    print(str_glue("Calculating stats + Saving"))
    
    
    # *****************
    
    
    # WITHOUT BOOTSTRAP (ALL)
    
    l_s_wl %>% 
      {do.call(c, c(., along = "time"))} %>%
      
      st_apply(c(1,2), 
               fn_statistics, index = seq_len(dim(.)[3]),
               FUTURE = T,
               .fname = "stats") %>% 
      aperm(c(2,3,1)) -> s_wl_stats_reg_all
    
    saveRDS(s_wl_stats_reg_all, str_glue("{dir_tiles}/s_wl_stats_reg_all_{wl_}_{r}.rds"))
    
    
    # *****************
    
    
    # WITHOUT BOOTSTRAP (RegCM)
    
    l_s_wl %>% 
      .[str_which(tb_models$rcm, "RegCM")] %>%
      {do.call(c, c(., along = "time"))} %>%
      
      st_apply(c(1,2), 
               fn_statistics, index = seq_len(dim(.)[3]),
               FUTURE = T,
               .fname = "stats") %>% 
      aperm(c(2,3,1)) -> s_wl_stats_reg_regcm
    
    saveRDS(s_wl_stats_reg_regcm, str_glue("{dir_tiles}/s_wl_stats_reg_regcm_{wl_}_{r}.rds"))
    
    
    # *****************
    
    
    # WITHOUT BOOTSTRAP (REMO)
    
    l_s_wl %>% 
      .[str_which(tb_models$rcm, "REMO")] %>% 
      {do.call(c, c(., along = "time"))} %>%
      
      st_apply(c(1,2), 
               fn_statistics, index = seq_len(dim(.)[3]),
               FUTURE = T,
               .fname = "stats") %>% 
      aperm(c(2,3,1)) -> s_wl_stats_reg_remo
    
    saveRDS(s_wl_stats_reg_remo, str_glue("{dir_tiles}/s_wl_stats_reg_remo_{wl_}_{r}.rds"))
    
    
    # *****************
    
    
    # WITH BOOTSTRAP (ONLY REMO)
    
    l_s_wl %>% 
      .[str_which(tb_models$rcm, "REMO")] %>% 
      {do.call(c, c(., along = "time"))} %>%
      
      st_apply(c(1,2), 
               function(x){
                 
                 map_df(1:1000, function(...){
                   fn_statistics(x, sample(length(x), length(x), replace = T), j = T)
                 }) %>% 
                   summarize(across(everything(), ~mean(.x) %>% round())) %>% 
                   unlist()
                 
               },
               
               FUTURE = T,
               .fname = "stats") %>% 
      aperm(c(2,3,1)) -> s_wl_stats_boot_remo
    
    saveRDS(s_wl_stats_boot_remo, str_glue("{dir_tiles}/s_wl_stats_boot_remo_{wl_}_{r}.rds"))
    
    
    
  }) # end of wl loop
  
} # end of tiles loop



# MOSAICK -----------------------------------------------------------------------------------------


c("s_wl_stats_reg_all",
  "s_wl_stats_reg_regcm",
  "s_wl_stats_reg_remo",
  "s_wl_stats_boot_remo") %>% 
  
  walk(function(experiment){
    
    print(str_glue("Mosaicking {experiment}"))
    
    walk(c("1.0", "3.0"), function(wl_){
      
      dir_tiles %>% 
        list.files(full.names = T) %>% 
        str_subset(experiment) %>% 
        str_subset(str_glue("_{wl_}_")) %>% 
        
        map(readRDS) %>% 
        map(function(s){
          
          s %>% 
            st_set_dimensions("lon", st_get_dimension_values(s, "lon") %>% round(1)) %>% 
            st_set_dimensions("lat", st_get_dimension_values(s, "lat") %>% round(1)) %>%
            st_set_crs(4326)
          
        }) %>% 
        map(as, "SpatRaster") %>% 
        do.call(terra::merge, .) %>% 
        st_as_stars() %>%
        split("band") %>%
        st_set_dimensions(which = c(1,2), name = c("lon", "lat")) -> mos
      
      
      experiment %>% 
        str_sub(6) %>% 
        str_replace_all("_", "-") -> name_experiment
      
      wl_ %>% 
        str_replace("[.]", "p") -> name_wl
      
      str_glue("/mnt/bucket_mine/results/global_remo_bootstrap/days-above-32C_{name_experiment}_{name_wl}C.nc") -> name_file
      
      func_write_nc_notime(mos, 
                           name_file)
      
      
    })
    
  })


unlink(dir_tiles, recursive = T)


































hist(ts)
dbinom(0:10, 365, mean(ts)/365) %>% plot()

1-pbinom(4, 365, mean(ts/365)) # prob of having x days > thres or more 
pbinom(0:15, 365, mean(ts/365)) %>% plot() # theorical
ecdf(ts) %>% plot() # empirical

qbinom(0.95, 365, mean(ts/365)) # days > thres that happen x*100 of times
qbinom(seq(0.01,0.99,0.05), 365, mean(ts/365)) %>% plot() # theoretical
quantile(ts, seq(0.1, 0.99, 0.05)) %>% plot() # empirical


rbinom(126, 365, mean(ts)/365) %>% hist()


































