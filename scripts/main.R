
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




# MODELS TABLE ------------------------------------------------------------------------------------

tb_models <- 
  dir_raw_data %>% 
  list.files() %>% 
  str_split("_", simplify = T) %>% 
  .[,c(6,3)] %>% 
  as_tibble(.name_repair = ~c("rcm", "gcm")) %>%
  {unique(.[c("rcm", "gcm")])}


tb_models <- 
  tb_models %>% 
  mutate(calendar = case_when(str_detect(gcm, "Had") ~ 360,
                              str_detect(rcm, "REMO") & str_detect(gcm, "Had", negate = T) ~ 365.25,
                              str_detect(rcm, "RegCM") & str_detect(gcm, "MPI") ~ 365.25,
                              str_detect(rcm, "RegCM") & str_detect(gcm, "Nor") ~ 365))






# CDO PRE-PROCESS ---------------------------------------------------------------------------------

set_units(32, degC) %>% 
  set_units(K) %>% 
  drop_units() -> lim_k


dir_derived <- str_glue("{dir_pers_disk}/derived")

{
  
  # dir.create(dir_derived)
  # 
  # future_pwalk(tb_models[3:6,], function(rcm, gcm, ...){
  #   
  #   ff <- 
  #     dir_raw_data %>% 
  #     list.files() %>% 
  #     str_subset(rcm) %>% 
  #     str_subset(gcm) %>% 
  #     {str_glue("{dir_raw_data}/{.}")} %>% 
  #     str_flatten(" ")
  #   
  #   outfile <- 
  #     str_glue("{dir_derived}/days-gec-32deg_{rcm}_{gcm}.nc")
  #   
  #   str_glue("cdo -yearsum -gec,{lim_k} -cat {ff} {outfile}") %>% 
  #     system()
  #   
  # })
  
}




# READ DATA ---------------------------------------------------------------------------------------

# (read in the order of the table)

l_s <- 
  pmap(tb_models, function(rcm, gcm, ...){
    
    dir_derived %>% 
      list.files(full.names = T) %>%
      str_subset(rcm) %>% 
      str_subset(gcm) %>% 
      read_ncdf(proxy = F)
    
  }) %>% 
  map(function(s){
    
    s %>% 
      st_set_dimensions("time",
                        values = st_get_dimension_values(s, "time") %>% 
                          year()) %>% 
      drop_units() %>% 
      setNames("n_days")
  })




# TILE DATA ---------------------------------------------------------------------------------------

# f <-
#   str_glue("{dir_pers_disk}/raw_data") %>%
#   list.files(full.names = T) %>%
#   str_subset("REMO") %>% 
#   .[1]
# 
# source("scripts/tiling.R")




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





# *****

fn_statistics <- function(ts, index = seq_along(ts), j = F){
  
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






walk(c("1.0", "2.0", "3.0"), function(wl_){
  
  # wl_ <- "2.0"
  
  print(str_glue(" "))
  print(str_glue("Processing WL {wl_}"))
  
  
  
  # SLICE WL
  
  print(str_glue("Slicing"))
  
  l_s_wl <- 
    future_map2(l_s, seq(1,6), function(s, i){
      
      rcm <- tb_models$rcm[i]
      gcm <- tb_models$gcm[i]
      
      thres_val <- 
        thresholds %>% 
        filter(Model == gcm) %>% 
        filter(wl == wl_) %>% 
        pull(value)
      
      s %>% 
        filter(time >= thres_val - 10,
               time <= thres_val + 10)
      
    })
  
  
  
  # CALCULATE STATS
  
  name_wl <- 
    wl_ %>% 
    str_replace("[.]", "p")
  
  
  # WHOLE ENSEMBLE - NO BOOTSTRAP
  
  l_s_wl %>% 
    {do.call(c, c(., along = "time"))} %>%
    
    st_apply(c(1,2), 
             fn_statistics,
             FUTURE = T,
             .fname = "stats") %>% 
    aperm(c(2,3,1)) %>% 
    split("stats") -> s_result
  
  
  func_write_nc_notime(s_result,
                       str_glue("{dir_bucket_mine}/results/global_remo_bootstrap/days-above-32C_stats-reg-all_{name_wl}C.nc"))
  
  
  # *****************
  
  
  # RegCM ENSEMBLE - NO BOOTSTRAP
  
  l_s_wl %>% 
    .[str_which(tb_models$rcm, "RegCM")] %>%
    {do.call(c, c(., along = "time"))} %>%
    
    st_apply(c(1,2), 
             fn_statistics,
             FUTURE = T,
             .fname = "stats") %>% 
    aperm(c(2,3,1)) %>% 
    split("stats") -> s_result
  
  func_write_nc_notime(s_result, 
                       str_glue("{dir_bucket_mine}/results/global_remo_bootstrap/days-above-32C_stats-reg-regcm_{name_wl}C.nc"))
  
  
  
  # *****************
  
  
  # REMO ENSEMBLE - NO BOOTSTRAP
  
  l_s_wl %>% 
    .[str_which(tb_models$rcm, "REMO")] %>%
    {do.call(c, c(., along = "time"))} %>%
    
    st_apply(c(1,2), 
             fn_statistics,
             FUTURE = T,
             .fname = "stats") %>% 
    aperm(c(2,3,1)) %>% 
    split("stats") -> s_result
  
  func_write_nc_notime(s_result, 
                       str_glue("{dir_bucket_mine}/results/global_remo_bootstrap/days-above-32C_stats-reg-remo_{name_wl}C.nc"))
  
  
  
  # *****************
  
  
  # REMO ENSEMBLE - STD BOOTSTRAP
  
  l_s_wl %>% 
    .[str_which(tb_models$rcm, "REMO")] %>%
    {do.call(c, c(., along = "time"))} %>%
    
    st_apply(c(1,2), 
             function(x){
               
               map_df(1:300, function(...){
                 fn_statistics(x, sample(length(x), length(x), replace = T), j = T)
               }) %>% 
                 summarize(across(everything(), ~mean(.x) %>% round())) %>% 
                 unlist()
               
             },
             FUTURE = T,
             .fname = "stats") %>% 
    aperm(c(2,3,1)) %>% 
    split("stats") -> s_result # ~20 min 300 iterations
  
  func_write_nc_notime(s_result, 
                       str_glue("{dir_bucket_mine}/results/global_remo_bootstrap/days-above-32C_stats-bootstd-remo_{name_wl}C.nc"))
  
  
  
  # *****************
  
  
  # REMO ENSEMBLE - PARAM BOOTSTRAP
  c(360, 365.25, 365.25) %>% mean()
  
  fn_summarize <- function(x){
    
    if(all(is.na(x))){
      c(`5%` = NA, mean = NA, `95%` = NA)
    } else {
      c(mean = mean(x),
        quantile(x, c(0.05, 0.95))) %>% 
        .[c(2,1,3)]
    }
    
  }
  
  
  tictoc::tic()
  l_s_wl %>% 
    .[str_which(tb_models$rcm, "REMO")] %>%
    {do.call(c, c(., along = "time"))} %>%
    
    st_apply(c(1,2), 
             function(x){
               
               map_df(1:300, function(...){
                 fn_statistics(rbinom(length(x), 365, mean(x/365)))
               }) %>% 
                 summarize(across(everything(), ~fn_summarize(.x) %>% round())) %>% 
                 unlist()
               
             },
             FUTURE = T,
             .fname = "stats") %>% 
    aperm(c(2,3,1)) %>% 
    split("stats") -> s_result # 30 min
  tictoc::toc()
  
  func_write_nc_notime(s_result, 
                       str_glue("{dir_bucket_mine}/results/global_remo_bootstrap/days-above-32C_stats-bootparam-remo_{name_wl}C.nc"))
  
  
  
})



l_s_wl %>% 
  .[str_which(tb_models$rcm, "REMO")] %>%
  {do.call(c, c(., along = "time"))} %>% 
  .[,200, seq(20, 100, length.out = 10), ] %>% 
  as_tibble() -> tb

saveRDS(tb, "tb_tmp.rds")





foo$n_days %>% hist()
rbinom(21*3, 365, mean(foo$n_days/365)) %>% hist()




# PRE-PROCESS FOR BOOTSTRAP -----------------------------------------------------------------------

# set_units(32, degC) %>% 
#   set_units(K) %>% 
#   drop_units() -> lim_k


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






























l_s_wl[[1]][,200,20,] %>% pull() %>% as.vector() -> ts
foo$n_days -> ts

hist(ts)
dbinom(0:10, 365, mean(ts/365)) %>% plot()

1-pbinom(4, 365, mean(ts/365)) # prob of having x days > thres or more 
pbinom(seq(min(ts), max(ts), length.out = 100), 365, mean(ts/365)) %>% plot() # theorical
ecdf(ts) %>% plot() # empirical

qbinom(0.95, 365, mean(ts/365)) # days > thres that happen x*100 of times
qbinom(seq(0.01,0.99,0.05), 365, mean(ts/365)) %>% plot() # theoretical
quantile(ts, seq(0.1, 0.99, 0.05)) %>% plot() # empirical


rbinom(126, 365, mean(ts/365)) %>% hist()




tibble(d = seq(min(ts), max(ts), length.out = 100),
       p = pbinom(d, 365, mean(ts/365)),
       e = ecdf(ts)(d)) -> tb_modeled


tibble(d = ts,
       p = ecdf(ts)(ts)) -> tb_observed

ggplot() +
  geom_line(data = tb_modeled, aes(x = d, y = p)) +
  geom_line(data = tb_observed, aes(x = d, y = p))
       
       
tb_modeled %>% 
  ggplot(aes(x = d)) +
  geom_line(aes(y = e)) +
  geom_line(aes(y = p), linetype = "3232")




tibble(q = seq(0.01,0.99,0.05),
       d = qbinom(seq(0.01,0.99,0.05), 365, mean(ts/365))) %>% 
  
  ggplot()




tb %>% 
  mutate(lat = round(lat, 1)) %>% 
  filter(lat == 26.1) -> foo 

ggplot(foo, aes(n_days)) +
  stat_ecdf() +
# stat_function(fun = "dbinom", args = list(size = 365, prob = mean(foo$n_days/365)))
# facet_wrap(~lat, ncol = 2, scales = "free")














ggplot(foo, aes(n_days)) +        # Draw hist & density with count on y-axis
  geom_histogram(binwidth = 5) +
  geom_density(aes(y = ..density.. * (nrow(foo) * 5)))


seq(min(ts)-10, 
    max(ts)+10,
    length.out = 100) %>%
  round() %>% 
  dbinom(size = 365,
         prob = mean(ts/365))

dbinom(((min(ts))-extension):((max(ts))+extension), 365, mean(ts/365))

plot(density(ts, bw = 5)$y * length(ts) * 5)

range(ts) %>% diff() %>% {./8} %>% round() -> ext
tibble(d = seq((min(ts)-ext),(max(ts)+ext)),
       m = dbinom(d, 365, mean(ts/365)),
       e = density(ts, n = length(d), from = min(d), to = max(d))$y) %>%
  
  ggplot(aes(x = d)) +
  geom_line(aes(y = m)) +
  geom_line(aes(y = e), linetype = "2222")



tibble(t = ts) %>% 
  ggplot(aes(t)) + geom_density()
