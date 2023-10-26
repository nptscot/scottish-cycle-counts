## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message=FALSE,
  warning=FALSE
)


## -----------------------------------------------------------------------------
library(tidyverse)


## ----eval=FALSE---------------------------------------------------------------
#> dir.create("data_raw")
#> system("gh release download 1 --dir data_raw")


## ----eval=FALSE---------------------------------------------------------------
#> zipped_data = list.files(path = "data_raw",pattern = "\\.zip$",full.names = T)
#> zipped_data


## ----echo=FALSE, eval=FALSE---------------------------------------------------
#> tail(zipped_data, 1)


## ----eval=FALSE---------------------------------------------------------------
#> unzip(zipped_data, exdir = "data_raw")


## -----------------------------------------------------------------------------
files_csv = list.files("data_raw", pattern = "\\.csv$", full.names = TRUE)
files_csv


## -----------------------------------------------------------------------------
# library(data.table)
# counts = data.frame(data.table::rbindlist(lapply(files_csv,data.table::fread))) #DT's quick way to read the files 
counts = map_dfr(files_csv, read_csv,show_col_types = FALSE)
dim(counts)
counts


## -----------------------------------------------------------------------------
counts_monthly = counts |>
  mutate(
    year = year(endTime),
    month = month(endTime, label = TRUE),
    # date rounded to nearest month in 2020-01-01 format:
    date = lubridate::floor_date(endTime, unit = "month")
  ) |>
  group_by(date, area) |>
  summarise(
    count = sum(count)
  )
# Add column with names for most common areas:
# Most common areas are:
area_counts = counts_monthly |>
  group_by(area) |>
  summarise(
    count = sum(count)
  )

top_5_areas = area_counts |> 
  slice_max(count, n = 5)

# Add column that is area name if in top 5, else "Other":
counts_monthly_top = counts_monthly |>
  mutate(
    Area = ifelse(
      area %in% top_5_areas$area,
      area,
      "Other"
    )
  ) |>
  group_by(date, Area) |>
  summarise(
    count = sum(count)
  )


## -----------------------------------------------------------------------------
counts_monthly_top |>
  ggplot(aes(x = date, y = count, colour = Area)) +
  geom_line() +
  # Add log y-axis:
  scale_y_log10()


## -----------------------------------------------------------------------------
counts |>
  group_by(siteID,location) |> 
  unique() |> 
  count() |> 
  ggplot(aes(n))+geom_histogram(bins = 30)


## -----------------------------------------------------------------------------
range(counts$startTime)
difftime(range(counts$startTime)[2],range(counts$startTime)[1],units = "days")


## -----------------------------------------------------------------------------
repeated_sites = counts |>
  group_by(siteID) |> 
  filter(n()>378) |> 
  select(siteID) |>
  unique()

repeated_sites$siteID  |> length()

fewer_sites = counts |>
  group_by(siteID) |> 
  filter(n()<300) |> 
  select(siteID) |>
  unique()

fewer_sites$siteID  |> length()



## -----------------------------------------------------------------------------

clean_counts = counts |>
  filter(startTime < as.Date("2023-06-01")) |> 
  anti_join(repeated_sites,by =join_by(siteID)) |>
  anti_join(fewer_sites,by =join_by(siteID)) |> 
  filter(n()==365,.by = siteID)

clean_counts |>
  group_by(siteID,location) |> 
  unique() |> 
  count() |> summary()



## -----------------------------------------------------------------------------
AADF_sites = clean_counts |> 
  summarise(across(count,list(mean = mean,
                              median = median,
                              min = min,
                              max = max)),.by = siteID)

AADF_sites


## -----------------------------------------------------------------------------
counts_per_area = counts |> select(siteID,area) |> unique()

AADF_sites |>
  left_join(counts_per_area,by = join_by(siteID)) |> 
  mutate(
    area = ifelse(
      area %in% top_5_areas$area,
      area,
      "Other"
    )) |> 
  ggplot(aes(x=fct_reorder(area,count_mean,.desc = T),
             y=count_mean))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(col = area),alpha = 0.2,shape = 19,show.legend = F)+
  coord_cartesian(ylim = c(0,1000))+
  theme_minimal()+
  labs(x = "Area",
       y = "AADF")


## -----------------------------------------------------------------------------
scot_bank_holidays = as.Date(c("2022/06/03",
                               "2022/06/02",
                               "2022/08/1",
                               "2022/11/30",
                               "2022/12/25",
                               "2022/12/26",
                               "2022/12/27",
                               "2023/01/02",
                               "2023/01/03",
                               "2023/04/07",
                               "2023/05/01",
                               "2023/05/08",
                               "2023/05/29"))


## -----------------------------------------------------------------------------
ADF_dtype = clean_counts |> 
  mutate(d.type = case_when(as.Date(startTime) %in% scot_bank_holidays~"Bank Holiday",
                            wday(startTime,week_start = 1)<6~"Weekday",
                            TRUE~"Weekend")) |> 
  summarise(across(count,list(mean = mean,
                              median = median,
                              min = min,
                              max = max)),
            .by = c(siteID,d.type))

ADF_dtype


## -----------------------------------------------------------------------------
counts_per_area = counts |> select(siteID,area) |> unique()

ADF_dtype |>
  left_join(counts_per_area,by = join_by(siteID)) |> 
  mutate(
    area = ifelse(
      area %in% top_5_areas$area,
      area,
      "Other"
    )) |> 
  ggplot(aes(x=fct_reorder(area,count_mean,.desc = T),
             y=count_mean))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(col = area),alpha = 0.2,shape = 19,show.legend = F)+
  coord_cartesian(ylim = c(0,1000))+
  facet_wrap(d.type~.,ncol = 2)+
  theme_minimal()+
  labs(x = "Area",
       y = "ADF")




## -----------------------------------------------------------------------------
library(sf)
library(tmap)


## -----------------------------------------------------------------------------
rnet_commute = read_rds("../npt/outputs/rnet_commute.Rds")
rnet_commute


## -----------------------------------------------------------------------------
sf_counts = clean_counts |>
  select(siteID,latitude,longitude,provider,location) |>
  unique() |>
  left_join(AADF_sites,by = "siteID") |>
  filter(count_mean > 0) |> 
  st_as_sf(coords = c("longitude","latitude"),crs = 4326)
sf_counts


## -----------------------------------------------------------------------------
rnet_buffer20 = rnet_commute |> st_union() |>  st_buffer(dist = 20)

sf_counts_selected = sf_counts[rnet_buffer20,]


## -----------------------------------------------------------------------------
val_app1 = st_join(sf_counts_selected,rnet_commute,join = st_nearest_feature)
val_app1


## -----------------------------------------------------------------------------
tmap_mode("view")
tm_shape(rnet_commute)+
  tm_lines(col = "bicycle",lwd = "bicycle",lwd.legend = tm_legend_combine("col"))+
  tm_shape(val_app1)+
  tm_dots(col = "count_mean")


## -----------------------------------------------------------------------------
val_app1 |> 
  st_drop_geometry() |>
  mutate(ratio = bicycle/count_mean) |> 
  ggplot(aes(ratio))+
  geom_histogram()


## -----------------------------------------------------------------------------
val_app1 |>
  st_drop_geometry() |>
  ggplot(aes(x = count_mean,
             y = bicycle)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = 'y ~ x',
              se = F) +
  coord_fixed(xlim = c(0, max(
    c(val_app1$count_mean, val_app1$bicycle)
  )),
  ylim = c(0, max(
    c(val_app1$count_mean, val_app1$bicycle)
  )))


## -----------------------------------------------------------------------------
lm_app1 = lm(bicycle ~ count_mean+0,data = val_app1)
summary(lm_app1)


## -----------------------------------------------------------------------------
sel_counts_buf30 = st_buffer(sf_counts_selected,dist = 30)


## -----------------------------------------------------------------------------
counts_overlap = st_intersects(sel_counts_buf30, sel_counts_buf30)


## -----------------------------------------------------------------------------
tmap_mode("view")
tm_shape(sel_counts_buf30[sel_counts_buf30$siteID %in% c("EDH0040","EDH0041"),])+
  tm_polygons(alpha = 0.5)+
  tm_shape(sf_counts_selected)+
  tm_dots()
  


## -----------------------------------------------------------------------------
tmap_mode("view")
tm_shape(sel_counts_buf30[sel_counts_buf30$siteID %in% c("EDH0042","EDH0043","EDH0044","EDH0045"),])+
  tm_polygons(alpha=0.3)+
  tm_shape(sf_counts_selected)+
  tm_dots()


## -----------------------------------------------------------------------------
grouped_counts =
  do.call(rbind,
          lapply(unique(counts_overlap),
                 function(x) {
                   tmp_group = sf_counts_selected[x,]
                   
                   # Count aggreagation
                   simp_data = tmp_group |>
                     # Removing the direction from the location string
                     mutate(location = str_remove(location,
                                                  "\\s\\w*bound")) |>
                     st_drop_geometry() |>
                     # Extracting the first value for the siteID,
                     # provider and adds up the counts for sites with
                     # the same 'location'
                     summarise(across(c("siteID", "provider"),
                                      \(x) head(x, n = 1)),
                               across(starts_with("count_"), sum),
                               .by =  "location")
                   
                   simp_group = tmp_group |>
                     select(siteID) |>
                     filter(siteID %in% simp_data$siteID)
                   
                   simp_counts = simp_group |>
                     left_join(simp_data, by = "siteID") |>
                     relocate(location, .after = provider) |>
                     relocate(geometry, .after = count_max)
                   
                   return(simp_counts)
                 }))


## -----------------------------------------------------------------------------
grouped_counts$nearest_edge = st_nearest_feature(grouped_counts,
                                                 rnet_commute,
                                                 check_crs = T)


val_app2 = cbind(grouped_counts,st_drop_geometry(rnet_commute)[grouped_counts$nearest_edge,])


## -----------------------------------------------------------------------------
tm_shape(grouped_counts)+
  tm_dots()+
  tm_shape(rnet_commute[grouped_counts$nearest_edge,])+
  tm_lines()


## -----------------------------------------------------------------------------
val_app2 |>
  st_drop_geometry() |>
  ggplot(aes(x = count_mean,
             y = bicycle)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = 'y ~ x',
              se = F) +
  coord_fixed(xlim = c(0, max(
    c(val_app2$count_mean, val_app2$bicycle)
  )),
  ylim = c(0, max(
    c(val_app2$count_mean, val_app2$bicycle)
  )))


## -----------------------------------------------------------------------------
lm_app2 = lm(bicycle ~ count_mean+0,data = val_app2)
summary(lm_app2)


## ----eval=FALSE, warning=FALSE------------------------------------------------
#> oldwd <- setwd("../npt/outputdata")
#> system("gh release download v2023-08-18-10-42-44_commit_cbb84b024550d638dbca066c5850d1b03d55fc66 --clobber")
#> system("gh release download v2023-09-10-17-43-21.109279_commit_86ae338b12f523c27fcc290f48105f2e5dbdcab7 --clobber")
#> setwd(oldwd)


## -----------------------------------------------------------------------------
networks_files = list.files(path = "../npt/outputdata/",pattern = "^rnet.*_list",full.names = TRUE)
networks_files              


## -----------------------------------------------------------------------------
library(purrr)

rnet_val = function(rnet_path,counts){
  
  # Loading the network RDS file
  rnet_nested_list = read_rds(rnet_path)
  
  # Detecting the name of the main purpose
  main_rnet_name = str_extract(rnet_path,"(commute|school)")
  
  # Flattening the list
  rnet_flat_list = rnet_nested_list |> list_flatten()
  
  
  lst_names = paste(main_rnet_name,names(rnet_flat_list),sep = ".")
  
  # Assigning names to the lists
  names(rnet_flat_list) = lst_names
  
  lm_rnet = 
    lapply(lst_names,function(rnet_name){
    rnet = rnet_flat_list[[rnet_name]]
    
  # Creating buffer
  rnet_buffer20 = rnet |>
    st_union() |>
    st_buffer(dist = 20)
  
  # Subsetting counts based on buffer
  counts_selected = counts[rnet_buffer20,]
  
  # Creating buffer around counts 30 m
  counts_buf30 = st_buffer(counts_selected,dist = 30)
  
  # Finding overlapping counts
  counts_overlap_30 = st_intersects(counts_buf30, counts_buf30)
  
  # Processing overlaps and aggregating if possible
  aggregated_counts =
    do.call(rbind,
            lapply(unique(counts_overlap_30),
                   function(x) {
                     tmp_group = counts_selected[x, ]
                     
                     # Count aggreagation
                     simp_data = tmp_group |>
                       # Removing the direction from the location string
                       mutate(location = str_remove(location,
                                                    "\\s\\w*bound")) |>
                       st_drop_geometry() |>
                       # Extracting the first value for the siteID,
                       # provider and adds up the counts for sites with
                       # the same 'location'
                       summarise(across(c("siteID", "provider"),
                                        \(x) head(x, n = 1)),
                                 across(starts_with("count_"), sum),
                                 .by =  "location")
                     
                     simp_group = tmp_group |>
                       select(siteID) |>
                       filter(siteID %in% simp_data$siteID)
                     
                     simp_counts = simp_group |>
                       left_join(simp_data, by = "siteID") |>
                       relocate(location, .after = provider) |>
                       relocate(geometry, .after = count_max)
                     
                     return(simp_counts)
                   }))
  
  aggregated_counts$nearest_edge = st_nearest_feature(aggregated_counts,
                                                      rnet,
                                                      check_crs = T)
  
  val_counts = cbind(aggregated_counts,
                     st_drop_geometry(rnet)[aggregated_counts$nearest_edge,])
  
  
  
  
  lm_counts = lm(bicycle ~ count_mean + 0,data = val_counts)
  
  return(lm_counts)
    
  })
  names(lm_rnet) = lst_names
  
  return(lm_rnet)
}


## ----eval=FALSE---------------------------------------------------------------
#> val_results = lapply(networks_files,
#>                      rnet_val,
#>                      counts = sf_counts
#>                      )


## ----include=FALSE,eval=FALSE-------------------------------------------------
#> dir.create("interim_results")
#> saveRDS(val_results,
#>         file = "interim_results/val_results.Rds")


## ----include=FALSE------------------------------------------------------------
val_results = readRDS("interim_results/val_results.Rds")


## -----------------------------------------------------------------------------
library(broom)
val_results_flat = val_results |> list_flatten(name_spec = "{inner}")

tbl_val_results = tibble(network = names(val_results_flat),
                         intercept = 0,
                         coef = vapply(val_results_flat,coef,FUN.VALUE = 0),
                         R2 = vapply(val_results_flat,\(x) summary(x)$adj.r.squared,FUN.VALUE = 0)
       )

tbl_val_results


## -----------------------------------------------------------------------------
rm(val_results,val_results_flat)


## -----------------------------------------------------------------------------
library(GWmodel)

rnet_val_gwr = function(rnet_path, counts) {
  # Loading the network RDS file
  rnet_nested_list = read_rds(rnet_path)
  
  # Detecting the name of the main purpose
  main_rnet_name = str_extract(rnet_path, "(commute|school)")
  
  # Flattening the list
  rnet_flat_list = rnet_nested_list |> list_flatten()
  
  
  lst_names = paste(main_rnet_name, names(rnet_flat_list), sep = ".")
  
  # Assigning names to the lists
  names(rnet_flat_list) = lst_names
  
  reg_rnet =
    lapply(lst_names, function(rnet_name) {
      rnet = rnet_flat_list[[rnet_name]]
      
      # Creating buffer
      rnet_buffer20 = rnet |>
        st_union() |>
        st_buffer(dist = 20)
      
      # Subsetting counts based on buffer
      counts_selected = counts[rnet_buffer20, ]
      
      # Creating buffer around counts 30 m
      counts_buf30 = st_buffer(counts_selected, dist = 30)
      
      # Finding overlapping counts
      counts_overlap_30 = st_intersects(counts_buf30, counts_buf30)
      
      # Processing overlaps and aggregating if possible
      aggregated_counts =
        do.call(rbind,
                lapply(unique(counts_overlap_30),
                       function(x) {
                         tmp_group = counts_selected[x,]
                         
                         # Count aggreagation
                         simp_data = tmp_group |>
                           # Removing the direction from the location string
                           mutate(location = str_remove(location,
                                                        "\\s\\w*bound")) |>
                           st_drop_geometry() |>
                           # Extracting the first value for the siteID,
                           # provider and adds up the counts for sites with
                           # the same 'location'
                           summarise(across(c("siteID", "provider"),
                                            \(x) head(x, n = 1)),
                                     across(starts_with("count_"), sum),
                                     .by =  "location")
                         
                         simp_group = tmp_group |>
                           select(siteID) |>
                           filter(siteID %in% simp_data$siteID)
                         
                         simp_counts = simp_group |>
                           left_join(simp_data, by = "siteID") |>
                           relocate(location, .after = provider) |>
                           relocate(geometry, .after = count_max)
                         
                         return(simp_counts)
                       }))
      
      aggregated_counts$nearest_edge = st_nearest_feature(aggregated_counts,
                                                          rnet,
                                                          check_crs = T)
      
      val_counts = cbind(aggregated_counts,
                         st_drop_geometry(rnet)[aggregated_counts$nearest_edge, ])
      
      
      
      
      bw <- bw.gwr(
        formula = bicycle ~ count_mean,
        approach = "AIC",
        adaptive = T,
        data = as_Spatial(val_counts)
      )
      gwr.mod <- gwr.basic(
        formula = bicycle ~ count_mean,
        adaptive = T,
        data = as_Spatial(val_counts),
        bw = bw
      )
      
      return(gwr.mod)
      
    })
  names(reg_rnet) = lst_names
  
  return(reg_rnet)
}


## ----eval=FALSE,echo=FALSE----------------------------------------------------
#> val_results_gwr = lapply(networks_files,
#>                      rnet_val_gwr,
#>                      counts = sf_counts
#>                      )
#> 


## ----include=FALSE,eval=FALSE-------------------------------------------------
#> saveRDS(val_results_gwr,
#>         file = "interim_results/val_results_gwr.Rds")


## ----include=FALSE------------------------------------------------------------
val_results_gwr = readRDS("interim_results/val_results_gwr.Rds")


## -----------------------------------------------------------------------------
library(broom)
val_results_gwr_flat = val_results_gwr |>
  list_flatten(name_spec = "{inner}")

tbl_val_gwr_results = tibble(network = names(val_results_gwr_flat),
       mean_intercept = vapply(val_results_gwr_flat,\(x) mean(x$SDF@data[, 1]),FUN.VALUE = 0),
       median_intercept = vapply(val_results_gwr_flat,\(x) median(x$SDF@data[, 1]),FUN.VALUE = 0),
       mean_coef = vapply(val_results_gwr_flat,\(x) mean(x$SDF@data[, 2]),FUN.VALUE = 0),
       median_coef = vapply(val_results_gwr_flat,\(x) median(x$SDF@data[, 2]),FUN.VALUE = 0),
       bw = vapply(val_results_gwr_flat,\(x) x$GW.arguments$bw,FUN.VALUE = 0),
       R2 = vapply(val_results_gwr_flat,\(x) x$GW.diagnostic$gw.R2,FUN.VALUE = 0)
       )

tbl_val_gwr_results


## -----------------------------------------------------------------------------
library(kableExtra)
library(knitr)
tbl_val_gwr_results |>
  left_join(tbl_val_results,by="network",suffix = c("",".")) |> 
  separate_wider_delim(network,delim = ".",names = c("Purpose","Type")) |> 
  kbl(digits=4) |>
  kable_classic_2("hover", full_width = T) |>
  add_header_above(c(" " = 2, "GWR" = 6, "Global" = 3)) |> 
  column_spec(1, bold = T) |> 
  collapse_rows(columns = 1:2, valign = "top") |>
  as_image(width = 7,file = "README_files/figure-gfm/model_table.png")


## ----eval=FALSE---------------------------------------------------------------
#> oldwd <- setwd("../npt/outputdata")
#> system("gh release download v2023-10-18-08-46-09.827274_commit_57d2dcdd739f665bb8d3fb5540cb6bbd10fcb08f -p combined_network_tile.zip")
#> unzip("combined_network_tile.zip")
#> setwd(oldwd)


## -----------------------------------------------------------------------------
combined_network = st_read(dsn = "../npt/outputdata/combined_network_tile.geojson")
combined_network


## -----------------------------------------------------------------------------
tmap_mode("plot")
qtm(combined_network)


## -----------------------------------------------------------------------------
cnet_val_gwr = function(cnet, counts) {
  lst_names = names(cnet)[grepl("_", names(cnet))]
  
  # Creating buffer
  rnet_buffer20 = cnet |>
    st_union() |>
    st_buffer(dist = 20)
  
  # Subsetting counts based on buffer
  counts_selected = counts[rnet_buffer20,]
  
  # Creating buffer around counts 30 m
  counts_buf30 = st_buffer(counts_selected, dist = 30)
  
  # Finding overlapping counts
  counts_overlap_30 = st_intersects(counts_buf30, counts_buf30)
  
  # Processing overlaps and aggregating if possible
  aggregated_counts =
    do.call(rbind,
            lapply(unique(counts_overlap_30),
                   function(x) {
                     tmp_group = counts_selected[x, ]
                     
                     # Count aggreagation
                     simp_data = tmp_group |>
                       # Removing the direction from the location string
                       mutate(location = str_remove(location,
                                                    "\\s\\w*bound")) |>
                       st_drop_geometry() |>
                       # Extracting the first value for the siteID,
                       # provider and adds up the counts for sites with
                       # the same 'location'
                       summarise(across(c("siteID", "provider"),
                                        \(x) head(x, n = 1)),
                                 across(starts_with("count_"), sum),
                                 .by =  "location")
                     
                     simp_group = tmp_group |>
                       select(siteID) |>
                       filter(siteID %in% simp_data$siteID)
                     
                     simp_counts = simp_group |>
                       left_join(simp_data, by = "siteID") |>
                       relocate(location, .after = provider) |>
                       relocate(geometry, .after = count_max)
                     
                     return(simp_counts)
                   }))
  
  aggregated_counts$nearest_edge = st_nearest_feature(aggregated_counts,
                                                      cnet,
                                                      check_crs = T)
  
  
  val_counts = cbind(aggregated_counts,
                     st_drop_geometry(cnet)[aggregated_counts$nearest_edge,])
  
  rnet_name = lst_names[1]
  
  val_counts = cbind(aggregated_counts,
                     st_drop_geometry(cnet)[aggregated_counts$nearest_edge, lst_names]) |> as_Spatial()
  
  
  reg_rnet =
    lapply(lst_names, function(rnet_name) {
      
      form = as.formula(paste(rnet_name, "~ count_mean"))
      
      bw <- bw.gwr(
        formula = form,
        approach = "AIC",
        adaptive = T,
        data = val_counts
      )
      gwr.mod <- gwr.basic(
        formula = form,
        adaptive = T,
        data = val_counts,
        bw = bw
      )
      
      return(gwr.mod)
      
    })
  names(reg_rnet) = lst_names
  
  return(reg_rnet)
}


## ----eval=FALSE---------------------------------------------------------------
#> val_results_gwr_last = cnet_val_gwr(cnet = combined_network,counts = sf_counts)


## ----eval=FALSE,include=FALSE-------------------------------------------------
#> saveRDS(val_results_gwr_last,
#>         file = "interim_results/val_results_gwr_last.Rds")


## ----include=FALSE------------------------------------------------------------
val_results_gwr_last = read_rds("interim_results/val_results_gwr_last.Rds")


## -----------------------------------------------------------------------------
tbl_val_gwr_results_last = tibble(network = names(val_results_gwr_last),
                                  mean_intercept = vapply(val_results_gwr_last,\(x) mean(x$SDF@data[, 1]),FUN.VALUE = 0),
                                  median_intercept = vapply(val_results_gwr_last,\(x) median(x$SDF@data[, 1]),FUN.VALUE = 0),
                                  mean_coef = vapply(val_results_gwr_last,\(x) mean(x$SDF@data[, 2]),FUN.VALUE = 0),
                                  median_coef = vapply(val_results_gwr_last,\(x) median(x$SDF@data[, 2]),FUN.VALUE = 0),
                                  bw = vapply(val_results_gwr_last,\(x) x$GW.arguments$bw,FUN.VALUE = 0),
                                  R2 = vapply(val_results_gwr_last,\(x) x$GW.diagnostic$gw.R2,FUN.VALUE = 0)
       )

tbl_val_gwr_results_last |>
  kbl(digits=4) |>
  kable_classic_2("hover", full_width = T) |>
  as_image(width = 7,file = "README_files/figure-gfm/val_last_table.png")


## ----echo=FALSE---------------------------------------------------------------
# Convert README.Rmd to counts.R:  
knitr::purl("README.Rmd", "counts.R")

