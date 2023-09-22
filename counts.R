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
  st_as_sf(coords = c("longitude","latitude"),crs = 4326)
sf_counts


## -----------------------------------------------------------------------------
rnet_buffer20 = rnet_commute |> st_buffer(dist = 20)

sf_counts_selected = sf_counts[rnet_buffer20,]



## -----------------------------------------------------------------------------
sf_counts_joined = st_join(sf_counts_selected,rnet_commute,join = st_nearest_feature)
sf_counts_joined


## -----------------------------------------------------------------------------
val_app1 = sf_counts_joined |> left_join(AADF_sites,by = "siteID") |> filter(count_mean > 0)


## -----------------------------------------------------------------------------
tm_shape(rnet_commute)+
  tm_lines(col = "bicycle",lwd = "bicycle")+
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
             y = bicycle))+
  geom_point()+
  geom_smooth(method = "lm",
              formula = 'y ~ x',
              se = F)


## ----echo=FALSE---------------------------------------------------------------
# Convert README.Rmd to counts.R:  
knitr::purl("README.Rmd", "counts.R")

