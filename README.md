
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scottish-cycle-counts

<!-- badges: start -->
<!-- badges: end -->

The goal of Scottish-cycle-counts is to read-in a process data on
cycling volumes in Scotland.

``` r
library(tidyverse)
```

First of all, a copy of the original files can be obtained with the code
below. This code is to be run just once.

``` r
dir.create("data_raw")
system("gh release download 1 --dir data_raw")
```

The input dataset is a single .zip file:

``` r
zipped_data = list.files(path = "data_raw",pattern = "\\.zip$",full.names = T)
zipped_data
```

We can unzip it as follows:

``` r
unzip(zipped_data, exdir = "data_raw")
```

``` r
files_csv = list.files("data_raw", pattern = "\\.csv$", full.names = TRUE)
files_csv
#>  [1] "data_raw/Aberdeen City Council.csv"                            
#>  [2] "data_raw/Aberdeenshire Council.csv"                            
#>  [3] "data_raw/City of Edinburgh Council.csv"                        
#>  [4] "data_raw/Comhairle nan Eilean Siar (Western Isles Council).csv"
#>  [5] "data_raw/East Ayrshire Council.csv"                            
#>  [6] "data_raw/East Dunbartonshire Council.csv"                      
#>  [7] "data_raw/Glasgow City Council.csv"                             
#>  [8] "data_raw/John Muir Way.csv"                                    
#>  [9] "data_raw/National Cycle Nework (Scotland) - Sustrans.csv"      
#> [10] "data_raw/National Monitoring Framework - Cycling Scotland.csv" 
#> [11] "data_raw/North Ayrshire Council.csv"                           
#> [12] "data_raw/Perth and Kinross Council.csv"                        
#> [13] "data_raw/Scotland North East Trunk Roads.csv"                  
#> [14] "data_raw/Scotland North West Trunk Roads.csv"                  
#> [15] "data_raw/Scotland South East Trunk Roads.csv"                  
#> [16] "data_raw/Scotland South West Trunk Roads.csv"                  
#> [17] "data_raw/South Ayrshire Council.csv"                           
#> [18] "data_raw/South Lanarkshire Council.csv"                        
#> [19] "data_raw/Stirling Council.csv"                                 
#> [20] "data_raw/The Highland Council.csv"
```

We can read this file in R as follows:

``` r
# library(data.table)
# counts = data.frame(data.table::rbindlist(lapply(files_csv,data.table::fread))) #DT's quick way to read the files 
counts = map_dfr(files_csv, read_csv,show_col_types = FALSE)
dim(counts)
#> [1] 230307     10
counts
#> # A tibble: 230,307 × 10
#>    area    count endTime             latitude location longitude provider siteID
#>    <chr>   <dbl> <dttm>                 <dbl> <chr>        <dbl> <chr>    <chr> 
#>  1 Aberde…    61 2023-06-14 00:00:00     57.1 A956 We…     -2.09 Aberdee… ABE65…
#>  2 Aberde…    26 2023-06-14 00:00:00     57.2 Dyce Dr…     -2.20 Aberdee… ABE918
#>  3 Aberde…    18 2023-06-14 00:00:00     57.2 Dyce Dr…     -2.22 Aberdee… ABE397
#>  4 Aberde…   158 2023-06-14 00:00:00     57.2 Ellon R…     -2.09 Aberdee… ABE920
#>  5 Aberde…   157 2023-06-14 00:00:00     57.1 Shell C…     -2.10 Aberdee… ABE229
#>  6 Aberde…   132 2023-06-14 00:00:00     57.2 Tillydr…     -2.11 Aberdee… ABE899
#>  7 Aberde…   395 2023-06-14 00:00:00     57.1 Deeside…     -2.10 Aberdee… ABE235
#>  8 Aberde…    20 2023-06-14 00:00:00     57.2 Farburn…     -2.17 Aberdee… ABE745
#>  9 Aberde…    91 2023-06-14 00:00:00     57.2 Parkway…     -2.10 Aberdee… ABE894
#> 10 Aberde…    39 2023-06-14 00:00:00     57.2 Site B …     -2.16 Aberdee… ABE288
#> # ℹ 230,297 more rows
#> # ℹ 2 more variables: startTime <dttm>, usmart_id <chr>
```

``` r
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
```

``` r
counts_monthly_top |>
  ggplot(aes(x = date, y = count, colour = Area)) +
  geom_line() +
  # Add log y-axis:
  scale_y_log10()
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

A quick look at the data to check if there are sites with missing data
or duplicated records:

``` r
counts |>
  group_by(siteID,location) |> 
  unique() |> 
  count() |> 
  ggplot(aes(n))+geom_histogram(bins = 30)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
range(counts$startTime)
#> [1] "2022-06-01 UTC" "2023-06-13 UTC"
difftime(range(counts$startTime)[2],range(counts$startTime)[1],units = "days")
#> Time difference of 377 days
```

Each site should have a maximum of 378 days in the dataset. The
following code detects the sites with some type of duplication and the
ones with fewer records.

``` r
repeated_sites = counts |>
  group_by(siteID) |> 
  filter(n()>378) |> 
  select(siteID) |>
  unique()

repeated_sites$siteID  |> length()
#> [1] 22

fewer_sites = counts |>
  group_by(siteID) |> 
  filter(n()<300) |> 
  select(siteID) |>
  unique()

fewer_sites$siteID  |> length()
#> [1] 133
```

A subset of the clean sites is produced, so we can do some AADT
analysis. Records after 2023-06-01 are filtered out to have only one
year of data for each site

``` r

clean_counts = counts |>
  filter(startTime < as.Date("2023-06-01")) |> 
  anti_join(repeated_sites,by =join_by(siteID)) |>
  anti_join(fewer_sites,by =join_by(siteID)) |> 
  filter(n()==365,.by = siteID)

clean_counts |>
  group_by(siteID,location) |> 
  unique() |> 
  count() |> summary()
#>     siteID            location               n      
#>  Length:383         Length:383         Min.   :365  
#>  Class :character   Class :character   1st Qu.:365  
#>  Mode  :character   Mode  :character   Median :365  
#>                                        Mean   :365  
#>                                        3rd Qu.:365  
#>                                        Max.   :365
```

Here we calculate some statistics for the whole year including mean
(Average Annual Daily Flow), median daily flow, minimum and maximum
daily flows,

``` r
AADF_sites = clean_counts |> 
  summarise(across(count,list(mean = mean,
                              median = median,
                              min = min,
                              max = max)),.by = siteID)

AADF_sites
#> # A tibble: 383 × 5
#>    siteID  count_mean count_median count_min count_max
#>    <chr>        <dbl>        <dbl>     <dbl>     <dbl>
#>  1 ABE6535       19.9           18         1        50
#>  2 ABE896        44.1           43         0       195
#>  3 ABE6545       34.3           32         2        84
#>  4 ABE895        86.0           74         0       288
#>  5 ABE918        15.6           15         0        39
#>  6 ABE122        16.6           16         0        57
#>  7 ABE123        13.8           12         0        49
#>  8 ABE400        43.5           33         2       150
#>  9 ABE288        12.5            9         0        81
#> 10 ABE893        43.3           43         0        93
#> # ℹ 373 more rows
```

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

We create a vector to store the bank holidays in Scotland extracted from
the [mygov.scot web](https://www.mygov.scot/scotland-bank-holidays)

``` r
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
```

We can calculate the same summary statistics by type of day: bank
holidays, weekends and weekdays (AAWDF).

``` r
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
#> # A tibble: 1,149 × 6
#>    siteID  d.type  count_mean count_median count_min count_max
#>    <chr>   <chr>        <dbl>        <dbl>     <dbl>     <dbl>
#>  1 ABE6535 Weekday       23.6           22         6        50
#>  2 ABE896  Weekday       46.4           46         2       101
#>  3 ABE6545 Weekday       40.3           40         2        84
#>  4 ABE895  Weekday       78.1           65         1       247
#>  5 ABE918  Weekday       19.5           20         0        39
#>  6 ABE122  Weekday       17             16         0        50
#>  7 ABE123  Weekday       12.7           11         0        45
#>  8 ABE400  Weekday       40.1           31         2       110
#>  9 ABE288  Weekday       11.8            8         0        62
#> 10 ABE893  Weekday       52.6           52         6        93
#> # ℹ 1,139 more rows
```

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

## Spatial Analysis

``` r
library(sf)
library(tmap)
```

The following code reads the network with the estimated commute trips
from the npt repository. Each edge/link of the network has four
attributes. We will focus on `bicycle` which is the estimated number of
daily commute trips in both directions.

``` r
rnet_commute = read_rds("../npt/outputs/rnet_commute.Rds")
rnet_commute
#> Simple feature collection with 1709 features and 4 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -3.33653 ymin: 55.89548 xmax: -3.12415 ymax: 55.98413
#> Geodetic CRS:  WGS 84
#> First 10 features:
#>    bicycle bicycle_go_dutch Gradient Quietness                       geometry
#> 1        0                1        0        40 LINESTRING (-3.1547 55.9749...
#> 2        0                1        0        60 LINESTRING (-3.27797 55.911...
#> 3        0                1        2        70 LINESTRING (-3.27672 55.910...
#> 4        0                1        4        70 LINESTRING (-3.14531 55.969...
#> 5        0                1        0        80 LINESTRING (-3.27741 55.911...
#> 6        0                1        0        80 LINESTRING (-3.15282 55.976...
#> 7        0                1        1        80 LINESTRING (-3.17113 55.961...
#> 8        0                1        0        90 LINESTRING (-3.1547 55.9749...
#> 9        0                1        1        90 LINESTRING (-3.16781 55.961...
#> 10       0                1        1        90 LINESTRING (-3.16555 55.962...
```

An `sf` object is created from the `clean_counts` data frame. AADF for
each counts are joined using the `siteID`

``` r
sf_counts = clean_counts |>
  select(siteID,latitude,longitude,provider,location) |>
  unique() |>
  left_join(AADF_sites,by = "siteID") |>
  filter(count_mean > 0) |> 
  st_as_sf(coords = c("longitude","latitude"),crs = 4326)
sf_counts
#> Simple feature collection with 338 features and 7 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -7.307251 ymin: 54.91291 xmax: -1.15021 ymax: 60.1511
#> Geodetic CRS:  WGS 84
#> # A tibble: 338 × 8
#>    siteID  provider         location count_mean count_median count_min count_max
#>  * <chr>   <chr>            <chr>         <dbl>        <dbl>     <dbl>     <dbl>
#>  1 ABE6535 Aberdeen City C… A96 Auc…       19.9           18         1        50
#>  2 ABE896  Aberdeen City C… F&B Way…       44.1           43         0       195
#>  3 ABE6545 Aberdeen City C… A9119 Q…       34.3           32         2        84
#>  4 ABE895  Aberdeen City C… Beach E…       86.0           74         0       288
#>  5 ABE918  Aberdeen City C… Dyce Dr…       15.6           15         0        39
#>  6 ABE122  Aberdeen City C… Maidenc…       16.6           16         0        57
#>  7 ABE123  Aberdeen City C… Maidenc…       13.8           12         0        49
#>  8 ABE400  Aberdeen City C… Seaton …       43.5           33         2       150
#>  9 ABE288  Aberdeen City C… Site B …       12.5            9         0        81
#> 10 ABE893  Aberdeen City C… Welling…       43.3           43         0        93
#> # ℹ 328 more rows
#> # ℹ 1 more variable: geometry <POINT [°]>
```

A subset of the counts are taken based on a buffer of the `rnet_commute`
object.

``` r
rnet_buffer20 = rnet_commute |> st_union() |>  st_buffer(dist = 20)

sf_counts_selected = sf_counts[rnet_buffer20,]
```

### Approach A

The nearest feature is joined to each point location

``` r
val_app1 = st_join(sf_counts_selected,rnet_commute,join = st_nearest_feature)
val_app1
#> Simple feature collection with 21 features and 11 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -3.24015 ymin: 55.92337 xmax: -3.12779 ymax: 55.9566
#> Geodetic CRS:  WGS 84
#> # A tibble: 21 × 12
#>    siteID  provider         location count_mean count_median count_min count_max
#>  * <chr>   <chr>            <chr>         <dbl>        <dbl>     <dbl>     <dbl>
#>  1 EDH0030 City of Edinbur… Dalry R…      1009.          331        40     27670
#>  2 EDH0037 City of Edinbur… Whiteho…       437.          457        54       753
#>  3 EDH0045 City of Edinbur… Melvill…       128.          128        13       266
#>  4 EDH0022 City of Edinbur… North M…       885.          888        99      1839
#>  5 EDH0035 City of Edinbur… A90 Dea…       415.          361        36       857
#>  6 EDH0039 City of Edinbur… Bruntsf…       218.          236         0       481
#>  7 EDH0043 City of Edinbur… Melvill…      1099.         1109       242      1813
#>  8 EDH0041 City of Edinbur… Mayfiel…       152.            0         0       683
#>  9 EDH0029 City of Edinbur… A8 Cors…       100.            0         0       829
#> 10 EDH0044 City of Edinbur… Melvill…       251.          257        32       446
#> # ℹ 11 more rows
#> # ℹ 5 more variables: geometry <POINT [°]>, bicycle <dbl>,
#> #   bicycle_go_dutch <dbl>, Gradient <dbl>, Quietness <dbl>
```

``` r
tmap_mode("view")
tm_shape(rnet_commute)+
  tm_lines(col = "bicycle",lwd = "bicycle",lwd.legend = tm_legend_combine("col"))+
  tm_shape(val_app1)+
  tm_dots(col = "count_mean")
```

The following plot compares the observed counts and the paired estimate
flows.

``` r
val_app1 |> 
  st_drop_geometry() |>
  mutate(ratio = bicycle/count_mean) |> 
  ggplot(aes(ratio))+
  geom_histogram()
```

![](README_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

A linear regression is used to evaluate the fit of the estimates, it is
assumed that the proportion of commute trips is constant across counts
(intercept = `0`).

``` r
lm_app1 = lm(bicycle ~ count_mean+0,data = val_app1)
summary(lm_app1)
#> 
#> Call:
#> lm(formula = bicycle ~ count_mean + 0, data = val_app1)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -4.8149 -0.5083  0.6479  4.4314  8.0276 
#> 
#> Coefficients:
#>            Estimate Std. Error t value Pr(>|t|)    
#> count_mean 0.006192   0.001606   3.855 0.000988 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 4.054 on 20 degrees of freedom
#> Multiple R-squared:  0.4263, Adjusted R-squared:  0.3976 
#> F-statistic: 14.86 on 1 and 20 DF,  p-value: 0.0009876
```

### Approach B

The previous approach assigned a count site to each road link based. So
far, it has not been addressed the fact that counts might be reporting
uni-directional flows along specific links.

Using the `sf_counts_selected` object, we produce a buffer of 30 metres
from each count site.

``` r
sel_counts_buf30 = st_buffer(sf_counts_selected,dist = 30)
```

All overlaps within the `sel_counts_buf30` are identified using the
following code.

``` r
counts_overlap = st_intersects(sel_counts_buf30, sel_counts_buf30)
```

Some overlaps might reveal count sites reporting flows of different
directions on the same edge/link. For example, the sites EDH0040 and
EDH0041 which are located on Mayfield Road.

``` r
tmap_mode("view")
tm_shape(sel_counts_buf30[sel_counts_buf30$siteID %in% c("EDH0040","EDH0041"),])+
  tm_polygons(alpha = 0.5)+
  tm_shape(sf_counts_selected)+
  tm_dots()
  
```

A more complex instance is the overlap of sensors on Melville Dr, two of
the sensors report flows on the main road and the other has data of the
adjacend path.

``` r
tmap_mode("view")
tm_shape(sel_counts_buf30[sel_counts_buf30$siteID %in% c("EDH0042","EDH0043","EDH0044","EDH0045"),])+
  tm_polygons(alpha=0.3)+
  tm_shape(sf_counts_selected)+
  tm_dots()
```

The following code aggregates some of the overlapping counts using the
`location` attribute as aggregation criteria.

``` r
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
```

As in the previous approach, the count sites are joined to the nearest
feature in the `rnet_commute` network.

``` r
grouped_counts$nearest_edge = st_nearest_feature(grouped_counts,
                                                 rnet_commute,
                                                 check_crs = T)


val_app2 = cbind(grouped_counts,st_drop_geometry(rnet_commute)[grouped_counts$nearest_edge,])
```

The following code shows the counts with the corresponding network
edges/links.

``` r
tm_shape(grouped_counts)+
  tm_dots()+
  tm_shape(rnet_commute[grouped_counts$nearest_edge,])+
  tm_lines()
```

The figure below compares the counts and the estimated flows for the
current approach

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

As in the previous approach, a linear regression is used to have a
high-level assessment of the estimations. With this approach, although
there is not a significant change in the estimate, there is a slight
improvement in the R<sup>2</sup>.

``` r
lm_app2 = lm(bicycle ~ count_mean+0,data = val_app2)
summary(lm_app2)
#> 
#> Call:
#> lm(formula = bicycle ~ count_mean + 0, data = val_app2)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -4.6021 -0.7120  0.2573  4.3372  8.0632 
#> 
#> Coefficients:
#>            Estimate Std. Error t value Pr(>|t|)   
#> count_mean 0.005965   0.001523   3.917  0.00111 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 3.916 on 17 degrees of freedom
#> Multiple R-squared:  0.4744, Adjusted R-squared:  0.4435 
#> F-statistic: 15.35 on 1 and 17 DF,  p-value: 0.001109
```

## Validation of the complete network builds

### Loading the data

Downloading the data from latest builds. RDS objects with the list of
networks have been selected for this analysis

``` r
oldwd <- setwd("../npt/outputdata")
system("gh release download v2023-08-18-10-42-44_commit_cbb84b024550d638dbca066c5850d1b03d55fc66 --clobber")
system("gh release download v2023-09-10-17-43-21.109279_commit_86ae338b12f523c27fcc290f48105f2e5dbdcab7 --clobber")
setwd(oldwd)
```

We select the ones that contain the naming pattern with a regular
expression.

``` r
networks_files = list.files(path = "../npt/outputdata/",pattern = "^rnet.*_list",full.names = TRUE)
networks_files              
#> [1] "../npt/outputdata/rnet_commute_list.Rds"
#> [2] "../npt/outputdata/rnet_school_list.Rds"
```

### Global Linear Regression Model

The following function is prepared to run the analysis as per the second
approach. Note that the linear model used for the evaluation does not
include an intercept.

``` r
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
```

The following code runs the validation for each network

``` r
val_results = lapply(networks_files,
                     rnet_val,
                     counts = sf_counts
                     )
```

Summarising the high-level results for all networks:

``` r
library(broom)
val_results_flat = val_results |> list_flatten(name_spec = "{inner}")

tbl_val_results = tibble(network = names(val_results_flat),
                         intercept = 0,
                         coef = vapply(val_results_flat,coef,FUN.VALUE = 0),
                         R2 = vapply(val_results_flat,\(x) summary(x)$adj.r.squared,FUN.VALUE = 0)
       )

tbl_val_results
#> # A tibble: 12 × 4
#>    network                   intercept    coef     R2
#>    <chr>                         <dbl>   <dbl>  <dbl>
#>  1 commute.fastest                   0 0.0626  0.0912
#>  2 commute.balanced                  0 0.0895  0.0954
#>  3 commute.quietest                  0 0.0953  0.0751
#>  4 commute.ebike                     0 0.0553  0.0814
#>  5 school.Primary_fastest            0 0.00229 0.0643
#>  6 school.Primary_balanced           0 0.00217 0.0410
#>  7 school.Primary_quietest           0 0.00270 0.0459
#>  8 school.Primary_ebike              0 0.00185 0.0468
#>  9 school.Secondary_fastest          0 0.00220 0.0274
#> 10 school.Secondary_balanced         0 0.00263 0.0321
#> 11 school.Secondary_quietest         0 0.00250 0.0227
#> 12 school.Secondary_ebike            0 0.00202 0.0227
```

``` r
rm(val_results,val_results_flat)
```

### Geographical Weighted Regression

Given that the proportion of trips by purpose might vary spatially,
e.g. urban areas might have a higher proportion of commuting trips. A
weighted geographical model is implemented. The default bi-square
function is used to determine the kernel bandwidth.

``` r
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
```

As in the previous attempt, the following code runs the validation using
a GWR model.

Summarising the high-level results for all networks for the GWR model:

``` r
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
#> # A tibble: 12 × 7
#>    network    mean_intercept median_intercept mean_coef median_coef    bw     R2
#>    <chr>               <dbl>            <dbl>     <dbl>       <dbl> <dbl>  <dbl>
#>  1 commute.f…          14.3             8.54    5.60e-2  0.0122       140 0.223 
#>  2 commute.b…          21.8            17.6     7.87e-2  0.0259       150 0.204 
#>  3 commute.q…          26.4            20.6     8.44e-2  0.0386       150 0.166 
#>  4 commute.e…          15.2             8.70    4.49e-2  0.00896      140 0.209 
#>  5 school.Pr…           1.44            1.47    3.78e-3  0.00148       40 0.258 
#>  6 school.Pr…           1.98            1.89   -1.41e-4 -0.000296     136 0.0590
#>  7 school.Pr…           2.30            1.94   -9.89e-5  0.000105     135 0.0851
#>  8 school.Pr…           1.44            1.43    2.05e-4  0.00000476   163 0.0677
#>  9 school.Se…           1.44            0.772  -4.50e-3  0.000617      31 0.548 
#> 10 school.Se…           1.59            0.966   1.48e-3  0.00142       31 0.510 
#> 11 school.Se…           1.59            1.43    2.80e-4  0.000265     125 0.136 
#> 12 school.Se…           1.37            1.05   -2.48e-3  0.000511     118 0.144
```

### Results

The following table compares the results of the two approaches.

``` r
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
```

<img src="README_files/figure-gfm/model_table.png" width="672" />

## Release 18 October 2023

### Downloading release

``` r
oldwd <- setwd("../npt/outputdata")
system("gh release download v2023-10-18-08-46-09.827274_commit_57d2dcdd739f665bb8d3fb5540cb6bbd10fcb08f -p combined_network_tile.zip")
unzip("combined_network_tile.zip")
setwd(oldwd)
```

### Loading data

``` r
combined_network = st_read(dsn = "../npt/outputdata/combined_network_tile.geojson")
#> Reading layer `combined_network_tile' from data source 
#>   `C:\Users\ts18jpf\OneDrive - University of Leeds\03_PhD\00_Misc_projects\npt\outputdata\combined_network_tile.geojson' 
#>   using driver `GeoJSON'
#> Simple feature collection with 563466 features and 26 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -7.54343 ymin: 54.67982 xmax: -0.8615 ymax: 60.76383
#> Geodetic CRS:  WGS 84
combined_network
#> Simple feature collection with 563466 features and 26 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -7.54343 ymin: 54.67982 xmax: -0.8615 ymax: 60.76383
#> Geodetic CRS:  WGS 84
#> First 10 features:
#>    all_fastest_bicycle all_fastest_bicycle_ebike all_fastest_bicycle_go_dutch
#> 1                    0                        14                            0
#> 2                    0                         0                            0
#> 3                    0                         3                            0
#> 4                    0                         1                            0
#> 5                    0                         0                            0
#> 6                    0                         0                            0
#> 7                    0                         0                            0
#> 8                    0                         0                            0
#> 9                    0                         0                            0
#> 10                   0                         0                            0
#>    all_quietest_bicycle all_quietest_bicycle_ebike
#> 1                     0                          0
#> 2                     0                          0
#> 3                     0                          0
#> 4                     0                          0
#> 5                     0                          0
#> 6                     0                          0
#> 7                     0                          0
#> 8                     0                          0
#> 9                     0                          0
#> 10                    0                          0
#>    all_quietest_bicycle_go_dutch commute_fastest_bicycle
#> 1                              0                       0
#> 2                              0                       0
#> 3                              0                       0
#> 4                              0                       0
#> 5                              0                       0
#> 6                              0                       0
#> 7                              0                       0
#> 8                              0                       0
#> 9                              0                       0
#> 10                             0                       0
#>    commute_fastest_bicycle_ebike commute_fastest_bicycle_go_dutch
#> 1                              0                                0
#> 2                              0                                0
#> 3                              0                                0
#> 4                              0                                0
#> 5                              0                                0
#> 6                              0                                0
#> 7                              0                                0
#> 8                              0                                0
#> 9                              0                                0
#> 10                             0                                0
#>    commute_quietest_bicycle commute_quietest_bicycle_ebike
#> 1                         0                              0
#> 2                         0                              0
#> 3                         0                              0
#> 4                         0                              0
#> 5                         0                              0
#> 6                         0                              0
#> 7                         0                              0
#> 8                         0                              0
#> 9                         0                              0
#> 10                        0                              0
#>    commute_quietest_bicycle_go_dutch primary_fastest_bicycle
#> 1                                  0                       0
#> 2                                  0                       0
#> 3                                  0                       0
#> 4                                  0                       0
#> 5                                  0                       0
#> 6                                  0                       0
#> 7                                  0                       0
#> 8                                  0                       0
#> 9                                  0                       0
#> 10                                 0                       0
#>    primary_fastest_bicycle_ebike primary_fastest_bicycle_go_dutch
#> 1                              0                                0
#> 2                              0                                0
#> 3                              3                                0
#> 4                              1                                0
#> 5                              0                                0
#> 6                              0                                0
#> 7                              0                                0
#> 8                              0                                0
#> 9                              0                                0
#> 10                             0                                0
#>    primary_quietest_bicycle primary_quietest_bicycle_ebike
#> 1                         0                              0
#> 2                         0                              0
#> 3                         0                              0
#> 4                         0                              0
#> 5                         0                              0
#> 6                         0                              0
#> 7                         0                              0
#> 8                         0                              0
#> 9                         0                              0
#> 10                        0                              0
#>    primary_quietest_bicycle_go_dutch secondary_fastest_bicycle
#> 1                                  0                         0
#> 2                                  0                         0
#> 3                                  0                         0
#> 4                                  0                         0
#> 5                                  0                         0
#> 6                                  0                         0
#> 7                                  0                         0
#> 8                                  0                         0
#> 9                                  0                         0
#> 10                                 0                         0
#>    secondary_fastest_bicycle_ebike secondary_fastest_bicycle_go_dutch
#> 1                               14                                  0
#> 2                                0                                  0
#> 3                                0                                  0
#> 4                                0                                  0
#> 5                                0                                  0
#> 6                                0                                  0
#> 7                                0                                  0
#> 8                                0                                  0
#> 9                                0                                  0
#> 10                               0                                  0
#>    secondary_quietest_bicycle secondary_quietest_bicycle_ebike
#> 1                           0                                0
#> 2                           0                                0
#> 3                           0                                0
#> 4                           0                                0
#> 5                           0                                0
#> 6                           0                                0
#> 7                           0                                0
#> 8                           0                                0
#> 9                           0                                0
#> 10                          0                                0
#>    secondary_quietest_bicycle_go_dutch Gradient Quietness
#> 1                                    0        6        40
#> 2                                    0        7        40
#> 3                                    0        0        85
#> 4                                    0        0        50
#> 5                                    0        2        70
#> 6                                    0        6        70
#> 7                                    0        0        50
#> 8                                    0        0        70
#> 9                                    0        6        70
#> 10                                   0        0        80
#>                          geometry
#> 1  LINESTRING (-5.03374 54.902...
#> 2  LINESTRING (-5.03983 54.906...
#> 3  LINESTRING (-5.00087 54.904...
#> 4  LINESTRING (-5.00484 54.906...
#> 5  LINESTRING (-5.1116 54.8585...
#> 6  LINESTRING (-5.1116 54.8585...
#> 7  LINESTRING (-4.93467 54.861...
#> 8  LINESTRING (-4.90534 54.869...
#> 9  LINESTRING (-4.89788 54.868...
#> 10 LINESTRING (-4.90534 54.869...
```

``` r
tmap_mode("plot")
qtm(combined_network)
```

![](README_files/figure-gfm/unnamed-chunk-54-1.png)<!-- --> \### GWR

The following code runs the GWR model for all columns in the
`combined_network` object

``` r
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
```

``` r
val_results_gwr_last = cnet_val_gwr(cnet = combined_network,counts = sf_counts)
```

``` r
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
```

<img src="README_files/figure-gfm/val_last_table.png" width="672" />

## Applying AADT factors

The number of trips in the latest build has not been expanded to AADT
counts. Using the factors estimated in
(here)\[<https://github.com/nptscot/TT_-Scottish_Household_Survey/commit/2869f2ac845924bf51f3d5262c3f7579f362bdfa>\],
the values will be adjusted.

Also, the number of trips considered for each count will now include the
links within a 20 m buffer. It is worth noting that this comparison
assumes the intercept as 0.

### Loading the data

``` r
combined_network = st_read(dsn = "../npt/outputdata/combined_network_tile.geojson")
#> Reading layer `combined_network_tile' from data source 
#>   `C:\Users\ts18jpf\OneDrive - University of Leeds\03_PhD\00_Misc_projects\npt\outputdata\combined_network_tile.geojson' 
#>   using driver `GeoJSON'
#> Simple feature collection with 563466 features and 26 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -7.54343 ymin: 54.67982 xmax: -0.8615 ymax: 60.76383
#> Geodetic CRS:  WGS 84
```

Redefining the function to apply the expansion factor and the buffer for
the counts.

``` r
cnet_val_gwr_AADT = function(cnet, counts) {
  lst_names = names(cnet)[grepl("_", names(cnet))]
  
  # Creating buffer
  rnet_buffer30 = cnet |>
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
  
  
  # Buffer for the aggregated counts
  aggr_countrs_buffer = st_buffer(aggregated_counts, dist = 20)
  
  # Subsetting the network based on the buffer
  cnet_selected = cnet[aggr_countrs_buffer,]
  
  # Spatial join to detect the network links within the buffer
    Agg_network_siteID = cnet_selected |>
    st_join(aggr_countrs_buffer) |>
    st_drop_geometry() |>
    select(all_fastest_bicycle:secondary_quietest_bicycle_go_dutch,
           siteID) |>
    summarise(across(contains("bicycle"), sum), .by = siteID)
  
  
  # Join with the aggregated counts and factoring for AADT
  val_counts = aggregated_counts |>
    left_join(Agg_network_siteID, by = "siteID") |>
    # Hard-coded AADT factors
    mutate(across(contains("commute"), \(x) x * 1.7373),
           across(contains("ary"), \(x) x * 1.5947)) |>
    mutate(across(
      contains('all'),
      ~ get(str_replace(cur_column(), "^all", "commute")) +
        get(str_replace(cur_column(), "^all", "primary")) +
        get(str_replace(cur_column(), "^all", "secondary"))
    )) |>
    drop_na(all_fastest_bicycle:secondary_quietest_bicycle_go_dutch) |>
    as_Spatial()
  
  reg_rnet =
    lapply(lst_names, function(rnet_name) {
      
      form = as.formula(paste(rnet_name, "~ count_mean + 0"))
      
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
```

``` r
val_results_gwr_last_AADT = cnet_val_gwr_AADT(cnet = combined_network,counts = sf_counts)
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 321.4634 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 321.5054 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 321.4634 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 381.5359 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 381.5944 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 381.5359 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 384.0953 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 384.1464 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 384.0953 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 327.8732 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 328.0741 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 327.8732 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 399.2125 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 399.3746 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 399.2125 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 391.7051 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 391.8658 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 391.7051 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 321.1632 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 321.2059 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 321.1632 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 381.0076 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 381.0666 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 381.0076 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 382.6936 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 382.7468 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 382.6936 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 327.5555 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 327.7577 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 327.5555 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 398.5721 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 398.7332 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 398.5721 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 390.2169 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 390.3757 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 390.2169 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 104.6033 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 104.6524 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 213.5064 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 213.528 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 213.5064 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 237.8506 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 237.8561 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 112.5676 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 112.5958 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 235.536 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 235.5831 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 235.536 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 255.3882 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 255.4336 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 255.3882 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 150.975 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 150.9997 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 244.7891 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 244.8312 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 272.2366 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 272.2693 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 158.7162 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 158.7603 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 158.7162 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 256.5197 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 256.6043 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 256.5197 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 278.3833 
#> Adaptive bandwidth (number of nearest neighbours): 19 AICc value: 278.466 
#> Adaptive bandwidth (number of nearest neighbours): 18 AICc value: 278.3833
```

``` r
tbl_val_gwr_results_last_AADT = tibble(
  network = names(val_results_gwr_last_AADT),
  mean_coef = vapply(val_results_gwr_last_AADT, \(x) mean(x$SDF@data[, 1]), FUN.VALUE = 0),
  median_coef = vapply(val_results_gwr_last_AADT, \(x) median(x$SDF@data[, 1]), FUN.VALUE = 0)
  ,
  bw = vapply(val_results_gwr_last_AADT, \(x) x$GW.arguments$bw, FUN.VALUE = 0),
  R2 = vapply(val_results_gwr_last_AADT, \(x) x$GW.diagnostic$gw.R2, FUN.VALUE = 0)
)

tbl_val_gwr_results_last_AADT |>
  kbl(digits=4) |>
  kable_classic_2("hover", full_width = T) |>
  as_image(width = 7,file = "README_files/figure-gfm/val_last_table_AADT.png")
```

<img src="README_files/figure-gfm/val_last_table_AADT.png" width="672" />
The `all_fastest_bicycle` network seems to provide the best fit compared
to the observed counts, and considering that some trip purposes have not
been modelled. Also the intercept is relatively low. Nonetheless, a
closer look at the distribution of the coefficients `count_mean` (the
ratio of estimated trips over the observed counts) shows that in
Edinburgh the estimations are greater than the observed, while in
Glasgow the estimates are lower than observed.

``` r
best_network = val_results_gwr_last_AADT$all_fastest_bicycle$SDF |> st_as_sf()
tmap_mode("plot")
mymap = tm_shape(best_network)+
  tm_dots(col = "count_mean",palette = "RdYlGn",
          midpoint=1,size = 0.3,alpha=0.4)+
  tm_basemap()+
  tm_layout(legend.position = c("left", "top")) 

  tmap_save(mymap,filename = "README_files/figure-gfm/coef_maps.png",units = "cm",width = 6)
  
```

<img src="README_files/figure-gfm/coef_maps.png" width="708" />

    #>   |                                                          |                                                  |   0%  |                                                          |                                                  |   1%                      |                                                          |.                                                 |   2% [unnamed-chunk-131]                      |                                                          |..                                                |   3% [unnamed-chunk-132]  |                                                          |..                                                |   4%                      |                                                          |..                                                |   5% [unnamed-chunk-133]  |                                                          |...                                               |   5%                      |                                                          |...                                               |   6% [unnamed-chunk-134]  |                                                          |...                                               |   7%                      |                                                          |....                                              |   8% [unnamed-chunk-135]                      |                                                          |.....                                             |   9% [unnamed-chunk-136]  |                                                          |.....                                             |  10%                      |                                                          |.....                                             |  11% [unnamed-chunk-137]  |                                                          |......                                            |  12%                     [unnamed-chunk-138]  |                                                          |.......                                           |  13%                      |                                                          |.......                                           |  14% [unnamed-chunk-139]  |                                                          |.......                                           |  15%                      |                                                          |........                                          |  15% [unnamed-chunk-140]  |                                                          |........                                          |  16%                      |                                                          |........                                          |  17% [unnamed-chunk-141]  |                                                          |.........                                         |  18%                     [unnamed-chunk-142]  |                                                          |..........                                        |  19%                      |                                                          |..........                                        |  20% [unnamed-chunk-143]  |                                                          |..........                                        |  21%                      |                                                          |...........                                       |  22% [unnamed-chunk-144]                      |                                                          |............                                      |  23% [unnamed-chunk-145]  |                                                          |............                                      |  24%                      |                                                          |............                                      |  25% [unnamed-chunk-146]  |                                                          |.............                                     |  25%                      |                                                          |.............                                     |  26% [unnamed-chunk-147]  |                                                          |.............                                     |  27%                      |                                                          |..............                                    |  28% [unnamed-chunk-148]                      |                                                          |...............                                   |  29% [unnamed-chunk-149]  |                                                          |...............                                   |  30%                      |                                                          |...............                                   |  31% [unnamed-chunk-150]  |                                                          |................                                  |  32%                     [unnamed-chunk-151]  |                                                          |.................                                 |  33%                      |                                                          |.................                                 |  34% [unnamed-chunk-152]  |                                                          |.................                                 |  35%                      |                                                          |..................                                |  35% [unnamed-chunk-153]  |                                                          |..................                                |  36%                      |                                                          |..................                                |  37% [unnamed-chunk-154]  |                                                          |...................                               |  38%                     [unnamed-chunk-155]  |                                                          |....................                              |  39%                      |                                                          |....................                              |  40% [unnamed-chunk-156]  |                                                          |....................                              |  41%                      |                                                          |.....................                             |  42% [unnamed-chunk-157]                      |                                                          |......................                            |  43% [unnamed-chunk-158]  |                                                          |......................                            |  44%                      |                                                          |......................                            |  45% [unnamed-chunk-159]  |                                                          |.......................                           |  45%                      |                                                          |.......................                           |  46% [unnamed-chunk-160]  |                                                          |.......................                           |  47%                      |                                                          |........................                          |  48% [unnamed-chunk-161]                      |                                                          |.........................                         |  49% [unnamed-chunk-162]  |                                                          |.........................                         |  50%                      |                                                          |.........................                         |  51% [unnamed-chunk-163]  |                                                          |..........................                        |  52%                     [unnamed-chunk-164]  |                                                          |...........................                       |  53%                      |                                                          |...........................                       |  54% [unnamed-chunk-165]  |                                                          |...........................                       |  55%                      |                                                          |............................                      |  55% [unnamed-chunk-166]  |                                                          |............................                      |  56%                      |                                                          |............................                      |  57% [unnamed-chunk-167]  |                                                          |.............................                     |  58%                     [unnamed-chunk-168]  |                                                          |..............................                    |  59%                      |                                                          |..............................                    |  60% [unnamed-chunk-169]  |                                                          |..............................                    |  61%                      |                                                          |...............................                   |  62% [unnamed-chunk-170]                      |                                                          |................................                  |  63% [unnamed-chunk-171]  |                                                          |................................                  |  64%                      |                                                          |................................                  |  65% [unnamed-chunk-172]  |                                                          |.................................                 |  65%                      |                                                          |.................................                 |  66% [unnamed-chunk-173]  |                                                          |.................................                 |  67%                      |                                                          |..................................                |  68% [unnamed-chunk-174]                      |                                                          |...................................               |  69% [unnamed-chunk-175]  |                                                          |...................................               |  70%                      |                                                          |...................................               |  71% [unnamed-chunk-176]  |                                                          |....................................              |  72%                     [unnamed-chunk-177]  |                                                          |.....................................             |  73%                      |                                                          |.....................................             |  74% [unnamed-chunk-178]  |                                                          |.....................................             |  75%                      |                                                          |......................................            |  75% [unnamed-chunk-179]  |                                                          |......................................            |  76%                      |                                                          |......................................            |  77% [unnamed-chunk-180]  |                                                          |.......................................           |  78%                     [unnamed-chunk-181]  |                                                          |........................................          |  79%                      |                                                          |........................................          |  80% [unnamed-chunk-182]  |                                                          |........................................          |  81%                      |                                                          |.........................................         |  82% [unnamed-chunk-183]                      |                                                          |..........................................        |  83% [unnamed-chunk-184]  |                                                          |..........................................        |  84%                      |                                                          |..........................................        |  85% [unnamed-chunk-185]  |                                                          |...........................................       |  85%                      |                                                          |...........................................       |  86% [unnamed-chunk-186]  |                                                          |...........................................       |  87%                      |                                                          |............................................      |  88% [unnamed-chunk-187]                      |                                                          |.............................................     |  89% [unnamed-chunk-188]  |                                                          |.............................................     |  90%                      |                                                          |.............................................     |  91% [unnamed-chunk-189]  |                                                          |..............................................    |  92%                     [unnamed-chunk-190]  |                                                          |...............................................   |  93%                      |                                                          |...............................................   |  94% [unnamed-chunk-191]  |                                                          |...............................................   |  95%                      |                                                          |................................................  |  95% [unnamed-chunk-192]  |                                                          |................................................  |  96% [unnamed-chunk-193]  |                                                          |................................................  |  97%                      |                                                          |................................................. |  98% [unnamed-chunk-194] [unnamed-chunk-195]  |                                                          |..................................................|  99%                      |                                                          |..................................................| 100% [unnamed-chunk-196]
    #> [1] "counts.R"
