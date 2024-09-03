pacman::p_load(leaflet, mapview, viridisLite, 
               leaflet.extras2, tidyverse, sf, patchwork)

rm(list = ls())

setwd("G:/My Drive/GIS")

# LSOAs --------------------------------------------

lsoas <- read_sf("LSOA_2011_EW_BSC_V4.shp")

# leeds imd score ----------------------------------

leeds_imd <- read_csv("imd_scores.csv") %>% 
  filter(lad_name == "Leeds")

leeds_imd <- lsoas %>% 
  right_join(leeds_imd, by = c("LSOA11CD" = "lsoa_code"))

# shapefile of postcode districts -----------------------------------

referrals <- read_csv("G:/My Drive/Ageing Better HIS evaluation/Leeds/referrals_postcode.csv")

districts <- read_sf("G:/My Drive/Ageing Better HIS evaluation/Leeds/GB_Postcodes/PostalDistrict.shp")

leeds_districts <- districts %>% 
  right_join(referrals,
             by = c("PostDist" = "postcode_district"))

leeds_districts <- st_transform(leeds_districts, 
                                crs = st_crs(leeds_imd))

# intersecting polygons ----------------------------------------

inters <- st_intersection(leeds_imd, leeds_districts)

inters <- inters %>% 
  mutate(id = str_c(LSOA11CD, ".", PostDist),
         .before = 1)

inters %>% 
  group_by(PostDist) %>% 
  summarise(mean_imd = mean(imd_score, na.rm = T),
            .groups = "drop") %>% 
  ggplot(aes(fill = mean_imd)) +
  geom_sf(alpha = .8) +
  scale_fill_viridis_c() +
  labs(fill = "IMD")

inters %>%
  group_by(PostDist) %>% 
  summarise(referrals = mean(Total, na.rm = T)) %>% 
  ggplot(aes(fill = referrals)) +
  geom_sf(alpha = .8) +
  scale_fill_viridis_c() +
  labs(fill = "Referrals")

# interactive map -------------------------------------------

Leeds <- inters %>% 
  group_by(PostDist) %>% 
  summarise(IMD = mean(imd_score, na.rm = T),
            Referrals = mean(Total, na.rm = T),
            .groups = "drop")

Leeds <- st_transform(Leeds, 4326)

pal <- colorRampPalette(viridis(n=9))
at1 <- seq(min(Leeds$IMD),max(Leeds$IMD),length.out = 8)
at2 <- seq(min(Leeds$Referrals),max(Leeds$Referrals),length.out = 8)

m1 <- mapview(Leeds, zcol = "IMD", map.types = "CartoDB.Positron",
              col.regions = pal, at = at1)
m2 <- mapview(Leeds, zcol = "Referrals", map.types = "CartoDB.Positron",
              col.regions = pal, at = at2)

m1 | m2
