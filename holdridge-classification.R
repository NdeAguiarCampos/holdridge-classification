#-------------------------------------------------------------------#
# THESIS CHAPTER: 1. Systematic review
# SCRIPT 1: Assigning Holdridge Life Zones to the study sites
# AUTHOR: Natalia de Aguiar Campos
# DATE: 15 April 2024
#-------------------------------------------------------------------#

library(readxl)
library(tidyverse)
library(stars)
library(sp)
library(raster)
library(Ternary)
library(climenv)
library(writexl)

#---------------------------#
# PART 1/2 OF THE SCRIPT ----
#---------------------------#

# Step 1 - Load functions ----
calculate_biotemperature <- function(df) {
  df[df < 0] <- 0
  df[df > 30] <- 0
  
  biotemperature <- rowMeans(df)
  
  return(biotemperature)
}

assign_preliminary_altitudinal_belts <- function(biotemperature) {
  alt_boundaries <- list(
    Alvar = c(0.000, 1.500),
    Alpine = c(1.501, 3.000),
    Subalpine = c(3.001, 6.000),
    Montane = c(6.001, 12.000),
    `Lower montane` = c(12.001, 18.000),
    Premontane = c(18.001, 24.000),
    Lowland = c(24.001, 30.000)
  )
  
  for (alt_belt in names(alt_boundaries)) {
    temp_range <- alt_boundaries[[alt_belt]]
    
    if (biotemperature >= temp_range[1] && biotemperature <= temp_range[2]) {
      return(alt_belt)
    }
  }
  
  return(NA)
}

assign_life_zones_humidity_provinces <- function(preliminary_altitudinal_belt, precipitation) {
  case_when(
    preliminary_altitudinal_belt == "Lowland" & precipitation <= 125 ~ "Superarid desert",
    preliminary_altitudinal_belt == "Lowland" & precipitation <= 250 ~ "Perarid desert scrub",
    preliminary_altitudinal_belt == "Lowland" & precipitation <= 500 ~ "Arid thorn woodland",
    preliminary_altitudinal_belt == "Lowland" & precipitation <= 1000 ~ "Semiarid very dry forest",
    preliminary_altitudinal_belt == "Lowland" & precipitation <= 2000 ~ "Subhumid dry forest",
    preliminary_altitudinal_belt == "Lowland" & precipitation <= 4000 ~ "Humid moist forest",
    preliminary_altitudinal_belt == "Lowland" & precipitation <= 8000 ~ "Perhumid wet forest",
    preliminary_altitudinal_belt == "Lowland" & precipitation <= 16000 ~ "Superhumid rain forest",
    (preliminary_altitudinal_belt %in% c("Premontane", "Lower montane")) & precipitation <= 125 ~ "Perarid desert",
    (preliminary_altitudinal_belt %in% c("Premontane", "Lower montane")) & precipitation <= 250 ~ "Arid desert scrub",
    (preliminary_altitudinal_belt %in% c("Premontane", "Lower montane")) & precipitation <= 500 ~ "Semiarid thorn steppe/woodland",
    (preliminary_altitudinal_belt %in% c("Premontane", "Lower montane")) & precipitation <= 1000 ~ "Subhumid dry forest",
    (preliminary_altitudinal_belt %in% c("Premontane", "Lower montane")) & precipitation <= 2000 ~ "Humid moist forest",
    (preliminary_altitudinal_belt %in% c("Premontane", "Lower montane")) & precipitation <= 4000 ~ "Perhumid wet forest",
    (preliminary_altitudinal_belt %in% c("Premontane", "Lower montane")) & precipitation <= 8000 ~ "Superhumid rain forest",
    preliminary_altitudinal_belt == "Montane" & precipitation <= 125 ~ "Arid desert",
    preliminary_altitudinal_belt == "Montane" & precipitation <= 250 ~ "Semiarid desert scrub",
    preliminary_altitudinal_belt == "Montane" & precipitation <= 500 ~ "Subhumid steppe",
    preliminary_altitudinal_belt == "Montane" & precipitation <= 1000 ~ "Humid moist forest",
    preliminary_altitudinal_belt == "Montane" & precipitation <= 2000 ~ "Perhumid wet forest",
    preliminary_altitudinal_belt == "Montane" & precipitation <= 4000 ~ "Superhumid rain forest",
    preliminary_altitudinal_belt == "Subalpine" & precipitation <= 125 ~ "Semiarid desert",
    preliminary_altitudinal_belt == "Subalpine" & precipitation <= 250 ~ "Subhumid dry scrub",
    preliminary_altitudinal_belt == "Subalpine" & precipitation <= 500 ~ "Humid moist forest",
    preliminary_altitudinal_belt == "Subalpine" & precipitation <= 1000 ~ "Perhumid wet forest",
    preliminary_altitudinal_belt == "Subalpine" & precipitation <= 2000 ~ "Superhumid rain forest",
    preliminary_altitudinal_belt == "Alpine" & precipitation <= 125 ~ "Subhumid dry tundra",
    preliminary_altitudinal_belt == "Alpine" & precipitation <= 250 ~ "Humid moist tundra",
    preliminary_altitudinal_belt == "Alpine" & precipitation <= 500 ~ "Perhumid wet tundra",
    preliminary_altitudinal_belt == "Alpine" & precipitation <= 1000 ~ "Superhumid rain tundra",
    preliminary_altitudinal_belt == "Alvar" & precipitation <= 500 ~ "Polar desert",
    TRUE ~ NA_character_
  )
}

assign_latitudinal_region <- function(sealevel_biotemperature) {
  lat_boundaries <- list(
    Polar = c(0.000, 1.500),
    Subpolar = c(1.501, 3.000),
    Boreal = c(3.001, 6.000),
    `Cool temperate` = c(6.001, 12.000),
    `Warm temperate` = c(12.001, 18.000),
    Subtropical = c(18.001, 24.000),
    Tropical = c(24.001, 30.000)
  )
  
  for (lat_belt in names(lat_boundaries)) {
    temp_range <- lat_boundaries[[lat_belt]]
    
    if (sealevel_biotemperature >= temp_range[1] && sealevel_biotemperature <= temp_range[2]) {
      return(lat_belt)
    }
  }
  
  return(NA)
}


# Step 2 - Import the data frame with site data ----

site_data <- read_xlsx(path = "./Input/PREPARING_SUBMISSION_Sankey_Data.xlsx",
                       sheet = "Geography") %>%
  filter(Use_site == "Yes") %>% 
  dplyr::select(
    Site_key, 
    Long,
    Lat,
    Elevation, 
    `Mean annual rainfall`,
    `Mean annual temperature`,
    `Vegetation type (revised)`
  ) %>%
  drop_na(Lat, Long) %>%
  rename("prec_review" = `Mean annual rainfall`,
         "temp_review" = `Mean annual temperature`,
         "vegetation_review" = `Vegetation type (revised)`) %>%
  data.frame()

# Step 3 - Transform site coordinates into simple features ----

site_coord <- st_as_sf(site_data, coords = c(2, 3))

# Step 4 - Extract environmental data at 30-arcsecond resolution from the site coordinates ----

env_list <- ce_extract(
  path = "./Input/Environmental_data/worldclim_monthly_1970-2000/30_arcseconds/",
  location = site_coord,
  location_g = "Site_key",
  c_source = "WorldClim",
  var = "all"
)

# Step 5 - Classify each site coordinate into its Life Zone classification based on the env. data ----

site_classification <- as.data.frame(calculate_biotemperature(env_list$tavg_m)) %>% 
  rename_with(~ "biotemperature", 1) %>% 
  rownames_to_column(., var = "Site_key") %>% 
  full_join(., site_data, by = "Site_key") %>% 
  relocate(biotemperature, .after = "vegetation_review") %>%
  mutate(precipitation = rowSums(env_list$prec_m),
         biotemperature = round(biotemperature, 3),
         elevation = env_list[["elev"]][,1]) %>% 
  mutate(sealevel_biotemperature = round(biotemperature + (0.006 * elevation), 3)) %>% 
  mutate(sealevel_biotemperature = if_else(sealevel_biotemperature <= 30.000,
                                           sealevel_biotemperature, 30)) %>% 
  mutate(preliminary_altitudinal_belt = 
           apply(.[, "biotemperature", drop = FALSE],
                 1, function(x) assign_preliminary_altitudinal_belts(x))) %>% 
  mutate(life_zone = assign_life_zones_humidity_provinces(preliminary_altitudinal_belt, precipitation)) %>% 
  mutate(humidity_province = tolower(word(life_zone, 1)),
         life_zone = str_remove(life_zone, '(\\w+\\s+){1}')) %>% 
  relocate(humidity_province, .before = "life_zone") %>% 
  mutate(latitudinal_region = apply(.[, "sealevel_biotemperature", drop = FALSE],
                                    1, function(x) assign_latitudinal_region(x))) %>% 
  relocate(latitudinal_region, .before = "preliminary_altitudinal_belt") %>% 
  mutate(sealevel_altitudinal_belt = apply(.[, "sealevel_biotemperature", drop = FALSE],
                                           1, function(x) assign_preliminary_altitudinal_belts(x))) %>% 
  mutate(altitudinal_belt = tolower(if_else(preliminary_altitudinal_belt == sealevel_altitudinal_belt,
                                            "Lowland", preliminary_altitudinal_belt))) %>% 
  dplyr::select(-preliminary_altitudinal_belt, -sealevel_altitudinal_belt) %>% 
  relocate(altitudinal_belt, .after = "latitudinal_region")

# Step 6 - Save the results in a CSV file ----

write_xlsx(site_classification, "./Results/Excel tables/Holdridge classification by study site_Updated.xlsx")


#---------------------------#
# PART 2/2 OF THE SCRIPT ----
#---------------------------#

library(readxl)
library(tidyverse)
library(ggsankey)

site_data <- read_xlsx("./Results/Excel tables/Holdridge classification by study site_Updated.xlsx") %>% 
  data.frame()

site_class <- site_data %>% 
  dplyr::select(
    latitudinal_region, 
    altitudinal_belt, 
    # humidity_province, 
    life_zone
  ) %>% 
  rename(
    `Latitudinal region` = "latitudinal_region",
    `Altitudinal belt` = "altitudinal_belt",
    # `Humidity province` = "humidity_province",
    `Life zone` = "life_zone"
  )

site_class_long <- site_class %>% 
  make_long(
    `Latitudinal region`,
    `Altitudinal belt`,
    # `Humidity province`,
    `Life zone`
  ) %>% 
  mutate(
    x_node_key = paste0(x, "_", node)
  )

node_numbers <- site_class_long %>%
  group_by(x_node_key) %>%
  summarise(n = n()) %>%
  mutate(percentage = 100*n/70) %>% # 70 is the total number of study sites
  ungroup()

site_sankey <- merge(site_class_long, node_numbers, by.x = 'x_node_key', by.y = 'x_node_key', all.x = T)

my_labels <- c("Latitudinal region",
               "Altitudinal belt",
               # "Humidity\nprovince",
               "Life zone")

site_sankey$node <- factor(site_sankey$node, levels = c(
  "Tropical",
  "Subtropical",
  "lowland",
  "premontane",
  "lower montane",
  "rain forest",
  "wet forest",
  "moist forest",
  "dry forest",
  "very dry forest",
  "thorn woodland"
))

sankey_gg <- 
  ggplot(site_sankey, aes(x = x, 
                         next_x = next_x, 
                         node = node, 
                         next_node = next_node, 
                         fill = factor(node), 
                         label = paste0(node, " (", round(percentage, 2),"%)"),
                         label = node
  )) +
  geom_sankey(flow.alpha = .7,
              node.color = "black",
              color = "black",) +
  geom_sankey_label(
    size = 3.5, 
    color = "black", 
    fill = "white", 
  ) +
  scale_fill_viridis_d(option = "turbo")+
  theme_sankey(base_size = 12) +
  labs(x = NULL) +
  scale_x_discrete(labels = my_labels,
                   position = "top") + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black"))

sankey_gg

ggsave(sankey_gg,
      filename = "./Results/Holdridge_Site classification_With label numbers.jpeg",
      width = 18,
      height = 12,
      units = c("cm"),
      dpi = 600)

#-----------------------------------------------------------------------------------#

# PLOTTING THE HOLDRIDGE LIFE ZONE CHART

# STEP 1 - Import the compiled site data ----

site_data <- read_xlsx(path = "./Input/PREPARING_SUBMISSION_Sankey_Data.xlsx",
                       sheet = "Geography") %>%
  filter(Use_site == "Yes") %>% 
  dplyr::select(
    Site_key, 
    Long,
    Lat,
    Elevation, 
    `Mean annual rainfall`,
    `Mean annual temperature`,
    `Vegetation type (revised)`
  ) %>%
  drop_na(Lat, Long) %>%
  rename("prec_review" = `Mean annual rainfall`,
         "temp_review" = `Mean annual temperature`,
         "vegetation_review" = `Vegetation type (revised)`) %>%
  data.frame()

site_coord <- SpatialPointsDataFrame(site_data[,2:3], site_data)
crs(site_coord) <- "+init=epsg:4326"
plot(site_coord)


# STEP 2 - Import data from WorldClim and Zomer et al. (2022) ----

temp <- raster("./Input/Environmental_data/worldclim_monthly_1970-2000/30_arcseconds/wc2.1_30s_bio/wc2.1_30s_bio_1.tif") # WC
prec <- raster("./Input/Environmental_data/worldclim_monthly_1970-2000/30_arcseconds/wc2.1_30s_bio/wc2.1_30s_bio_12.tif") # WC
pet <- raster("./Input/Environmental_data/ET0_Aridity/Global-AI_ET0_annual_v3/Global-AI_ET0_v3_annual/et0_v3_yr.tif") # Zomer

temp_wc <- as.data.frame(raster::extract(temp, site_coord))
names(temp_wc) <- "temp_wc"

prec_wc <- as.data.frame(raster::extract(prec, site_coord))
names(prec_wc) <- "prec_wc"

pet_zomer <- as.data.frame(raster::extract(pet, site_coord))
names(pet_zomer) <- "pet_zomer"


# STEP 3 - Fill the NAs in precipitation and temperature with WorldClim, attach PET data ----

site_coord_env <- as.data.frame(cbind(site_coord, temp_wc, prec_wc, pet_zomer)) %>%
  dplyr::select(-Long.1, -Lat.1) %>% 
  mutate(pet_ratio = pet_zomer/prec_wc) %>% # The PET ratio is ONLY used in the function HoldridgePlot()
  drop_na(prec_wc, temp_wc, pet_ratio)

# STEP 4 - Plot the Holdridge diagram without the spectrum legend ----

setwd('..')

# UPPER
# jpeg(file = "./Figures/Holdrige by site_70 sites_Upper part_LARGER2.jpeg",
#      width = 30, height = 20, units = "cm", res = 300)

par(mar = c(0, 0, 0, 0))

{
  lifezoneLabels_Upper <- c(
    "NA",
    "NA",
    "NA",
    "desert",
    "desert",
    "desert",
    "NA",
    "dry\ntundra",
    "moist\ntundra",
    "wet\ntundra",
    "rain\ntundra",
    "NA",
    "desert",
    "dry\nscrub",
    "moist\nforest",
    "wet\nforest",
    "rain\nforest",
    "NA",
    "desert",
    "desert\nscrub",
    "steppe",
    "moist\nforest",
    "wet\nforest",
    "rain\nforest",
    "NA",
    "desert",
    "desert\nscrub",
    "thorn\nsteppe/\nwoodland",
    "dry\nforest",
    "moist\nforest",
    "wet\nforest",
    "rain\nforest",
    "NA",
    "desert",
    "desert\nscrub",
    "thorn\nwoodland",
    "very\ndry\nforest",
    "dry\nforest",
    "moist\nforest",
    "wet\nforest",
    "rain\nforest"
  ) 
}

HoldridgePlot(hex.labels = lifezoneLabels_Upper,
              hex.cex = 0.8,
              lab.cex = 2,
              axis.cex = 1.2) # Uses PET from the literature and MAP from wc
HoldridgeBelts()

HoldridgePoints(site_coord_env$pet_ratio, site_coord_env$prec_wc,
                col = hcl.colors(25)[abs(site_coord_env$Lat) + 1],
                # lwd = 8,
                cex = 1.5,
                pch = 19)

# dev.off()

# LOWER
# jpeg(file = "./Figures/Holdrige by site_70 sites_Lower part_LARGER2.jpeg",
#      width = 30, height = 20, units = "cm", res = 300)

par(mar = c(0, 0, 0, 0))

{
  lifezoneLabels_Lower <- c(
    "desert",
    "desert",
    "desert",
    "dry\ntundra",
    "moist\ntundra",
    "wet\ntundra",
    "rain\ntundra",
    "desert",
    "dry\nscrub",
    "moist\nforest",
    "wet\nforest",
    "rain\nforest",
    "desert",
    "desert\nscrub",
    "steppe",
    "moist\nforest",
    "wet\nforest",
    "rain\nforest",
    "desert",
    "desert\nscrub",
    "thorn\nsteppe/\nwoodland",
    "dry\nforest",
    "moist\nforest",
    "wet\nforest",
    "rain\nforest",
    "desert",
    "desert\nscrub",
    "thorn\nwoodland",
    "very dry\nforest",
    "dry\nforest",
    "moist\nforest",
    "wet\nforest",
    "rain\nforest"
  ) 
}

HoldridgePlot(hex.labels = lifezoneLabels_Lower,
              hex.cex = 0.8,
              lab.cex = 2,
              axis.cex = 1.2) # Uses PET from the literature and MAP from wc
HoldridgeBelts()

HoldridgePoints(site_coord_env$pet_ratio, site_coord_env$prec_wc,
                col = hcl.colors(25)[abs(site_coord_env$Lat) + 1],
                # lwd = 8,
                cex = 1.5,
                pch = 19)

# dev.off()

# STEP 5 - Plot the spectrum legend by itself

# jpeg(file = "./Figures/Spectrum legend.jpeg",
#      width = 18, height = 12, units = "cm", res = 300)

plot.new()

PlotTools::SpectrumLegend(
  "topright", bty = "n", # No box
  horiz = TRUE, # Horizontal
  x.intersp = -0.5, # Squeeze in X direction
  legend = paste0(seq(0, 25, 5), "°"),
  palette = hcl.colors(25),
  title = "Latitude"
)

# dev.off()

# Later, I assembled the upper and lower parts, along with the spectrum legend in PowerPoint.










