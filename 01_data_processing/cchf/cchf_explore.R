## EXPLORING THE CCHF DATA

# Load packages 
library(dplyr)
library(lubridate)
library(sf)
library(readxl)
library(ggplot2)

## View spatially in maps 

# Load CCHF data 
cchf_clean <- read_excel("cchf_clean.xlsx")

# Import place names key to match to GADM 
districts <- read_delim("02_data_harmonisation/irq_places_key/district-Table 1.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)
provinces <- read_delim("02_data_harmonisation/irq_places_key/province-Table 1.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)

districts <- districts %>% rename("District" = Distinct) # spelling error in my metadata

# Merge to my province names 
cchf <- left_join(cchf_clean, provinces, by = "Province")

# Subset the year with no province names 
df_na <- cchf %>% filter(is.na(Province))
cchf <- cchf %>% filter(!is.na(Province))

# Merge to the district names 
df_na <- df_na %>% select(-GID_1, -NAME_1)
df_na <- left_join(df_na, districts, by = "District")

# Rebind back to original data 
df_na <- df_na %>% select(-GID_2, -NAME_2)
cchf <- rbind(cchf, df_na)

# Merge to map and plot 
irq <- read_sf("gadm41_IRQ_1.shp", stringsAsFactors = F)
map <- cchf %>% group_by(GID_1, NAME_1) %>% tally()
map <- left_join(irq, map)
ggplot(map) +
  geom_sf(aes(fill = n)) +
  geom_sf_text(aes(label = NAME_1), size = 5, fontface = "bold", color = "#ff6f61") +   # add labels
  theme_minimal() +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1,
    na.value = "#f5f0e6"   # light beige
  ) +
  labs(fill = "Total Cases") +
  theme(
    text = element_text(face = "bold"),
    legend.position = "top"
  )

# Merge to map and plot for ADM2 
map2 <- left_join(cchf, districts, by = "District")
map2 <- map2 %>% group_by(GID_2, NAME_2) %>% tally()
irq <- read_sf("gadm41_IRQ_2.shp", stringsAsFactors = F)
map2 <- left_join(irq, map2)
ggplot(map2) + geom_sf(aes(fill = n)) + theme_minimal() + 
  scale_fill_distiller(palette = "Reds", direction = 1) + 
  labs(fill = "Total Cases") + 
  theme(text = element_text(face = "bold"))

## View temporally in a time-series  

# Separate dates into columns
cchf <- cchf %>%
  mutate(
    year    = year(Dateofonset),
    month   = month(Dateofonset),
    day     = day(Dateofonset),
    epiweek = epiweek(Dateofonset)
    )

# Create time-series data at different temporal scales 
time_series <- cchf %>% group_by(year, month) %>% mutate(monthly_count = n()) %>% ungroup() %>% 
  group_by(year, month, day) %>% mutate(daily_count = n()) %>% ungroup() %>% 
  group_by(year, epiweek) %>% mutate(weekly_count = n()) %>% ungroup()

# Plot series 
ggplot(time_series, aes(x = Dateofonset, y = monthly_count)) + geom_line(color = "dodgerblue") + 
  theme_minimal() + labs(y = "Monthly Cases") + theme(text = element_text(face = "bold"))
ggplot(time_series, aes(x = Dateofonset, y = daily_count)) + geom_line(color = "dodgerblue") + 
  theme_minimal() + labs(y = "Daily Cases") + theme(text = element_text(face = "bold"))

ggplot(time_series, aes(x = Dateofonset, y = weekly_count)) + 
  geom_bar(stat = "identity", color = "dodgerblue") + 
  theme_minimal() + labs(y = "Weekly Cases", x = "Date of Onset") + 
  theme(text = element_text(face = "bold"))

# Accounting for gender 
time_series <- cchf %>% group_by(year, month, Gender) %>% mutate(monthly_count = n()) %>% ungroup()
  
# Plot series 
ggplot(time_series, aes(x = Dateofonset, y = monthly_count, color = Gender)) + geom_line() + 
  theme_minimal() + labs(y = "Monthly Cases") + theme(text = element_text(face = "bold"))

# Accounting for age
cchf <- cchf %>% mutate(age_group = case_when(Ageyear < 16 ~ "Children",
                                              Ageyear > 15 & Ageyear < 60 ~ "Adults",
                                              Ageyear > 59 ~ "Older Adults")) 

time_series <- cchf %>% group_by(year, month, age_group) %>% mutate(monthly_count = n()) %>% ungroup()

# Plot series 
time_series <- na.omit(time_series)
ggplot(time_series, aes(x = Dateofonset, y = monthly_count, color = age_group)) + geom_line() + 
  theme_minimal() + labs(y = "Monthly Cases") + theme(text = element_text(face = "bold")) + 
  labs(color = "Age Group")

## View age/gender in bars 

bar1 <- cchf %>% group_by(Gender) %>% tally()
bar2 <- cchf %>% group_by(age_group) %>% tally()
bar2 <- na.omit(bar2)

p1 <- ggplot(bar1, aes(x = Gender, y = n)) + geom_bar(stat = "identity", fill = "#FFB6C1", color = "#FFB6C1") + 
  labs(y = "Total Cases") + theme(text = element_text(face = "bold")) + theme_minimal() + 
  theme(text = element_text(face = "bold"))
bar2$age_group <- factor(bar2$age_group, levels=c("Children", "Adults", "Older Adults"))
p2 <- ggplot(bar2, aes(x = age_group , y = n)) + geom_bar(stat = "identity", fill = "#ffae42", color = "#ffae42") + 
  labs(y = "Total Cases", x = "Age Groups") + theme(text = element_text(face = "bold")) + theme_minimal() + 
  theme(text = element_text(face = "bold"))

ggpubr::ggarrange(p1, p2)


