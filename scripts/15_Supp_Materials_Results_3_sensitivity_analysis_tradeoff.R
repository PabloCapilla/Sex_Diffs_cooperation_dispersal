###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2023-02-10
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### libraries and functions #####
##
##
pacman::p_load(dplyr, lubridate, ggplot2, 
               lme4, data.table, glmmTMB, DHARMa, 
               extrafont, gtsummary, gt)

source("../prep_repository/R_library/FUNCTION_drop1_output.R")
source("../prep_repository/R_library/FUNCTION_gather_dist_data.R")

#####


##
##
##### PROSPECTING DATA #####
##
##
foray <- readRDS(file = "./data/08_prospecting_rate_250m_15s.RDS")

##
## converting yday to date
foray$date <- ymd(as.Date(foray$day, origin = "2016-12-31"))
head(foray)

###
### retain forays from provisioning periods
foray_prov <- foray %>% 
  filter(prov == "YES") %>% 
  group_by(bird_id, date) %>% 
  summarise(day_foray = sum(forays),
            day_foray_duration = sum(total_duration)) %>% 
  mutate(bird_date = paste(bird_id, date, sep = "_"))

head(foray_prov)
length(unique(foray_prov$bird_id)) # data for 18 birds while provisioning

#####

##
##
##### PROVISIONING DATA #####
##
##
prov00 <- prov <- readRDS(file = "./data/01_provisioning_analysis_data.RDS")
prov00 <- prov00 %>% mutate(key = paste0(ind_id, group_id, date))
head(prov00)

##
##
##### PROV DATA FOR DAYS AND BIRDS WITH PROSPECTING RECORDED #####
##
##

## creating key columns

# bird and days with prospecting data
prospect_bird_day <- unique(foray_prov$bird_date)

# creating key column in prov
prov <- prov %>% 
  mutate(bird_date = paste(ind_id, date, sep = "_"))

# filtering prov data
prov_prov <- prov[prov$bird_date %in% prospect_bird_day,]

# merging data sets
length(unique(prov_prov$ind_id)) 

# birds in prov
unique(prov_prov$ind_id)

# birds in foray
unique(foray_prov$bird_id)

# missing in one?
unique(foray_prov$bird_id)[!unique(foray_prov$bird_id) %in% unique(prov_prov$ind_id)]


prov00[prov00$ind_id == "4H-53068",] # no prov data for this bird

prov00[prov00$ind_id == "4H-55327",]
foray_prov[foray_prov$bird_id == "4H-55327",] # no prov data for this clutch in the prov

prov00[prov00$ind_id == "4H-55361",]
foray_prov[foray_prov$bird_id == "4H-55361",] # no prov data for this clutch in the prov

# combining data sets
df_analysis <- prov_prov
df_analysis <- left_join(x = prov_prov, 
                         y = foray_prov %>% select(-date), 
                         by = "bird_date") %>% 
  mutate(foray_made = ifelse(day_foray == 0, 0, 1))

length(unique(df_analysis$bird_id))
head(df_analysis)
table(df_analysis$foray_made)
table(df_analysis$ind_id)

#####

##
##
##### Calculate prospecting in the 7 days before provisioning periods #####
##
## 
vec_days <- seq(1, 7, 1)
tbl_df <- data.frame(vec_days = vec_days,
                     prosp_coef = NA,
                     se_prosp_coef = NA)

for(d in 1:nrow(tbl_df)){
  
  # list of prospecting birds
  df_analysis$forays_before <- as.numeric(NA)
  df_analysis$foray_rate_before <- as.numeric(NA)
  df_analysis$foray_duration_before <- as.numeric(NA)

    for (i in 1:nrow(df_analysis)){
    ind_loop <- df_analysis$ind_id[i] 
    date_loop <- df_analysis$date[i]
    clutch_loop <- df_analysis$clutch_id[i]
    
    # prospecting 5 days before prospecting
    df_loop <-  foray %>% 
      filter(bird_id == ind_loop) %>% 
      filter(date <= date_loop & date > date_loop - tbl_df$vec_days[d])
    
    days_obs <- length(unique(df_loop$date))
    total_forays <- sum(df_loop$forays)
    duration_forays <- sum(df_loop$total_duration)
    
    df_analysis$forays_before[i] <- total_forays
    df_analysis$foray_rate_before[i] <- total_forays / days_obs
    df_analysis$foray_duration_before[i] <- duration_forays
  }
  
  # number of morning forays continuous
  m2 <- glmer(provisioning_feeds ~ 
                foray_rate_before +
                nestlings +
                ba_day + 
                sex +
                (1|group_id) + 
                (1|ind_id),
              data = df_analysis,
              offset = log(df_analysis$session_duration),
              family = 'poisson',
              na.action = "na.fail")
  
  tbl_df$prosp_coef[d] <- {summary(m2)}$coefficients[2,1]
  tbl_df$se_prosp_coef[d] <- {summary(m2)}$coefficients[2,2]

}

#####

##
##
##### Plot #####
##
##
tradeoff_sensi <- ggplot(data = tbl_df,
                         aes(x = vec_days, 
                             y = prosp_coef)) +
  geom_errorbar(aes(ymin = prosp_coef - se_prosp_coef,
                    ymax = prosp_coef + se_prosp_coef),
                width = 0) +
  geom_point(size = 5,
             color = "black",
             fill = "#2b8cbe",
             shape = 21) +
  geom_hline(yintercept = 0,
             linetype = 3) +
  theme_bw() +
  labs(y = expression(atop("Effect of individual prospecting rate (foray / day)", 
                           "on provisioning rate (feeds / hour)")), 
       x = expression(atop("Number of days for calculation", 
                           "of prospecting rates (forays / day)"))) + 
  theme(axis.title.x = element_text(family = "Arial", color = "black", size = 10),
        axis.title.y = element_text(family = "Arial", color = "black", size = 10),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = vec_days, labels = vec_days)

#####

ggsave(filename = "./plots/Figure_S3.png", 
       plot = tradeoff_sensi, 
       units = "mm",
       device = "jpeg", 
       width = 125,
       height = 95)


