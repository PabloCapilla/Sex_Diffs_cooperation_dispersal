###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2023-09-05
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
pacman::p_load(dplyr, lubridate, ggplot2, cowplot,
               lme4, data.table, glmmTMB, DHARMa, 
               extrafont, gtsummary, gt)

source("../prep_repository/R_library/FUNCTION_drop1_output.R")
source("../prep_repository/R_library/FUNCTION_gather_dist_data.R")

#####

##
##
##### data #####
##
##

##
## bird id list
b_id <- read.table("../prep_repository/raw_data/bird_ID_list.txt", header = TRUE)
b_id %>% 
  arrange(bird_id)

b_id <- b_id %>%
  filter(bird_id != "4H-48951") %>% # dom not present in table
  filter(bird_id != "4H-55183") %>% # dom not present in table
  filter(bird_id != "4H-76107") %>% # fledgy not present in table
  filter(bird_id != "4H-55343") %>% # not present in table
  filter(bird_id != "4H-68764") %>% # not present in table
  filter(bird_id != "4H-68762") %>% # not present in table
  filter(bird_id != "4H-55338") %>% # not present in table
  filter(bird_id != "4H-68644") # not present in table

#####

##
##
##### Set up loop to do sensitivity analysis #####
##
##

# distance
distance_value <- seq(from = 100,
                      to = 300, 
                      by = 25)

# collecting data
## forays for distance threshold 250 m and time threshold 15seg
path_files <- "../../5_Prospecting/R_analysis/data/4_flights/counts/"

# data frame to store results
tbl_df <- data.frame(vec_days = distance_value,
                     n_forays = NA,
                     sex_coef = NA,
                     se_sex_coef = NA,
                     p_sex = NA,
                     prov_coef = NA,
                     se_prov_coef = NA,
                     p_prov = NA)

# loop to run through different datasets/models
for (d in 1:length(distance_value)){
  
  # read data
  data <- gather_dist_data(distance_threshold = distance_value[d],
                           duration_threshold = 15,
                           path_to_files = path_files)
  ##
  ## amend typo in data
  data[data$bird_id == "4h-68783","bird_id"] <- "4H-68783"
  
  # filter to retain only good quality data #####
  data <- data[data$bird_id %in% b_id$bird_id,]
  
  ## there is a female in HY with missing prov info
  ## I checked it manually and include it now
  data$prov <- ifelse(is.na(data$prov), "NO", as.character(data$prov))
  
  ##
  ## include bird age
  data <- left_join(x = data, 
                    y = b_id %>% select(bird_id, age2), 
                    by = "bird_id") %>% 
    rename(age_days = age2)
  
  # n of forays
  n_forays <- sum(data$forays)
  
  ##
  ## summarise data per day
  data_days <- data %>% 
    group_by(bird_id, day, sex, prov, home) %>%
    summarise(total_foray_day = sum(forays), n_days = n(),
              age_days = mean(age_days))
  data_days$sex_name <- ifelse(data_days$sex == "FEMALE", "female", "male")
  data_days$obs_id <- 1:nrow(data_days)
  data_days$age_years <- data_days$age_days/365
  
  # model
  model_loop <- glmer(total_foray_day ~
                        sex_name +
                        prov +
                        age_years +
                        (1|obs_id) + 
                        (1|bird_id) + 
                        (1|home),
                      data = data_days,
                      family = "poisson",
                      na.action = "na.fail")

    tbl_df$n_forays[d] <- n_forays
  tbl_df$sex_coef[d] <- {summary(model_loop)}$coefficients[2,1]
  tbl_df$se_sex_coef[d] <- {summary(model_loop)}$coefficients[2,2]
  tbl_df$p_sex[d] <- {summary(model_loop)}$coefficients[2,4]
  tbl_df$prov_coef[d] <- {summary(model_loop)}$coefficients[3,1]
  tbl_df$se_prov_coef[d] <- {summary(model_loop)}$coefficients[3,2]
  tbl_df$p_prov[d] <- {summary(model_loop)}$coefficients[3,4]
}

##
##
##### Plot #####
##
##
tbl_df <- tbl_df %>% 
  mutate(sign_sex = ifelse(p_sex < 0.05, 'YES', 'NO'),
         sign_prov = ifelse(p_prov < 0.05, 'YES', 'NO'))
  
  
  # plot for n forays
n_forays <- ggplot(data = tbl_df,
                   aes(x = vec_days, 
                       y = n_forays/1000)) +
  geom_col(color = "black",
           fill = "#2b8cbe") +
  theme_bw() +
  labs(y = expression(atop('Number of forays detected', '(x 1000)')), 
       x = '\n') + 
  theme(axis.title.x = element_text(family = "Arial", color = "black", size = 15, vjust = -1),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq( 50,350,100), labels = seq( 50,350,100))

# plot for sex effect
sex_sensi <- ggplot(data = tbl_df,
                         aes(x = vec_days, 
                             y = sex_coef)) +
  geom_errorbar(aes(ymin = sex_coef - se_sex_coef,
                    ymax = sex_coef + se_sex_coef),
                width = 0) +
  geom_point(size = 3.5,
             color = "black",
             fill = "#2b8cbe",
             shape = 21) +
  geom_hline(yintercept = 0,
             linetype = 3) +
  theme_bw() +
  labs(y = expression(atop("Effect of sex on prospecting", 
                           "rate (forays / day)")), 
       x = '\n') + 
  theme(axis.title.x = element_text(family = "Arial", color = "black", size = 15,  vjust = -1),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(-0.60,0.60), breaks = seq(-0.60, 0.60, 0.20), labels = round(seq(-0.60, 0.60, 0.20), digits = 1)) +
  scale_x_continuous(breaks = seq( 50,350,100), labels = seq( 50,350,100))

# plot for prov effect
prov_sensi <- ggplot(data = tbl_df,
                    aes(x = vec_days, 
                        y = prov_coef )) +
  geom_errorbar(aes(ymin = prov_coef - se_prov_coef,
                    ymax = prov_coef + se_prov_coef),
                width = 0) +
  geom_point(size = 3.5,
             color = "black",
             fill = "#2b8cbe",
             shape = 21) +
  geom_hline(yintercept = 0,
             linetype = 3) +
  theme_bw() +
  labs(y = expression(atop("Effect of cooperative provisioning", 
                           "on prospecting rate (forays / day)")), 
       x = '\n') + 
  theme(axis.title.x = element_text(family = "Arial", color = "black", size = 15,  vjust = -1),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(-0.40,0.40), 
                     breaks = seq(-0.40, 0.40, 0.20), 
                     labels = round(seq(-0.40, 0.40, 0.20), digits = 1)) +
  scale_x_continuous(breaks = seq( 50,350,100), labels = seq( 50,350,100))



Figure_sensi <- ggdraw() +
  draw_plot(n_forays, x = 0.00, y = 0.00, width = 0.33, height = 1.00) +
  draw_plot(sex_sensi, x = 0.33, y = 0.00, width = 0.33, height = 1.00) +
  draw_plot(prov_sensi, x = 0.66, y = 0.00, width = 0.33, height = 1.00) +
  draw_plot_label(label = c("A", "B", 'C'), size = 15, x = c(0.00, 0.33, 0.66), y = c(1, 1, 1)) +
  draw_label('Distance threshold used to assign prospecting forays', 
             colour = "black", 
             size = 15, angle = 0, x = 0.5, y = 0.05)


#####

ggsave(filename = "./plots/Figure_S2.png", 
       plot = Figure_sensi, 
       units = "mm",
       device = "jpeg", 
       width = 200,
       height = 125)


