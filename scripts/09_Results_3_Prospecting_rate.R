###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2023-09-04
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
##### data #####
##
##
data <- readRDS(file = "./data/08_prospecting_rate_250m_15s.RDS")
head(data)

length(unique(data$bird_id))

## initial forays
sum(data$forays)

#####

##
##
##### Add weight of tagged ind #####
##
##
captures <- read.csv('../prep_repository/raw_data/All CAPTURES_adults_Maria.csv')
captures <- captures %>% 
  select(Date, BIRD.ID, Weight) %>% 
  mutate(Date = dmy(Date))
head(captures)

captures_tagged <- captures[captures$BIRD.ID %in% unique(data$bird_id),]
captures_tagged %>% 
  group_by(BIRD.ID) %>% 
  summarise(mean_weight = mean(Weight)) %>% 
  filter(!is.na(mean_weight)) %>% 
  mutate(prop_weight = 1.2/mean_weight) %>% 
  summarise(mean_prop_weight = mean(prop_weight),
            min_prop_weight = min(prop_weight),
            max_prop_weight = max(prop_weight))


#####

##
##
##### Summary info and plot #####
##
##

##
## 24h plot of rates
data %>% 
  group_by(bird_id, hour) %>%
  summarise(mean_foray_hour = mean(forays)) %>%
  ggplot(aes(x = hour, y = mean_foray_hour)) +
  geom_point(fill = "grey60") 
  

data %>% 
  group_by(bird_id, hour) %>%
  summarise(mean_foray_hour = mean(forays)) %>% 
  group_by(hour) %>% 
  summarise(n_points_hour = n())

# params to plot data
upper_lim <-0.15
rise <- 6
set <- 19

##
## plot 
forays_hour <- data %>% 
  group_by(bird_id, hour) %>%
  summarise(total_foray_hour = mean(forays)) %>%
  group_by(hour) %>%
  summarise(mean_rate = mean(total_foray_hour),
            se_rate = sd(total_foray_hour)/sqrt(n()),
            n = n()) %>%
  ggplot(aes(x = hour, y = mean_rate)) +
  geom_col(fill = "grey60") +
  geom_errorbar(aes(ymin = mean_rate-se_rate, ymax = mean_rate+se_rate),
                width = 0) +
  geom_polygon(data = data.frame(x = c(0.5,0.5,rise,rise,set,set,24.5,24.5), 
                                 y = c(0,upper_lim,upper_lim,0,0,upper_lim,upper_lim,0)),
               aes(x = x, y = y), 
               fill = "#3182bd",
               alpha = 0.25) +
  scale_x_continuous(limits = c(0,24.5), breaks = seq(2,24,4), labels = seq(2,24,4)) + 
  scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.15,0.05), labels = seq(0,0.15,0.05)) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(family = "Arial", color = "black", size = 15),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12), 
        plot.title = element_text("Arial", size = 6)) +
  labs(x = "\nHour of the day", y = expression(atop('Prospecting rate', '(forays / bird / hour)')))


saveRDS(file = "./plots/Figure 3b.RDS", 
        object = forays_hour)

ggsave(filename = "./plots/Figure_3b.png",
       plot = forays_hour, 
       units = "mm", 
       height = 89, 
       width = 89, 
       device = "png")

#####

##
##
##### Summarising data for the day #####
##
##
data_days <- data %>% 
  group_by(bird_id, day, sex, prov, home) %>%
  summarise(total_foray_day = sum(forays), n_days = n(),
            age_days = mean(age_days))
data_days$sex_name <- ifelse(data_days$sex == "FEMALE", "female", "male")
data_days$obs_id <- 1:nrow(data_days)
data_days$age_years <- data_days$age_days/365

##
## number of forays
sum(data_days$total_foray_day)

# binomial forays
data_days$foray_cat <- ifelse(data_days$total_foray_day != 0, 1, 0)
head(data_days)

#####

##
##
##### Summary of raw rates of prospecting #####
##
##
range(data_days$total_foray_day)
nrow(data_days)
length(unique(data_days$bird_id))
mean(data_days$total_foray_day)
median(data_days$total_foray_day) 

#####

##
##
##### Analysis for proportion of days with any foray #####
##
##

##
## age diff between sexes
summary(lm(age ~ sex_name, 
           data = data_days %>% 
             group_by(bird_id, sex_name) %>% 
             summarise(age = mean(age_days))))


# full model
prospecting_full_model <- glmer(total_foray_day ~
                                  sex_name : prov +
                                  sex_name +
                                  prov +
                                  age_years +
                                  (1|obs_id) + 
                                  (1|bird_id) + 
                                  (1|home),
                                data = data_days,
                                family = "poisson",
                                na.action = "na.fail")
drop1(prospecting_full_model, test = "Chisq")
summary(prospecting_full_model)

library(DHARMa)
res_poisson <- simulateResiduals(fittedModel = prospecting_full_model, n = 500)
testUniformity(res_poisson)
testZeroInflation(res_poisson)


##
## model without interaction
prospecting_model <- glmer(total_foray_day ~
                             sex_name +
                             prov +
                             age_years +
                             (1|obs_id) +
                             (1|bird_id) + 
                             (1|home),
                           data = data_days,
                           family = "poisson",
                           na.action = "na.fail")
drop1(prospecting_model, test = "Chisq")
summary(prospecting_model)

#####

##
##
##### Table of results #####
##
##

## base table
table_prospecting00 <- prospecting_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `sex_name` = "Subordinate sex", 
                   `prov` = "Provisioning phase",
                   `age_years` = "Subordinate age"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_prospecting <- table_prospecting00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=prospecting_model) %>% 
                          dplyr::select(variable = term, Chisq=statistic, df),
                        by = "variable")
    output$df <- ifelse(output$row_type == "label",  output$df, NA)
    output$Chisq <- ifelse(output$row_type == "label",  output$Chisq, NA)
    return(output)
  }) %>% 
  modify_fmt_fun(c(Chisq) ~ function(x) style_number(x, digits = 2)) %>%
  modify_fmt_fun(c(std.error) ~ function(x) style_number(x, digits = 3)) %>%
  modify_fmt_fun(c(p.value) ~ function(x) style_number(x, digits = 3)) %>%
  modify_table_body(~.x %>% dplyr::relocate(p.value, .after = df)) %>% 
  modify_header(label ~ "**Fixed effect**") %>% 
  modify_header(std.error ~ "**SE**") %>%
  modify_header(estimate ~ "**Estimate**") %>%
  modify_header(df ~ "**df**") %>% 
  modify_header(Chisq ~ html("<b>&chi;<sup>2</sup></b>")) %>% 
  as_gt() %>% 
  opt_footnote_marks(marks = "LETTERS")

##
## save table
gtsave(table_prospecting, "./tables/TABLE S17 - Prospecting rate.html")

#####

##
##
##### Plot with model predictions #####
##
##

##
## refit model to facilitate generating predictions
prospecting_model_plot <- glmer(total_foray_day ~
                                  sex_name +
                                  prov +
                                  age_years +
                                  (1|obs_id) +
                                  (1|bird_id) + 
                                  (1|home),
                                data = data_days,
                                family = "poisson",
                                na.action = "na.fail")
##
## plot predictions
df_predict <- expand.grid(
  sex_name = c("female", "male"),
  prov = c("NO", "YES"),
  age_years = mean(data_days$age_days)/365
)

# model predictions
df_predict$fit <- predict(prospecting_model_plot, 
                          df_predict, 
                          re.form = NA, 
                          type = "link")

mm <- model.matrix(~sex_name +
                     prov +
                     age_years,
                   data = df_predict)

pvar1 <- diag(mm %*% tcrossprod(vcov(prospecting_model_plot),mm))
cmult <- 1 ## could use 1.96
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)
df_predict$fit_resp <- exp(df_predict$fit)
df_predict$plow_resp <- exp(df_predict$plo)
df_predict$phi_resp <- exp(df_predict$phi)
data_days$sex_name <- ifelse(data_days$sex == "FEMALE", "female", "male")

##
## individual prob of helping
prospecting_plot <- ggplot(df_predict,
                 aes(x = prov, 
                     fill = sex_name,
                     color = sex_name,
                     y = fit_resp)) + 
  geom_errorbar(aes(ymin = (plow_resp), 
                    ymax = (phi_resp), 
                    group = sex_name), 
                width = 0,
                color= "black",
                position = position_dodge(width = 0.5)) +
  geom_point(size = 5, 
             shape = 21,
             color= "black",
             position = position_dodge(width = 0.5)) + 
  theme_bw() +
  labs(y = expression(atop("Prospecting rate", "(forays / day)")), 
       x = "\nProvisioning") + 
  theme(axis.title.x = element_text(family = "Arial", color = "black", size = 15),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#9970ab", "#5aae61")) +
  scale_color_manual(values = c("#9970ab", "#5aae61"), name = c("", ""))

saveRDS(file = "./plots/Figure 3e.RDS", 
        object = prospecting_plot)

ggsave(filename = "./plots/Figure_3e.png", 
       plot = prospecting_plot, 
       units = "mm",
       device = "png", 
       width = 89,
       height = 89)

