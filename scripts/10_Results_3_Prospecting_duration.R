###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2024-05-17
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
data <- readRDS(file = "./data/09_prospecting_list_250m_15s.RDS")
data$sex_name <- ifelse(data$sex == "FEMALE", "female", "male")
data$age_years <- data$age_days/365
length(unique(data$bird_id)) # number of birds
nrow(data)
head(data)


##
## foray duration
data %>% 
  summarise(mean_time_next_home = mean(time_next_home)/60,
            median_time_next_home = median(time_next_home)/60, 
            minimum_time_next_home = min(time_next_home)/60,
            maximum_time_next_home = max(time_next_home)/60)
##
## n forays
nrow(data)

#####

##
##
##### Analysis foray duration #####
##
##

# full model
foray_duration_full_model <- lmer(log(time_next_home/60) ~
                                    sex_name : prov +
                                    sex_name +
                                    prov +
                                    age_years +
                                    (1|bird_id) + 
                                    (1|home),
                                  data = data,
                                  na.action = "na.fail")
# check residuals
hist(residuals(foray_duration_full_model)) # looks fairly good

drop1(foray_duration_full_model, test = "Chisq")
summary(foray_duration_full_model)

##
## model without interaction
foray_duration_model <- lmer(log(time_next_home/60) ~
                               sex_name +
                               prov +
                               age_years +
                               (1|bird_id) + 
                               (1|home),
                             data = data,
                             na.action = "na.fail")
drop1(foray_duration_model, test = "Chisq")
summary(foray_duration_model)


#####

##
##
##### Table of results #####
##
##

## base table
table_distance00 <- foray_duration_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `sex_name` = "Subordinate sex", 
                   `prov` = "Provisioning phase",
                   `age_years` = "Subordinate age"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_distance <- table_distance00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=foray_duration_model) %>% 
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
gtsave(table_distance, "./tables/TABLE S18 - Prospecting duration.html")

#####

##
##
##### Plot #####
##
##

##
## plot predictions
df_predict <- expand.grid(
  sex_name = c("female", "male"),
  prov = c("NO", "YES"),
  age_years = mean(data$age_days)/365
)

# model predictions
df_predict$fit <- predict(foray_duration_model, 
                          df_predict, 
                          re.form = NA, 
                          type = "link")

mm <- model.matrix(~sex_name +
                     prov +
                     age_years,
                   data = df_predict)

pvar1 <- diag(mm %*% tcrossprod(vcov(foray_duration_model),mm))
cmult <- 1 ## could use 1.96
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)
df_predict$fit_resp <- exp(df_predict$fit)
df_predict$plow_resp <- exp(df_predict$plo)
df_predict$phi_resp <- exp(df_predict$phi)
data_days$sex_name <- ifelse(data$sex == "FEMALE", "female", "male")

##
## individual prob of helping
duration_plot <- ggplot(df_predict,
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
  labs(y = expression(atop("Foray duration", "(minutes)")), 
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

saveRDS(file = "./plots/Figure 3c.RDS", 
        object = duration_plot)

ggsave(filename = "./plots/Figure_3c.png", 
       plot = duration_plot, 
       units = "mm",
       device = "png", 
       width = 89,
       height = 89)

#####

##
##
##### data histograms #####
##
##
table(data$time_next_home/60 >150, data$sex)

data_plot <- data %>% 
  mutate(duration_plot2 = ifelse(time_next_home/60 >150, 151, time_next_home/60))


durations_histogram <- ggplot(data_plot,
       aes(x = duration_plot2, fill = sex_name)) +
  geom_histogram(binwidth = 10) +
  facet_grid(sex_name~.) +
  theme_bw() +
  theme(axis.title.x = element_text(family = "Arial", color = "black", size = 15),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12, angle = -45, hjust = 0),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(family = "Arial", color = "black", size = 12)) +
  scale_fill_manual(values = c("#9970ab", "#5aae61")) +
  scale_x_continuous(labels = c('0','25', '50', '75', '100', '125', '> 150'), breaks = seq(0,150, 25)) +
  labs(y = 'Count of forays\n', x = expression(atop('Foray duration', '(minutes)')))
