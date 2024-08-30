###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2024-08-30
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

#####

##
##
##### data #####
##
##
data <- readRDS(file = "./data/01_provisioning_analysis_data.RDS")
head(data)

##
##
##### sample sizes #####
##
##
nrow(data) # observations
length(unique(data$date)) # video sessions
length(unique(data$ind_id)) # individuals
length(unique(data$ind_id[data$sex == 'M'])) # males
length(unique(data$ind_id[data$sex == 'F'])) # females
length(unique(data$clutch_id)) # breeding attempts
length(unique(data$group_id)) # groups

data %>% 
  group_by(clutch_id, date) %>% 
  filter(row_number() == 1) %>% 
  group_by(clutch_id) %>% 
  summarise(n_days = n()) %>% 
  summarise(min = min(n_days), 
            max = max(n_days),
            med = median(n_days))

  



table(data$season)
table(data$ba_day)

#####

##
##
###### Model 1 - Negative binomial regression #####
##
##
data$ind_id <- as.factor(data$ind_id)
data$group_id <- as.factor(data$group_id)
data$clutch_id <- as.factor(data$clutch_id)
data$season <- as.factor(data$season)
data$sex <- as.factor(data$sex)
data$age_cat <- as.factor(data$age_cat)
data$ba_day <- as.factor(data$ba_day)


## full model with interaction
interaction_model_cooperation <- glmmTMB(provisioning_feeds ~ 
                                           age_cat : sex +
                                           
                                           age_cat +
                                           sex + 
                                           ba_day +
                                           nestlings +
                                           offset(log(session_duration)) +
                                           (1|group_id) + 
                                           (1|ind_id) + 
                                           (1|clutch_id) + 
                                           (1|season), 
                                         data=data,
                                         REML = F,
                                         family='nbinom2',
                                         ziformula = ~ 1)
# zero-inflated negative binomial distribution accouinting for potential overdispersion, and zero inflation

# model results
summary(interaction_model_cooperation)  

# LRT for the interaction
anova(interaction_model_cooperation, update(interaction_model_cooperation, .~.-age_cat : sex), test = "Chisq")

##
## model without interaction
full_model_cooperation <- update(interaction_model_cooperation, .~.-age_cat : sex)
summary(full_model_cooperation)
drop1(full_model_cooperation, test = "Chisq")

##
##
##### Table of results S2 #####
##
##

## base table
table_cooperation00 <- full_model_cooperation %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   age_cat = "Age", 
                   sex = "Sex",
                   ba_day = "Brood age",
                   nestlings = "Brood size"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_cooperation <- table_cooperation00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=full_model_cooperation) %>% 
                          dplyr::select(variable = term, Chisq=statistic, df),
                        by = "variable")
    output$df <- ifelse(output$row_type == "label",  output$df, NA)
    output$Chisq <- ifelse(output$row_type == "label",  output$Chisq, NA)
    return(output)
  }) %>% 
  modify_fmt_fun(c(Chisq) ~ function(x) style_number(x, digits = 2)) %>%
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
gtsave(table_cooperation, "./tables/TABLE S2 - Provisioning.html")


##
##
##### Plot with model predictions #####
##
##

## session durations for predictions
mean_duration_df <- data %>%  # video session for predictions
  group_by(age_cat, sex, nestlings) %>% 
  summarise(session_duration = mean(session_duration))

## make prediction data frame
df_predict <- expand.grid(
  sex = c("F", "M"),
  ba_day = as.factor(6:12),
  age_cat = as.factor(1:4),
  nestlings = c(1),
  group_id = NA,
  ind_id = NA,
  clutch_id = NA,
  season = NA
)

# complete data fram for preditions
df_predict <- left_join(x = df_predict,
                        y = mean_duration_df,
                        by = c("age_cat", "sex", "nestlings"))

# predictions
df_predict$fit <-  {predict(object = full_model_cooperation,
                            newdata = df_predict,
                            newparams = NULL,
                            se.fit = T,
                            re.form = NA,
                            type = c("response"))}$fit
df_predict$fit_se <-  {predict(object = full_model_cooperation,
                               newdata = df_predict,
                               newparams = NULL,
                               se.fit = T,
                               re.form = NA,
                               type = c("response"))}$se.fit
# get data ready
df_predict_ave <- df_predict %>%
  group_by(sex, age_cat) %>%
  summarise(fit = mean(fit, na.rm = T),
            fit_se = mean(fit_se, na.rm = T))
df_predict$sex_name <- ifelse(df_predict$sex == "F", "female", "male")
df_predict_ave$sex_name <- ifelse(df_predict_ave$sex == "F", "female", "male")

##
## plot
plot_cooperation_age <- ggplot(df_predict_ave,
                               aes(x = age_cat, 
                                   fill = sex_name,
                                   y = fit)) +
  geom_point(data = data,
             aes(y = provisioning_feeds/session_duration,
                 x = age_cat,
                 color = sex_name,
                 fill = sex_name),
             position = position_jitterdodge(jitter.width = 0.15, 
                                             jitter.height = 0.07,
                                             dodge.width = 0.5),
             shape = 21,
             size = 2.5,
             alpha = 0.10) +
  geom_errorbar(aes(ymin = fit - fit_se, ymax = fit + fit_se), 
                width = 0,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 5, 
             shape = 21,
             #fill= "black",
             color= "black",
             position = position_dodge(width = 0.5)) + 
  theme_bw() +
  labs(x = "Subordinate age (years)", 
       y = expression(atop("Cooperative provisioning rate", "(feeds / hour)"))) + 
  theme(axis.title.x = element_text(family = "Arial", color = "black", size = 15),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = c(0.75,0.90),
        #legend.position = "top",#c(0.5,0.85),
        legend.text = element_text(family = "Arial", size = 15),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0,20)) +
  scale_fill_manual(values = c("#9970ab", "#5aae61")) +
  scale_color_manual(values = c("#9970ab", "#5aae61")) +
  scale_x_discrete(labels = c("< 1", "1 - 2", "2 - 3", "> 3"))

saveRDS(file = "./plots/Figure 1b.RDS", 
        object = plot_cooperation_age)

ggsave(filename = "./plots/Figure 1b.png", 
       plot = plot_cooperation_age, 
       units = "mm",
       device = "png", 
       width = 85,
       height = 115)

#####