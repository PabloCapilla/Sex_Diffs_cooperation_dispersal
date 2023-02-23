###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 10/02/2023
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
data <- readRDS(file = "./data/04_prey_type_analysis_data.RDS")

##
## number of video sessions
nrow(data %>%
       group_by(clutch_id, date) %>%
       count())
##
## birds
length(unique(data$ind_id))
##
## clutches
length(unique(data$clutch_id))
##
## groups
length(unique(data$group_id))


##
##
###### Model size of prey item #####
##
data$ind_id <- as.factor(data$ind_id)
data$group_id <- as.factor(data$group_id)
data$clutch_id <- as.factor(data$clutch_id)
data$season <- as.factor(data$SEASON)
data$sex <- as.factor(data$sex)
data$age_cat <- as.factor(data$age_cat)
data$ba_day <- as.factor(data$ba_day)


interaction_model_prey_item <- glmer(prey_type ~ 
                                       age_cat : sex +
                                       age_cat +
                                       ba_day +
                                       nestlings+
                                       sex + 
                                       (1|group_id) + 
                                       (1|ind_id) + 
                                       (1|SEASON) +
                                       (1|clutch_id),
                                     glmerControl(optimizer = "bobyqa"),
                                     family = "binomial",
                                     na.action = "na.fail",
                                     data = data)
summary(interaction_model_prey_item)

# LRT for the interaction
anova(interaction_model_prey_item, update(interaction_model_prey_item, .~.-age_cat : sex), test = "Chisq")

##
## model without interaction
full_model_prey_item <- update(interaction_model_prey_item, .~.-age_cat : sex)
summary(full_model_prey_item)
check_model(full_model_prey_item)
drop1(full_model_prey_item, test = "Chisq")


##
##
##### Table of results S4 #####
##
##


## base table
table_prey_type00 <- full_model_prey_item %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   age_cat = "Age", 
                   sex = "Sex",
                   ba_day = "Brood age",
                   nestlings = "Brood size"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_prey_type <- table_prey_type00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=full_model_prey_item, family = "binomial") %>% 
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
gtsave(table_prey_type, "./tables/TABLE S4 - Prey Type.html")



