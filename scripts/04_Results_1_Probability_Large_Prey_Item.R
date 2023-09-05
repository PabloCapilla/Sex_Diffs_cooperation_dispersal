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

#####

##
##
##### data #####
##
##
data <- readRDS(file = "./data/03_prey_type_analysis_data.RDS")

##
## observations
nrow(data)

##
## number of video sessions
nrow(data %>%
       group_by(clutch_id, date) %>%
       count())
##
## birds
length(unique(data$ind_id))
length(unique(data$ind_id[data$sex == 'M']))
length(unique(data$ind_id[data$sex == 'F']))

##
## clutches
length(unique(data$clutch_id))
##
## groups
length(unique(data$GROUP))

table(is.na(data$nestlings))
table(is.na(data$nestlings))

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
# check model residuals
model_residuals <- simulateResiduals(fittedModel = interaction_model_prey_item, n = 1000)
DHARMa::testResiduals(model_residuals) # all tests look fairly okay

# model results
summary(interaction_model_prey_item)
anova(interaction_model_prey_item, update(interaction_model_prey_item, .~.-age_cat : sex), test = "Chisq")

##
## model without interaction
full_model_prey_item <- update(interaction_model_prey_item, .~.-age_cat : sex)
summary(full_model_prey_item)
drop1(full_model_prey_item, test = "Chisq")

#####

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

#####

##
##
##### Plot with model predictions #####
##
##

## make prediction data frame
df_predict <- expand.grid(
  sex = c("F", "M"),
  ba_day = as.factor(6:12),
  age_cat = as.factor(1:4),
  nestlings = c(1:3),
  group_id = NA,
  ind_id = NA,
  clutch_id = NA,
  SEASON = NA
)

# predictions
df_predict$fit <-  predict(object = full_model_prey_item,
                           newdata = df_predict,
                           newparams = NULL,
                           re.form = NA,
                           type = c("response"))

mm <- model.matrix(~
                     age_cat +
                     ba_day +
                     sex + 
                     nestlings,
                   data = df_predict)

pvar1 <- diag(mm %*% tcrossprod(vcov(full_model_prey_item),mm))
cmult <- 1 
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)
df_predict$plow_resp <- boot::inv.logit(df_predict$plo)
df_predict$phi_resp <- boot::inv.logit(df_predict$phi)
df_predict$fit_resp <- boot::inv.logit(df_predict$fit)
df_predict$sex_name <- ifelse(df_predict$sex == "F", "female", "male")
data$sex_name <- ifelse(data$sex == "F", "female", "male")

# get data ready
df_predict_ave <- df_predict %>%
  group_by(sex) %>%
  summarise(fit = mean(fit_resp, na.rm = T),
            fit_lo = mean(plow_resp, na.rm = T),
            fit_hi = mean(phi_resp, na.rm = T))
df_predict_ave$sex_name <- ifelse(df_predict_ave$sex == "F", "female", "male")


##
## plot
plot_large_item <- ggplot(df_predict_ave,
                             aes(x = sex_name, 
                                 fill = sex_name,
                                 y = fit)) +
  geom_point(data = data,
             aes(y = prey_type,
                 x = sex_name,
                 color = sex_name,
                 fill = sex_name),
             position = position_jitter(height = 0.05, width = 0.2),
             shape = 21,
             size = 2.5,
             alpha = 0.10) +
  geom_errorbar(data = df_predict_ave, 
                aes(ymin = fit_lo, ymax = fit_hi), 
                width = 0) +
  geom_point(size = 5, 
             shape = 21,
             #fill= "black",
             color= "black") + 
  theme_bw() +
  labs(x = "Subordinate sex", 
       y = expression(atop("Probability of delivering", "a large prey item"))) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",#c(0.5,0.85),
        legend.text = element_text(family = "Arial", size = 15),
        legend.title = element_blank()) +
  scale_y_continuous(breaks = c(0,0.5,1), labels = c(0,0.5,1)) +
  scale_fill_manual(values = c("#9970ab", "#5aae61")) +
  scale_color_manual(values = c("#9970ab", "#5aae61")) 

saveRDS(file = "./plots/Figure 1d.RDS", 
        object = plot_large_item)

ggsave(filename = "./plots/Figure 1d.png", 
       plot = plot_large_item, 
       units = "mm",
       device = "png", 
       width = 56,
       height = 100)

#####







