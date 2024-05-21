###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2024-05-21
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
data <- readRDS(file = "./data/10_prospecting_cooperation_individual.RDS")
head(data)

# sample sizes
length(unique(data$date))
length(unique(data$group_id))


#####

## n of birds tracked
length(unique(data$bird_id))
length(unique(data$bird_id[data$sex == 'M']))
length(unique(data$bird_id[data$sex == 'F']))

##
##
##### Models #####
##
##

## model 1: rate of morning forays continuous
tradeoff_full_model <- glmer(provisioning_feeds ~ 
                               foray_rate_before : sex +
                               foray_rate_before +
                               nestlings +
                               ba_day + 
                               sex +
                               (1|group_id) + 
                               (1|ind_id),
                             data = data,
                             offset = log(data$session_duration),
                             family = 'poisson',
                             glmerControl(optimizer = 'bobyqa'),
                             na.action = "na.fail")
# inspecting model residuals
simul_res <- simulateResiduals(fittedModel = tradeoff_full_model, n = 1000)
testResiduals(simul_res)

summary(tradeoff_full_model)
drop1(tradeoff_full_model, test = "Chisq")


## model without interaction

tradeoff_model_poisson <- glmer(provisioning_feeds ~ 
                                  foray_rate_before +
                                  nestlings +
                                  ba_day + 
                                  sex +
                                  (1|group_id) + 
                                  (1|ind_id),
                                data = data,
                                offset = log(data$session_duration),
                                glmerControl(optimizer = 'bobyqa'),
                                family = 'poisson',
                                na.action = "na.fail")
drop1(tradeoff_model_poisson, test = "Chisq")
summary(tradeoff_model_poisson)

plot(tradeoff_model_poisson)
hist(residuals(tradeoff_model_poisson))

#####

##
##
##### Table of results #####
##
##


## base table
table_tradeoff00 <- tradeoff_model_poisson %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `sex` = "Subordinate sex", 
                   `foray_rate_before` = "Prospecting rate (forays/day)",
                   `ba_day` = "Brood age",
                   `nestlings` = "Brood size"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_tradeoff <- table_tradeoff00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=tradeoff_model_poisson) %>% 
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
gtsave(table_tradeoff, "./tables/TABLE S11 - Trade off.html")

#####


##
##
##### Plot model predictions #####
##
##

# model predictions
df_predict <- expand.grid(sex = c("F", "M"),
                          ba_day = mean(data$ba_day), 
                          nestlings = mean(data$nestlings),
                          foray_rate_before = seq(min(data$foray_rate_before),
                                      max(data$foray_rate_before),
                                      0.01))
df_predict$session_duration <- median(data$session_duration)
df_predict$fit <- predict(object = tradeoff_model_poisson, 
                          df_predict, 
                          re.form = NA, 
                          type = "link")


mm <- model.matrix(delete.response(terms(tradeoff_model_poisson)),df_predict)
predvar <- diag(mm %*% vcov(tradeoff_model_poisson) %*% t(mm))
df_predict$fit_se <- sqrt(predvar)
df_predict$low_resp <- exp(df_predict$fit - df_predict$fit_se)
df_predict$high_resp <- exp(df_predict$fit + df_predict$fit_se)
df_predict$fit_resp <- exp(df_predict$fit)

# cut predictions to data range
#dfmales <- df_predict[df_predict$sex == "M" & df_predict$foray_rate_before <= max(data$foray_rate_before[data$sex == "M"]),] 
#dffemales <- df_predict[df_predict$sex == "F" & df_predict$foray_rate_before <= max(data$foray_rate_before[data$sex == "F"]),] 
#df_predict <- rbind(dfmales, dffemales)

df_predict <- df_predict %>% 
  group_by(foray_rate_before) %>% 
  summarise(fit = mean(fit_resp),
            low_resp = mean(low_resp),
            high_resp = mean(high_resp))

## plot
tradeoff_plot <- ggplot(df_predict,
                           aes(y = fit,
                               x = foray_rate_before)) +
  geom_point(data = data,
             aes(y = provisioning_rate, x = foray_rate_before),
             shape =21,
             size = 5, 
             alpha =1,
             color = 'black',
             fill = 'grey75',
             position = position_jitter(width = 0.05)) +
  geom_line(size = 1.5) + 
  geom_ribbon(aes(ymin = low_resp, 
                  ymax = high_resp),
              color = NA,
              alpha = 0.25) +
  theme_bw() +
  labs(x = "\nProspecting rate (forays / day)", 
       y = expression(atop("Cooperative provisioning rate", "(feeds / hour)"))) + 
  theme(axis.title.x = element_text(family = "Arial", color = "black", size = 15),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(labels = seq(0,3,0.5), breaks = seq(0,3,0.5))

saveRDS(file = "./plots/Figure 3f.RDS", 
        object = tradeoff_plot)

ggsave(filename = "./plots/Figure_3f.png", 
       plot = tradeoff_plot, 
       units = "mm",
       device = "png", 
       width = 89,
       height = 89)



#prospecting_prov_margu <- ggMarginal(p = tradeoff_plot, 
#                                     type="boxplot", 
#                                     groupFill = T, 
#                                     groupColour = T,
#                                     size = 7.5)
#
#
#ggsave(filename = "./plots/Figure_2d.png", 
#       plot = prospecting_prov_margu, 
#       units = "mm",
#       device = "png", 
#       width = 89,
#       height = 89)


