###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2023/09/05
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
pacman::p_load(dplyr, lubridate, ggplot2, ggpattern,
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
head(data)


##
## foray distance
data %>% 
  summarise(mean_max_distance = mean(max_distance)*2,
            median_max_distance = median(max_distance)*2, 
            minimum_max_distance = min(max_distance)*2,
            maximum_max_distance = max(max_distance)*2)
##
## n forays
nrow(data)

##
## foray distance
data %>% 
  summarise(mean_max_distance = mean(max_distance),
            median_max_distance = median(max_distance), 
            minimum_max_distance = min(max_distance),
            maximum_max_distance = max(max_distance))

#####

##
##
##### Analysis foray max distance #####
##
##

# full model
foray_distance_full_model <- lmer(log(max_distance_2) ~
                                    sex_name : prov +
                                    sex_name +
                                    prov +
                                    age_years +
                                    (1|bird_id) + 
                                    (1|home),
                                  data = data %>% mutate(max_distance_2 = max_distance * 2),
                                  na.action = "na.fail")
# check residuals
hist(residuals(foray_distance_full_model))

drop1(foray_distance_full_model, test = "Chisq")
summary(foray_distance_full_model)



##
## model without interaction
foray_distance_model <- lmer(log(max_distance_2) ~
                               sex_name +
                               prov +
                               age_years +
                               (1|bird_id) + 
                               (1|home),
                             data = data %>% mutate(max_distance_2 = max_distance * 2),
                             na.action = "na.fail")
drop1(foray_distance_model, test = "Chisq")
summary(foray_distance_model)


#####

##
##
##### Table of results #####
##
##

## base table
table_distance00 <- foray_distance_model %>%
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
                        y = drop1_output(x=foray_distance_model) %>% 
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
gtsave(table_distance, "./tables/TABLE S19 - Prospecting distance.html")

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
df_predict$fit <- predict(foray_distance_model, 
                          df_predict, 
                          re.form = NA, 
                          type = "link")

mm <- model.matrix(~sex_name +
                     prov +
                     age_years,
                   data = df_predict)

pvar1 <- diag(mm %*% tcrossprod(vcov(foray_distance_model),mm))
cmult <- 1 ## could use 1.96
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)
df_predict$fit_resp <- exp(df_predict$fit)
df_predict$plow_resp <- exp(df_predict$plo)
df_predict$phi_resp <- exp(df_predict$phi)
data$sex_name <- ifelse(data$sex == "FEMALE", "female", "male")

##
## individual prob of helping
distance_plot <- ggplot(df_predict,
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
  labs(y = expression(atop("Foray distance", "(meters)")), 
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

saveRDS(file = "./plots/Figure 3d.RDS", 
        object = distance_plot)

ggsave(filename = "./plots/Figure_3d.png", 
       plot = distance_plot, 
       units = "mm",
       device = "png", 
       width = 89,
       height = 89)

#####

##
##
##### data histogram #####
##
##
distance_plot <- ggplot(data,
                         aes(x = max_distance*2, fill = sex_name)) +
  geom_histogram(binwidth = 100) +
  geom_ribbon_pattern(data = data.frame(x = c(0,500)), aes(x=x, ymin=0, ymax=Inf), 
                      pattern_color = NA,
                      pattern_fill = "grey25",
                      pattern = "stripe", 
                      pattern_alpha = 0.5,
                      fill = NA) +
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
  scale_x_continuous(labels = seq(500, 2100, 500), breaks = seq(500, 2100, 500), limit = c(0, 2100)) +
  labs(y = 'Count of forays\n', x = expression(atop('Foray distance', '(meters)')))


