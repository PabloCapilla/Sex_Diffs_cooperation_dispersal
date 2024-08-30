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

data <- readRDS("./data/06_dispersal_analysis_data.RDS")
sum(data$disap) # 2 birds that EM after becoming dominants
head(data)

##
## preparing age categorical variable
b_points <- c((365*1)+1, (365*2)+1, (365*3)+1)
data$age_cat <- as.factor(ifelse(data$age < b_points[1], 1,                   # < 1y
                             ifelse(data$age < b_points[2], 2,            # 1y < age < 2y
                                    ifelse(data$age < b_points[3], 3,     # 2y < age < 3y
                                           4))))                                # > 3y

unique(data$age_cat)
length(unique(data$bird_id))
length(unique(data$natal_group))

#####

##
##
##### 2-sample test of equal proportions #####
##
##

##
## total number of males and females
table({data %>% 
    group_by(bird_id) %>% 
    filter(row_number() == 1)}$sex)

##
## dispersing individuals
table({data %>% 
    filter(disap == 1)}$sex)

data %>% 
  filter(disap == 1) %>% 
  group_by(bird_id) %>%  
summarise(n_obs_bird = n()) %>% 
  arrange(desc(n_obs_bird))

##
## test
prop.test(x = c(52, 30), 
          n = c(178,160),
          alternative = "two.sided",
          correct = T)

#####

##
##
##### Model for dispersal prob after 1y of life #####
##
##
interaction_model_dispersal <- glmer(disap ~ 
                                         age_cat : sex +
                                         age_cat +
                                         sex +
                                         (1|bird_id) +
                                         (1|natal_group) +
                                         (1|start_life_season),
                                       data = data, 
                                       na.action = "na.fail",
                                       family = "binomial", 
                                       glmerControl(optimizer = "bobyqa"))
# check model residuals
model_residuals <- simulateResiduals(fittedModel = interaction_model_dispersal, n = 1000)
DHARMa::testResiduals(model_residuals) # all tests look fairly okay

# model results
summary(interaction_model_dispersal)
drop1(interaction_model_dispersal, test = "Chisq")

# LRT for the interaction
anova(interaction_model_dispersal, update(interaction_model_dispersal, 
                                            .~. - age_cat : sex),
      test = "Chisq")

##
## model without interaction
full_model_dispersal <- update(interaction_model_dispersal, 
                               .~. - age_cat : sex)
summary(full_model_dispersal)
drop1(full_model_dispersal, test = "Chisq")

#####

##
##
##### Table of results S1 #####
##
##

## base table
table_dispersal00 <- full_model_dispersal %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `age_cat` = "Subordinate age class", 
                   sex = "Subordinate sex"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_dispersal_model <- table_dispersal00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=full_model_dispersal) %>% 
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
gtsave(table_dispersal_model, "./tables/TABLE S1 - Dispersal.html")

##### 

##
##
##### plotting age trajectories #####
##
##
df_predict_sta <- expand.grid(
  sex = c("F", "M"),
  age_cat = as.character(1:4)
)

df_predict_sta$fit <- predict(full_model_dispersal, 
                              df_predict_sta, 
                              re.form = NA, 
                              type = "link")
mm <- model.matrix(~
                     age_cat +
                     sex,
                   data = df_predict_sta)

pvar1 <- diag(mm %*% tcrossprod(vcov(full_model_dispersal),mm))
cmult <- 1 ## could use 1.96
df_predict_sta <- data.frame(
  df_predict_sta
  , plo = df_predict_sta$fit-cmult*sqrt(pvar1)
  , phi = df_predict_sta$fit+cmult*sqrt(pvar1)
)
df_predict_sta$plow_resp <- boot::inv.logit(df_predict_sta$plo)
df_predict_sta$phi_resp <- boot::inv.logit(df_predict_sta$phi)
df_predict_sta$fit_resp <- boot::inv.logit(df_predict_sta$fit)

df_predict_sta$sex_name <- ifelse(df_predict_sta$sex == "F", "female", "male")
data$sex_name <- ifelse(data$sex == "F", "female", "male")


##
## individual prob of helping
dispersal_plot <- ggplot(df_predict_sta,
                         aes(y = fit_resp, 
                             x = age_cat,
                             color = sex_name, 
                             fill = sex_name)) +
  geom_point(data = data,
             aes(y = disap,
                 x = age_cat,
                 color = sex_name),
             position = position_jitterdodge(jitter.width = 0.15, 
                                             jitter.height = 0.02,
                                             dodge.width = 0.25),
             size = 2.5,
             alpha = 0.10) +
  geom_errorbar(aes(ymin = (plow_resp), ymax = (phi_resp)), 
                width = 0,
                color = "black",
                position = position_dodge(width = 0.25)) +
  geom_point(size = 5,
             shape = 21,
             color = "black",
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  labs(x = "Subordinate age (years)", 
       y = expression(atop("Probability of dispersal", "from natal territory"))) + 
  theme(axis.title.x = element_text(family = "Arial", color = "black", size = 15),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        legend.text = element_text(family = "Arial", size = 10),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#9970ab", "#5aae61")) +
  scale_color_manual(values = c("#9970ab", "#5aae61")) +
  scale_x_discrete(labels = c("< 1", "1 - 2", "2 - 3", "> 3"))

saveRDS(file = "./plots/Figure 1a.RDS", 
        object = dispersal_plot)

ggsave(filename = "./plots/Figure 1a.png", 
       plot = dispersal_plot, 
       units = "mm",
       device = "png", 
       width = 85,
       height = 115)

