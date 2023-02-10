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
data <- readRDS(file = "data/02_weaving_analysis_data.RDS")
head(data)
table(data$grass_cat)

##
##
###### Model 1 - Probability of  #####
##
m1_fit <- glmer(grass_cat ~ 
                  as.factor(age_cat): sex_num +
                  
                  as.factor(age_cat) +
                  as.factor(ba_day) +
                  sex_num + 
                  nestlings +
                  (1|group_id) + 
                  (1|ind_id) + 
                  (1|clutch_id) + 
                  (1|season), 
                data = data, 
                family = "binomial",
                na.action = "na.fail",
                glmerControl(optimizer = "bobyqa"))
summary(m1_fit)

## assessing fit
RVAideMemoire::overdisp.glmer(m1_fit) # evidence for overdispersion
m1_fit_res <- simulateResiduals(fittedModel = m1_fit, 
                                 n = 1000)
testUniformity(m1_fit_res)    # again evidence for overdispersion


ms_table <- dredge(m1_fit, rank = "AIC", trace = 3)
subset(ms_table, delta < 6)
subset(ms_table, delta < 6 & !nested(.))

top_model <- get.models(ms_table, 1)[[1]]



##
##
##### Plot for weaving probability #####
##
##
mean_calc <- data_model %>% 
  group_by(sex) %>% 
  summarise(mean_prob = mean(grass_cat), 
            se_prob = sd(grass_cat)/sqrt(n()))


weaving_prob_plot <- ggplot(data = data_model, aes(x = sex, y = grass_cat, color = sex)) +
  geom_point(position = position_jitter(width = 0.05, height = 0.025),
             shape = 19,
             alpha = 0.15,
             size = 0.5) +
  geom_point(data = mean_calc, 
             aes(y = mean_prob, x = sex, color = sex),
             size = 1.5, 
             shape = 17) +
  theme_bw() +
  theme(axis.title.x = element_text(family = "Arial", size = 10),
        axis.title.y = element_text(family = "Arial", size = 10),
        axis.text.x = element_text(family = "Arial", size = 8),
        axis.text.y = element_text(family = "Arial", size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(family = "Arial", size = 8),
        legend.title = element_blank()) +
  scale_color_manual(values = c("#9970ab", "#5aae61")) +
  scale_x_discrete(labels = c("female", "male")) +
  labs(x = "Sex of the individual", y = "Probability of weaving")


ggsave(filename = "./plots/weaving_prob.jpeg", 
       plot = weaving_prob_plot, 
       units = "mm",
       device = "jpeg", 
       width = 89,
       height = 89)

## no raw data and model predict
df_predict_ini <- expand.grid(
  sex = c("F", "M"),
  ba_day = 10
)

m1_sta_top <- standardize::standardize(prov_cat ~ 
                                     sex +
                                     ba_day +
                                     (1|group_id) + 
                                     (1|ind_id) + 
                                     (1|clutch_id) + 
                                     (1|season), 
                                   family = "binomial",
                                   data = data_model)
df_predict <- predict(m1_sta_top, df_predict_ini, random = FALSE)


df_predict$fit <- predict(top_model, 
                          df_predict, 
                          re.form = NA, 
                          type = "link")

mm <- model.matrix(~ 
                     sex + 
                     ba_day,
                   data = df_predict)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(top_model),mm))
cmult <- 1 ## could use 1.96
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)
df_predict_ini$fit_resp <- boot::inv.logit(df_predict$fit)
df_predict_ini$plow_resp <- boot::inv.logit(df_predict$plo)
df_predict_ini$phi_resp <- boot::inv.logit(df_predict$phi)

df_predict_ini$sex_name <- ifelse(df_predict_ini$sex == "F", "females", "males")
data_model$sex_name <- ifelse(data_model$sex == "F", "females", "males")

weaving_plot <- ggplot(data = df_predict_ini, aes(x = sex, y = fit_resp, fill = sex)) +
  geom_errorbar(aes(ymin = plow_resp, ymax = phi_resp), 
                width = 0) +
  geom_point(size = 2.5, 
             shape = 21) +
  theme_bw() +
  theme(axis.title.x = element_text(family = "Arial", size = 10),
        axis.title.y = element_text(family = "Arial", size = 10),
        axis.text.x = element_text(family = "Arial", size = 8),
        axis.text.y = element_text(family = "Arial", size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(family = "Arial", size = 8),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#9970ab", "#5aae61")) +
  scale_y_continuous(limits = c(0, 0.035), labels = seq(0,0.04, 0.01)) +
  scale_x_discrete(labels = c("females", "males")) +
  labs(x = "Sex of the individual", y = "Probability of weaving")

ggsave(filename = "./plots/weaving_prob_pred.jpeg", 
       plot = weaving_plot, 
       units = "mm",
       device = "jpeg", 
       width = 89,
       height = 89)
## 
## 
##### number of feeds when more than 1 #####
##
##

# data with grass
sub_data <- data_model %>% filter(prov_cat == 1)

number_grass <- lmer(final_feeds ~ 
           age_cat : sex +
           
           age_cat +
           sex +
           ba_day +
           nestlings +
           (1|group_id) + 
           (1|ind_id) + 
           (1|clutch_id) + 
           (1|season), 
         data = sub_data)
summary(number_grass)


t_test <- lm(final_feeds ~ 
       sex,
     data = sub_data)
summary(t_test)


ggplot(data = sub_data, aes(x = sex, y = final_feeds)) +
  geom_point(position = position_jitter(width = 0.1))


