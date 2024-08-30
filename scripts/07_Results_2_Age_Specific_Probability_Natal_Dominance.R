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

## full dispersal events
data <- readRDS("./data/06_dispersal_full_table_analysis_data.RDS")
head(data)

## budding data
data_budding <- read.csv(file = "./data/07_dominance_via_budding.csv") # data

# dominant in natal terr?
data$dom_in_natal <- ifelse(is.na(data$dominant_territory), 0, 
                            ifelse(data$dominant_bird == 1 &
                                     data$dominant_territory == data$Start.Group, 1, 0))

#####

##
##
##### Data summaries #####
##
##
table(data$Fate)
table(data$dominant_bird, data$dominant_territory == data$Start.Group, data$SEX)

# how many birds remained home at age 2.75 per sex? (median age at dom acquisition)
table({data %>% 
    filter(age_fate > 365) %>%
    filter(dmy(BBBirthday) < "2014-06-01") %>%
    filter(dmy(LSeen.Date) < "2016-06-01")}$SEX)

# how many birds [older than 1y] remained home at age 2.75 per sex? 
table({data %>%
    filter(age_fate > 365/2) %>%
    filter(dmy(BBBirthday) < "2014-06-01") %>%
    filter(dmy(LSeen.Date) < "2016-06-01") %>%
    filter(age_fate > (2.75*365))}$SEX) 

#####

##
##
##### Dominance turnover in natal group + including budding: Individual probabilities #####
##
##

##
## data
natal_turnover_ind <- data %>%
  filter(SEX != "U") %>% 
  mutate(SEX = as.character(SEX)) %>% 
  mutate(dom_age = ymd(dominant_start) - dmy(BBBirthday))

length(natal_turnover_ind$BIRD.ID[natal_turnover_ind$SEX == 'F'])
length(natal_turnover_ind$BIRD.ID[natal_turnover_ind$SEX == 'M'])

#
## I create a loop to create a dataset for birds that are at least 6-month, 1-y, 2-y and 3-y old at the time of assessment
days_old <- c(365, 365*2, 365*3, 365*4)

list_predictions <- as.list(NA)
list_df <- as.list(NA)

for (d in 1:length(days_old)){
  ##
  ## data
  natal_turnover_ind_df_model <- natal_turnover_ind %>% 
    filter(ymd("2016-05-30") - (dmy(BBBirthday) + days_old[d]) > (365*2)) %>% # only include ind when there are two years after age limit
    filter(age_fate > days_old[d]) %>%  # birds that are at least d days
    mutate(age_cat = as.character(d))
  
  list_df[[d]] <- natal_turnover_ind_df_model
}
df <- rbindlist(list_df)

length(unique(df$BIRD.ID[df$SEX == 'F']))
length(unique(df$BIRD.ID[df$SEX == 'M']))

#####

##
##
##### model #####
##
##
natal_turnover_ind_model <- glmer(dom_in_natal ~
                                      age_cat:SEX +
                                      age_cat +
                                      SEX +
                                      (1|GROUP) +
                                      (1|Season),
                                    data = df, 
                                    na.action = "na.fail",
                                    family = "binomial")
summary(natal_turnover_ind_model)
drop1(natal_turnover_ind_model, test = 'Chisq')

## model withouth interaction
natal_turnover_ind_model_simple <- update(natal_turnover_ind_model, .~.-age_cat:SEX)
drop1(natal_turnover_ind_model_simple, test = 'Chisq') # model without interaction
summary(natal_turnover_ind_model_simple)

# check model fit
simul_res <- DHARMa::simulateResiduals(fittedModel = natal_turnover_ind_model, n = 1000)
testResiduals(simul_res) # they look fine
  
#####

##
##
##### table of results #####
##
##

## base table
table_natal_turn_over00 <- natal_turnover_ind_model_simple %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `age_cat` = 'Subordinate age',
                   `SEX` = "Sex"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_natal_turn_over <- table_natal_turn_over00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=natal_turnover_ind_model_simple) %>% 
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
table_name <- paste0("./tables/TABLE S6 - Dominance acquisition natal.html")
gtsave(table_natal_turn_over, table_name)

#####  

##
##
##### plot #####
##
##

# model predictions
df_predict <- expand.grid(SEX = c("M", "F"), 
                          age_cat = factor(1:4, ordered = T))
df_predict$fit <- predict(object = natal_turnover_ind_model_simple, 
                          df_predict, re.form = NA, type = "link")

mm <- model.matrix(~ SEX + age_cat,
                   data = df_predict)

pvar1 <- diag(mm %*% tcrossprod(vcov(natal_turnover_ind_model_simple),mm))
cmult <- 1 ## could use 1.96
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)
df_predict$plow_resp <- boot::inv.logit(df_predict$plo)
df_predict$phi_resp <- boot::inv.logit(df_predict$phi)
df_predict$fit_resp <- boot::inv.logit(df_predict$fit)

df_predict$sex_name <- ifelse(df_predict$SEX == "F", "female", "male")
natal_turnover_ind_df_model$sex_name <- ifelse(natal_turnover_ind_df_model$SEX == "F", "female", "male")
natal_turnover_ind_df_model$age_cat <- days_old[d]
df$age_cat <- factor(df$age_cat, ordered = T)



df$sex_name <- ifelse(df$SEX == "F", "female", "male")


dom_acq_natal <- ggplot(df, 
                        aes(x = age_cat, 
                            y = dom_in_natal, 
                            color = sex_name,
                            fill = sex_name)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, 
                                             jitter.height = 0.02,
                                             dodge.width = 0.25),
             shape = 21,
             size = 2.5,
             alpha = 0.10) +
  geom_errorbar(data = df_predict, 
                aes(ymin = plow_resp, 
                    ymax = phi_resp,
                    y = fit_resp, 
                    x = age_cat),
                position = position_dodge(width = 0.5),
                color = "black",
                width = 0) +
  geom_point(data = df_predict,
             aes(y = fit_resp, 
                 fill = sex_name),
             size = 5,
             position = position_dodge(width = 0.50),
             color = "black",
             shape = 21) +
  theme_bw() +
  labs(x = "Subordinate age (years)", 
       y = "Prob of dominance acquisition natal group") + 
  theme(axis.title.x = element_text(family = "Arial", size = 15),
        axis.title.y = element_text(family = "Arial", size = 15),
        axis.text.x = element_text(family = "Arial", size = 12),
        axis.text.y = element_text(family = "Arial", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        legend.text = element_text(family = "Arial", size = 10),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#9970ab", "#5aae61")) +
  scale_color_manual(values = c("#9970ab", "#5aae61")) +
  scale_x_discrete(labels = c("1", "2", "3", "4"))

ggsave(filename = "./plots/Figure 2.png", 
       plot = dom_acq_natal, 
       units = "mm",
       device = "png", 
       width = 75,
       height = 120)

#####


##
##### 2 -  Natal turn over + budding #####
##

##
## changing budded males to dominance acq in turnover table

# males that budded to become dominants
budded_males <- {data_budding %>% 
    filter(budded_dom_male == 1) %>% 
    select(X1st_dom_male)}$X1st_dom_male

# females that budded to become dominants
budded_females <- {data_budding %>% 
    filter(budded_dom_female == 1) %>% 
    select(X1st_dom_female)}$X1st_dom_female

# are they as 'no natal dominant' in turnover table?
table(budded_males %in% natal_turnover_ind$BIRD.ID) # 1 male budded before year 1
table(natal_turnover_ind$BIRD.ID %in% budded_males) # 5/6 males present in table as expected

table(budded_females %in% natal_turnover_ind$BIRD.ID) # all there
table(natal_turnover_ind$BIRD.ID %in% budded_females) # all females present in table as expected


# new data including budded groups as 'inheritance'
all_inh_ind <- natal_turnover_ind
all_inh_ind$dom_in_natal[all_inh_ind$BIRD.ID %in% budded_males] <- 1
all_inh_ind$dom_in_natal[all_inh_ind$BIRD.ID %in% budded_females] <- 1

# compare new tables against natal dom table
table(natal_turnover_ind$dom_in_natal)
table(all_inh_ind$dom_in_natal) # fine, 14 new 'natal dominants'

days_old <- c(365, 365*2, 365*3, 365*4)

list_predictions_bud <- as.list(NA)
list_df_bud <- as.list(NA)

for (d in 1:length(days_old)){
  
  ##
  ## data
  natal_turnover_ind_bud <- all_inh_ind %>%
    filter(ymd("2016-05-30") - (dmy(BBBirthday) + days_old[d]) > (365*2)) %>% # only include ind when there are two years after age limit
    filter(age_fate > days_old[d]) %>%  # birds that are at least d days
    mutate(age_cat = d)  
  
  list_df_bud[[d]] <- natal_turnover_ind_bud
}
df_bud <- rbindlist(list_df_bud)

#####

##
##
##### model including bud events #####
##
##
natal_turnover_ind_model_bud <- glmer(dom_in_natal ~
                                    age_cat:SEX +
                                    age_cat +
                                    SEX +
                                    (1|GROUP) +
                                    (1|Season),
                                  data = df_bud, 
                                  na.action = "na.fail",
                                  family = "binomial")
summary(natal_turnover_ind_model_bud)
drop1(natal_turnover_ind_model_bud, test = 'Chisq')

## model withouth interaction
natal_turnover_ind_model_bud_simple <- update(natal_turnover_ind_model_bud, .~.-age_cat:SEX)
drop1(natal_turnover_ind_model_bud_simple, test = 'Chisq') # model without interaction
summary(natal_turnover_ind_model_bud_simple)

# check model fit
simul_res <- DHARMa::simulateResiduals(fittedModel = natal_turnover_ind_model_bud, n = 1000)
testResiduals(simul_res) # they look fine

#####


##
##
##### table of results #####
##
##

## base table
table_natal_turn_over00 <- natal_turnover_ind_model_bud_simple %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `age_cat` = 'Subordinate age',
                   `SEX` = "Sex"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_natal_turn_over <- table_natal_turn_over00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=natal_turnover_ind_model_bud_simple) %>% 
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
table_name <- paste0("./tables/TABLE S7 - budding dominance acquisition natal.html")
gtsave(table_natal_turn_over, table_name)

#####
