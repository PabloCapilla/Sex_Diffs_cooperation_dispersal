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
data <- readRDS("./data/06_dispersal_full_table_analysis_data.RDS")
head(data)
table(data$Start.Context)

# dominant in natal terr?
data$dom_in_natal <- ifelse(is.na(data$dominant_territory), 0, 
                            ifelse(data$dominant_bird == 1 &
                                     data$dominant_territory == data$Start.Group, 1, 0))
dom_in_natal_prebudding <- table(data$dom_in_natal, data$SEX)


##
## including dominance via budding
data_budding <- read.csv(file = "./data/budding_data.csv") # data
data_budding <- data_budding %>% filter(DISCARDED == "no") # removing one female that budded being dominant already (from Dassie to Canary)

## changing budded males to dominance acq in turnover table

# males that budded to become dominants
budded_males <- {data_budding %>% 
    filter(Start.Context_dom_male == "HATCHED") %>% 
    filter(budded_dom_male == 1) %>% 
    select(X1st_dom_male)}$X1st_dom_male

# females that budded to become dominants
budded_females <- {data_budding %>% 
    filter(Start.Context_dom_female == "HATCHED") %>% 
    filter(budded_dom_female == 1) %>% 
    select(X1st_dom_female)}$X1st_dom_female

 new data including budded groups as 'inheritance'
all_inh_ind <- data
all_inh_ind$dom_in_natal[all_inh_ind$BIRD.ID %in% budded_males] <- 1
all_inh_ind$dom_in_natal[all_inh_ind$BIRD.ID %in% budded_females] <- 1

dom_in_natal_postbudding <- table(all_inh_ind$dom_in_natal, all_inh_ind$SEX)

#####

##
##
##### TENURE LENGTH #####
##
##

# stop tenures on 2016-06-01 the end of the census period
all_inh_ind$dominant_stop <- ifelse(all_inh_ind$dominant_stop == "2099-01-01", 
                                    "2016-06-01", 
                                    all_inh_ind$dominant_stop)

# tenure length
all_inh_ind$dom_tenure <- as.numeric(ymd(all_inh_ind$dominant_stop) - ymd(all_inh_ind$dominant_start))
hist(all_inh_ind$dom_tenure)

# data for model
data_ten <- all_inh_ind %>% 
  filter(dominant_bird == 1) %>% 
  mutate(dom_age = ymd(dominant_start) - dmy(BBBirthday)) %>% 
  filter(age_fate > (365/2))

# means
data_ten %>% 
  group_by(SEX) %>% 
  summarise(mean_tenure = mean(dom_tenure/365),
            se_tenure = sd(dom_tenure/365)/sqrt(n()))


# model
ten_model_full <- lmer(dom_tenure ~ 
                         SEX +
                         (1|GROUP) +
                         (1|Season),
                       data = data_ten,
                       na.action = "na.fail",
                       REML = F)
summary(ten_model_full)

#####

##
##
##### Table of results S9 #####
##
##

## base table
tenure_model00 <- ten_model_full %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `SEX` = "Sex"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
tenure_model <- tenure_model00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=ten_model_full) %>% 
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
gtsave(tenure_model, "./tables/TABLE S9 - Dominance tenure length.html")

#####



##
##
##### INCIDENCE OF TURNOVERS #####
##
##
table(data_ten$SEX, data_ten$dominant_bird)


##
##
##### TURNOVERS BY NATALS #####
##
##
# means
data_ten %>% 
  group_by(SEX) %>% 
  summarise(mean_natal = mean(dom_in_natal),
            se_natal = sd(dom_in_natal)/sqrt(n()))

natalturnover_model_full <- glmer(dom_in_natal ~ 
                                    SEX +
                                    (1|GROUP) +
                                    (1|Season),
                                  data = data_ten,
                                  na.action = "na.fail",
                                  family = "binomial")
summary(natalturnover_model_full)

## base table
natalturnover_model00 <- natalturnover_model_full %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `SEX` = "Sex"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
natalturnover_model00 <- natalturnover_model00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=natalturnover_model_full) %>% 
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
gtsave(tenure_model, "./tables/TABLE S10 - Dominance turnover.html")

#####

# proportion of natal turnovers
table(data_ten$dom_in_natal)
nrow(data_ten[data_ten$dom_in_natal==1,])/nrow(data_ten)



