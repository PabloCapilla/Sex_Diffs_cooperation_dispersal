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

source("../prep_repository/R_library/FUNCTION_drop1_output.R")

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
data$dom_in_anywhere <- ifelse(is.na(data$dominant_territory), 0, 1)
dom_in_natal_prebudding <- table(data$dom_in_natal, data$SEX)

##
## including dominance via budding
data_budding <- read.csv(file = "./data/07_dominance_via_budding.csv") # data
data_budding <- data_budding 

## changing budded males to dominance acq in turnover table

# males that budded to become dominants
budded_males <- {data_budding %>% 
    filter(budded_dom_male == 1) %>% 
    select(X1st_dom_male)}$X1st_dom_male

# females that budded to become dominants
budded_females <- {data_budding %>% 
    filter(budded_dom_female == 1) %>% 
    select(X1st_dom_female)}$X1st_dom_female

#new data including budded groups as 'inheritance'
all_inh_ind <- data
all_inh_ind$dom_in_natal_budding <- all_inh_ind$dom_in_natal

all_inh_ind$dom_in_natal_budding[all_inh_ind$BIRD.ID %in% budded_males] <- 1
all_inh_ind$dom_in_natal_budding[all_inh_ind$BIRD.ID %in% budded_females] <- 1

dom_in_natal_postbudding <- table(all_inh_ind$dom_in_natal_budding, all_inh_ind$SEX)

#####

##
##
##### contingency table for all dominance acquisition in natal group #####
##
##
all_inh_ind$dom_in_anywhere <- ifelse(all_inh_ind$dom_in_natal == 1, 1, all_inh_ind$dom_in_anywhere)

df_table_dom <- all_inh_ind %>% 
  filter(dom_in_anywhere == 1)
table(df_table_dom$dom_in_natal, df_table_dom$SEX)

prop.test(table(df_table_dom$dom_in_natal, df_table_dom$SEX),
          alternative = "two.sided",
          correct = T)





#####

##
##
##### Tenure length #####
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

table(data_ten$SEX)

# means
data_ten %>% 
  group_by(SEX) %>% 
  summarise(mean_tenure = mean(dom_tenure/365),
            se_tenure = sd(dom_tenure/365)/sqrt(n()),
            min_tenure = min(dom_tenure/365),
            max_tenure = max(dom_tenure/365))


t_test <- lm(dom_tenure ~ 
     SEX ,
   data = data_ten)
summary(t_test)

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
##### Table of results dominance tenure #####
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
#gtsave(tenure_model, "./tables/TABLE S14 - Dominance tenure length.html")

#####

##
##
##### What do male and female subordinate become when they leave home? do they take dominant or sub positions? #####
##
##
df_move <- all_inh_ind %>% 
  filter(Fate == 'EM') %>%  # birds that emigrated to known locations
  mutate(Gone.Date = dmy(Gone.Date),
         time_to_dom = as.numeric(ymd(dominant_start) - Gone.Date))
head(df_move)

##
## emigrate to become subordinate?
df_move$natal_move_to_sub <- ifelse(df_move$EM_group != df_move$dominant_territory, 1, 0)
df_move$natal_move_to_sub[is.na(df_move$natal_move_to_sub)] <- 1
initial_numbers <- table(df_move$natal_move_to_sub, df_move$SEX)

# refine for those that became dominant in the first group they first moved to but after a while
for(i in 1:nrow(df_move)) {
  moved_to_sub <- df_move$natal_move_to_sub[i]
  
  if(moved_to_sub == 1){
    next()
  } else {
    df_move$natal_move_to_sub[i] <- ifelse(df_move$time_to_dom[i] > 14, 1, 0)
  }
}
final_numbers <- table(df_move$natal_move_to_sub, df_move$SEX)

## analysis of probability to emigrate to subordinate position
model_em_sub <- glmer(natal_move_to_sub ~ 
                        SEX +
                        (1|Season) +
                        (1|GROUP),
                      data = df_move, 
                      family = 'binomial')
summary(model_em_sub)
drop1(model_em_sub, test = 'Chisq')

## base table
table_em_sub00 <- model_em_sub %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   SEX = "Sex"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_em_sub <- table_em_sub00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model_em_sub) %>% 
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
gtsave(table_em_sub, "./tables/TABLE S10 - Probability of emigrating to subordinate position.html")

