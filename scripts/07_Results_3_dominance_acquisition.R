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

## full dispersal events
data <- readRDS("./data/06_dispersal_full_table_analysis_data.RDS")
head(data)

## budding data
data_budding <- read.csv(file = "./data/budding_data.csv") # data
data_budding <- data_budding %>% filter(DISCARDED == "no") # removing one female that budded being dominant already (from Dassie to Canary)



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

# data
natal_turnover_ind <- data %>%
  filter(SEX != "U") %>% 
  mutate(SEX = as.character(SEX)) %>% 
  mutate(dom_age = ymd(dominant_start) - dmy(BBBirthday)) %>% 
  filter(age_fate > (365/2)) # birds that are at least 6m old
head(natal_turnover_ind)

table(natal_turnover_ind$dom_in_natal)

#####

##
##
##### 1 -  Natal turn over #####
##
##

# model
natal_turnover_ind_model <- glmer(dom_in_natal ~ 
                                    SEX +
                                    (1|GROUP) +
                                    (1|Season),
                                  data = natal_turnover_ind, 
                                  na.action = "na.fail",
                                  family = "binomial")
summary(natal_turnover_ind_model)

#####

##
##
##### Table of results S8a #####
##
##

## base table
table_natal_turn_over00 <- natal_turnover_ind_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
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
                        y = drop1_output(x=natal_turnover_ind_model) %>% 
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
gtsave(table_natal_turn_over, "./tables/TABLE S8a - Dominance acquisition natal.html")

#####


##
##
##### Plotting model predictions #####
##
##

# summary estimates
data_summary_plot <- natal_turnover_ind %>% 
  group_by(SEX) %>% 
  summarise(mean_prob = mean(dom_in_natal),
            se_prob = sd(dom_in_natal)/sqrt(n()))

# model predictions
df_predict <- expand.grid(SEX = c("M", "F"))
df_predict$fit <- predict(object = natal_turnover_ind_model, 
                          df_predict, re.form = NA, type = "link")

mm <- model.matrix(~ SEX,
                   data = df_predict)
pvar1 <- diag(mm %*% tcrossprod(vcov(natal_turnover_ind_model),mm))
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
natal_turnover_ind$sex_name <- ifelse(natal_turnover_ind$SEX == "F", "female", "male")




dom_acq_natal <- ggplot(natal_turnover_ind, aes(x = SEX, 
                                  y = dom_in_natal, 
                                  fill = SEX, color = SEX)) +
  geom_point(position = position_jitter(width = 0.03, 
                                        height = 0.02),
             size = 1.5,
             alpha = 0.15) +
  geom_errorbar(data = df_predict, 
                aes(ymin = plow_resp, 
                    ymax = phi_resp,
                    y = fit_resp),
                color = "black",
                width = 0) +
  geom_point(data = df_predict,
             aes(y = fit_resp),
             size = 4,
             color = "black",
             shape = 21) +
  theme_bw() +
  labs(x = "Individual sex", 
       y = "Prob. of dominance acquisition in natal territory") + 
  theme(axis.title.x = element_text(family = "Arial", size = 10),
        axis.title.y = element_text(family = "Arial", size = 10),
        axis.text.x = element_text(family = "Arial", size = 10),
        axis.text.y = element_text(family = "Arial", size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        legend.text = element_text(family = "Arial", size = 10),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#9970ab", "#5aae61")) +
  scale_color_manual(values = c("#9970ab", "#5aae61")) +
  scale_x_discrete(labels = c("female", "male")) 


ggsave(filename = "./1_Age_1_disp/plots/prob_dom_acq.jpeg", 
       plot = dom_acq_natal, 
       units = "mm",
       device = "jpeg", 
       width = 89,
       height = 89)

#####


##
##### 2 -  Natal turn over + budding #####
##

##
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


##
## model for inheritance of dom in natal group or territory
all_inh_ind_model <- glmer(dom_in_natal ~ 
                     SEX +
                     (1|GROUP) +
                     (1|Season),
                   data = all_inh_ind, 
                   na.action = "na.fail",
                   family = "binomial")
summary(all_inh_ind_model)

#####

##
##
##### Table of results S8b #####
##
##

## base table
table_budding_turn_over00 <- all_inh_ind_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `SEX` = "Sex"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_budding_turn_over <- table_budding_turn_over00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=natal_turnover_ind_model) %>% 
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
gtsave(table_budding_turn_over, "./tables/TABLE S8b - Dominance acquisition natal and budding.html")

#####


##
##
##### Dominance turnover in natal group + including budding: Fisher's exact test #####
##
##

##
##  natal turn over

# data for test
dom_birds_only <- data %>%
  filter(SEX != "U") %>% 
  mutate(SEX = as.character(SEX)) %>% 
  filter(dominant_bird == 1) %>%
  mutate(dom_age = ymd(dominant_start) - dmy(BBBirthday)) %>% 
  filter(age_fate > (365/2)) # birds that are at least 0.5y old


# natal turn overs table and F test
table(dom_birds_only$dom_in_natal, dom_birds_only$SEX)
fisher.test(dom_birds_only$dom_in_natal, dom_birds_only$SEX)

##
## including budded groups as 'inheritance' of natal territory


table(data_budding$budded_dom_male)
table(data_budding$budded_dom_female)

# natal turn over table
dom_turnover_natal <- rbind(c(23, 25),
                            c(9, 7))
rownames(dom_turnover_natal) <- c("0", "1")
colnames(dom_turnover_natal) <- c("F", "M")

# budding table
budding_table <- rbind(c(2, 3),
                       c(9, 6))
rownames(budding_table) <- c("0", "1")
colnames(budding_table) <- c("F", "M")

fisher.test(budding_table)

# combined table and test
dom_aq_total <- dom_turnover_natal + budding_table
fisher.test(dom_aq_total)

#####








