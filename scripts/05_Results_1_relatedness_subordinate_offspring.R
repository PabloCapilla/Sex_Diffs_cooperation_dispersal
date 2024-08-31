###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2024-08-31
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
data <- readRDS(file = "data/05_relatedness_analysis_data.RDS")
head(data)

## sample sizes
nrow(data) # total obs
table(data$sub_sex)

# sample size per sex
sex_table <- data %>% 
  group_by(sub_id) %>% 
  filter(row_number() == 1)

nrow(sex_table)
table(sex_table$sub_sex)

##
##
###### Relatedness to the offspring #####
##
##
model_relatedness <- lmer(rel_to_offspring ~ sub_sex + (1|clutch_id), 
                        data = data,
                        na.action = "na.fail")
summary(model_relatedness)
drop1(model_relatedness, test = 'Chisq')


## base table
table_relatedness00 <- model_relatedness %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   sub_sex = "Subordinate sex"),
                 estimate_fun = ~ style_number(.x, digits = 3)) 

## add features
table_relatedness <- table_relatedness00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=model_relatedness, family = "binomial") %>% 
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
gtsave(table_relatedness, "./tables/TABLE S3 - Relatedness Subordinate to offspring.html")

#####

##
##
##### plotting age trajectories #####
##
##
df_predict <- expand.grid(
  sub_sex = c("F", "M")
)

df_predict$fit <- predict(model_relatedness, 
                              df_predict, 
                              re.form = NA, 
                              type = "link")
mm <- model.matrix(~
                     sub_sex,
                   data = df_predict)

pvar1 <- diag(mm %*% tcrossprod(vcov(model_relatedness),mm))
cmult <- 1 ## could use 1.96
df_predict <- data.frame(
  df_predict
  , plo = df_predict$fit-cmult*sqrt(pvar1)
  , phi = df_predict$fit+cmult*sqrt(pvar1)
)


df_predict$sex_name <- ifelse(df_predict$sub_sex == "F", "female", "male")
data$sex_name <- ifelse(data$sub_sex == "F", "female", "male")


##
## individual prob of helping
plot_relatedness <- ggplot(df_predict,
                         aes(y = fit, 
                             x = sex_name,
                             color = sex_name, 
                             fill = sex_name)) +
  geom_point(data = data,
             aes(y = rel_to_offspring,
                 x = sex_name,
                 color = sex_name),
             position = position_jitter(width = 0.2),
             size = 2.5,
             alpha = 0.10) +
  geom_errorbar(aes(ymin = plo, ymax = phi), 
                width = 0,
                color = "black",
                position = position_dodge(width = 0.25)) +
  geom_point(size = 5,
             shape = 21,
             color = "black",
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  labs(x = " ", 
       y = expression(atop("Relatedness", "subordinate to offspring"))) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Arial", color = "black", size = 15),
        axis.text.x = element_text(family = "Arial", color = "black", size = 12),
        axis.text.y = element_text(family = "Arial", color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        legend.text = element_text(family = "Arial", size = 10),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#9970ab", "#5aae61")) +
  scale_color_manual(values = c("#9970ab", "#5aae61")) 

saveRDS(file = "./plots/Figure 1e.RDS", 
        object = plot_relatedness)

ggsave(filename = "./plots/Figure 1e.png", 
       plot = plot_relatedness, 
       units = "mm",
       device = "png", 
       width = 85,
       height = 115)




