library(dplyr)
library(ggplot2)
library(reshape2)
library(stats)
library(reticulate)
install.packages('rmdformats')
library(rmdformats)
data <- read.csv('~/Desktop/private/TAU/tau_stats/data/rhc.csv', header = T, sep = ';')

library(optmatch)
library(MatchIt)


data <- read.csv('~/Desktop/private/TAU/tau_stats/data/rhc.csv', header = T, sep = ';')
rhc_data <- data %>% 
  mutate(death = as.integer(gsub(',', '', death)),
         age = as.numeric(gsub(',', '', age))/(10^7),
         X = NULL,
         gender = factor(female, levels = c(0,1), labels = c('Male', 'Female')),
         labeled_trt = factor(treatment, levels = c(0,1), labels = c('Untreated', 'Treated')),
         meanbp1 = as.numeric(meanbp1),
         id = row_number())

10^7 == 10000000
n = 10000000
rhc_data %>%
  ggplot(aes(x=age)) + geom_histogram()

class(rhc_data$meanbp1)

set.seed(10)

calc_it <- function(x){
  # b1
  x = rnorm(n = 500, mean = 0, sd = 1)
  v = rnorm(n = 500, mean = 0, sd = 1)
  beta = list(`0` = 100, `1` = 2, `2` = 3.5, `3` = 5)
  sigma = 3
  eps = rnorm(n = 500, mean = 0, sd = sigma)
  
  y0 = beta$`0` + beta$`1`*0 + beta$`2`*x + beta$`3`*v + eps
  y1 = beta$`0` + beta$`1`*1 + beta$`2`*x + beta$`3`*v + eps
  
  # b2
  gamma = list(`0` = log(0.25), `1` = log(2))
  p = (exp(gamma$`0` + gamma$`1`*x)) / (1+exp(gamma$`0` + gamma$`1`*x))
  a = rbinom(n = 500, size = 1, prob = p)
  
  # b3
  y = rep(NA, 500)
  y[which(a==0)] = y0[which(a==0)]
  y[which(a==1)] = y1[which(a==1)]
  
  # c1
  ATE_naive = mean(y[which(a==1)]) - mean(y[which(a==0)])
  
  # c2
  df_sim = as_tibble(cbind(y, a, x, v))
  matching_opt_ps <- MatchIt::matchit(formula = a ~ x + v, method = 'optimal', distance = 'logit', data = df_sim)
  matched_data_opt_ps = match.data(matching_opt_ps)
  
  # c3
  ATE_matching_opt_ps = (matched_data_opt_ps %>% 
                           summarise(diff = mean(y[a==1]) - mean(y[a==0])))$diff
  
  # c4
  lm_ax = lm(formula = y ~ a + x, data = matched_data_opt_ps)
  ATE_lm_ax = coef(lm_ax)[2]
  
  # c5
  lm_axv = lm(formula = y ~ a + x + v, data = matched_data_opt_ps)
  ATE_lm_axv = coef(lm_axv)[2]
  
  return(c(ATE_naive, ATE_matching_opt_ps, ATE_lm_ax, ATE_lm_axv))
}

res = sapply(rep(0, 10), calc_it)
res_data = as_tibble(t(res))
names(res_data) = c('Naive', 'Matching Optimal, PS', 'LM: Y ~ A + X', 'LM: Y ~ A + X + V')
output = res_data %>% melt() %>% group_by(variable) %>% summarise(Mean = mean(value), SD = sd(value))


output



exact_match_11 <- MatchIt::matchit(formula = treatment ~ female + age + meanbp1,
                                   method = 'exact',
                                   data = rhc_data)
df_exact_match_11 <- MatchIt::match.data(exact_match_11)

present_first_3_rows <- function(method = 'glm', distance = 'exact'){
  match_instance <- MatchIt::matchit(
    formula = treatment ~ female + age + meanbp1,
    method = method,
    distance = distance,
    data = rhc_data)
  
  df <- MatchIt::match.data(exact_match_11) %>% 
    filter(id <= 3) %>% select(treatment, female, age, meanbp1)
  
  return(formattable(df,align = rep('c', ncol(df))))
}

set.seed(10)
match_instance <- MatchIt::matchit(
  formula = treatment ~ female + age + meanbp1,
  method = 'optimal',
  distance = 'mahalanobis',
  data = rhc_data)
  
df_matches <- MatchIt::match.data(match_instance) %>% 
  select(id, subclass, labeled_trt, gender, age, meanbp1) %>%
  group_by(subclass) %>% filter(n() > 1) %>% 
  mutate(rn = row_number(id))

df_matches %>% filter(rn == 1) %>%
  inner_join(df_matches %>% filter(rn == 2), by='subclass') %>%
  select(-c('rn.x', 'rn.y'))

match_instance <- MatchIt::matchit(
  formula = treatment ~ female + age + meanbp1,
  method = 'nearest',
  distance = 'mahalanobis',
  data = rhc_data %>% head(1000))

df_matches <- MatchIt::match.data(match_instance) %>% 
  select(id, subclass, labeled_trt, gender, age, meanbp1) %>%
  group_by(subclass) %>% filter(n() > 1) %>% 
  mutate(rn = row_number(id))

df_matches %>% filter(rn == 1) %>%
  inner_join(df_matches %>% filter(rn == 2), by='subclass') %>%
  select(-c('rn.x', 'rn.y'))

















set.seed(10)
match_instance <- MatchIt::matchit(
  formula = treatment ~ female + age + meanbp1,
  method = 'optimal',
  distance = 'logit',
  data = rhc_data %>% head(1000))

df_matches <- MatchIt::match.data(match_instance) %>% 
  select(id, subclass, labeled_trt, gender, age, meanbp1) %>%
  group_by(subclass) %>% filter(n() > 1) %>% 
  mutate(rn = row_number(id))

output = df_matches %>% filter(rn == 1) %>%
  inner_join(df_matches %>% filter(rn == 2), by='subclass') %>%
  select(-c('rn.x', 'rn.y')) %>% arrange(id.x)

formattable(output %>% filter(id <= 3), align = rep('c', ncol(output)))


matched_index <- as.integer(rownames(match_instance$match.matrix))
original_index <- as.integer(match_instance$match.matrix[,1])

df_map <- tibble::tibble(id=original_index, matched_id=matched_index) %>%
  arrange(id)

model_df <- rhc_data %>% select(id, labeled_trt, gender, age, meanbp1)

output <- model_df %>% inner_join(df_map, by='id') %>%
  inner_join(model_df, by = c('matched_id'='id'))

formattable(output %>% filter(id <= 3), align = rep('c', ncol(output)))









MatchIt::match.data(match_instance) %>% filter(subclass == 124)

matched_index <- as.integer(rownames(match_instance$match.matrix))
original_index <- as.integer(match_instance$match.matrix[,1])

df_map <- tibble::tibble(id=original_index, matched_id=matched_index) %>%
  arrange(id)

model_df <- rhc_data %>% select(id, labeled_trt, gender, age, meanbp1)

output <- model_df %>% inner_join(df_map, by='id') %>%
  inner_join(model_df, by = c('matched_id'='id'))

formattable(output %>% filter(id <= 3), align = rep('c', ncol(output)))





set.seed(10)
match_instance <- MatchIt::matchit(
  formula = treatment ~ female + age + meanbp1,
  method = 'nearest',
  distance = 'logit',
  data = rhc_data)

df_matches <- MatchIt::match.data(match_instance) %>% 
  select(id, subclass, labeled_trt, gender, age, meanbp1) %>%
  group_by(subclass) %>% filter(n() > 1) %>% 
  mutate(rn = row_number(id))

output = df_matches %>% filter(rn == 1) %>%
  inner_join(df_matches %>% filter(rn == 2), by='subclass') %>%
  select(-c('rn.x', 'rn.y'))

formattable(output %>% filter(id.x <= 3), align = rep('c', ncol(output)))





find_match <- function(method = 'exact', distance = NULL){
  match_instance <- MatchIt::matchit(
    formula = treatment ~ female + age + meanbp1,
    method = method,
    distance = distance,
    data = rhc_data)
  
  df_matches <- MatchIt::match.data(match_instance) %>% 
    select(id, subclass, labeled_trt, gender, age, meanbp1) %>%
    group_by(subclass) %>% filter(n() > 1) %>% 
    mutate(rn = row_number(id))
  
  output <- df_matches %>% filter(rn == 1) %>%
    inner_join(df_matches %>% filter(rn == 2), by='subclass') %>%
    select(-c('rn.x', 'rn.y'))
  table = formattable(output %>% filter(id.x <= 3), align = rep('c', ncol(output)))
  n_match = nrow(output)
  return(list(t=table, n=n_match))
}
find_match(method = 'nearest', distance='mahalanobis')$t









data <- read.csv('~/Desktop/private/TAU/tau_stats/data/rhc.csv', header = T, sep = ';')
rhc_data <- data %>% 
  mutate(death = as.integer(gsub(',', '', death)),
         meanbp1 = as.numeric(ifelse(is.na(meanbp1), NULL, meanbp1)),
         age = as.numeric(gsub(',', '', age))/(10^7),
         X = NULL,
         gender = factor(female, levels = c(0,1), labels = c('Male', 'Female')),
         labeled_trt = factor(treatment, levels = c(0,1), labels = c('Untreated', 'Treated')),
         id = row_number())

data %>% filter(is.na(meanbp1))
str(data)


if(length(new.packages)) install.packages(new.packages)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(reshape2)
library(formattable)
library(ggcharts)
library(optmatch)
library(MatchIt)

data <- read.csv('~/Desktop/private/TAU/tau_stats/data/rhc.csv', header = T, sep = ';')
rhc_data <- data %>% 
  mutate(death = as.integer(gsub(',', '', death)),
         meanbp1 = as.numeric(gsub(',', '', meanbp1)),
         age = as.numeric(gsub(',', '', age))/(10^7),
         X = NULL,
         gender = factor(female, levels = c(0,1), labels = c('Male', 'Female')),
         labeled_trt = factor(treatment, levels = c(0,1), labels = c('Untreated', 'Treated')),
         id = row_number())

find_match <- function(method = 'exact', distance = NULL){
  match_instance <- MatchIt::matchit(
    formula = treatment ~ female + age + meanbp1,
    method = method,
    distance = distance,
    data = rhc_data)
  
  raw_data <- MatchIt::match.data(match_instance) 
  
  df_matches <- raw_data %>% select(id, subclass, labeled_trt, gender, age, meanbp1) %>% 
    group_by(subclass) %>% filter(n() > 1) %>% mutate(rn = row_number(id))
  
  output <- df_matches %>% filter(rn == 1) %>%
    inner_join(df_matches %>% filter(rn == 2), by='subclass') %>%
    select(-c('rn.x', 'rn.y'))
  
  table <- formattable(output %>% filter(id.x <= 3), align = rep('c', ncol(output)))
  return(list(t=table, n=nrow(output), df=raw_data))
}

output2 = find_match(method = 'nearest', distance='mahalanobis')
output2$t

M = 2184
calcualtor <- function(m){
  output2$df %>% filter(subclass == m) %>% 
    summarise(ATT = mean(treatment == 1))
}

& (output2$df$treatment==1)

output2$df[(output2$df$subclass==51) & (output2$df$treatment==1), 'death']

