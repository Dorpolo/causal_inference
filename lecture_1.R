# loading packages
library(dplyr)
library(ggplot2)
library(reshape2)

data <- read.csv('data/rhc.csv', header = T, sep = ';')
N = data %>% nrow()

# Calculate Naive Estimator
df <- data %>% 
  group_by(treatment) %>% 
  summarise(avg_dead = mean(death),
            sum_dead = sum(death),
            var_dead = var(death)) %>%
  mutate(ATE = avg_dead[treatment == 1] - avg_dead[treatment == 0])

naive_ATE = df$ATE[1]
print("Naive ATE estimator:"); print(naive_ATE)

# Confidence interval using asymptotic variance
trt_data <- df %>% filter(treatment == 1)
n_trt <- trt_data$sum_dead
var_trt <- trt_data$var_dead

ctrl_data <- df %>% filter(treatment == 0)
n_ctrl <- ctrl_data$sum_dead
var_ctrl <- ctrl_data$var_dead

alpha = 0.05
Z = qnorm(1-(alpha/2))
naive_ATE_sd <- sqrt(var_trt/n_trt) + sqrt(var_ctrl/n_ctrl)
CI_asymp <- naive_ATE + c(-1, 1)*Z*naive_ATE_sd
print("Asymptotic CI for ATE is:"); print(CI_asymp)

# Confidence interval using bootstrap
B <- 1000
naive_vec <- rep(NA, B)

boot_it <- function(N){
  index = sample(x = 1:N, size = N, replace = T)
  tuncated_df <- data[index,] %>% 
    group_by(treatment) %>% 
    summarise(avg_dead = mean(death), .groups = 'drop') %>%
    mutate(ATE = avg_dead[treatment == 1] - avg_dead[treatment == 0])
  return((tuncated_df %>% head(1))$ATE)
}

boot_res = sapply(rep(N, B), boot_it)

# option 1
boot_sd = sd(boot_res)
CI_boot_option_1 = naive_ATE + c(1, -1)*Z*boot_sd
print("Bootstraped CI for ATE (option 1) is:"); print(CI_boot_option_1)

# option 2
lower = quantile(x = boot_res, probs = alpha/2)
upper = quantile(x = boot_res, probs = 1-(alpha/2))
print("Bootstraped CI for ATE (option 1) is:"); print(c(lower, upper))



# Calculate Naive Estimator
df_stand <- data %>% 
  group_by(treatment, sepsis) %>% 
  summarise(avg_dead = mean(death),
            sum_dead = sum(death),
            var_dead = var(death),
            n = n(),
            adj_var = var_dead/n) 

prop_sepsis <- mean(data$sepsis)

stand_m_trt <- prop_sepsis*df_stand[4, 3] + (1 - prop_sepsis)*df_stand[3, 3]
stand_m_ctrl <- prop_sepsis*df_stand[2, 3] + (1 - prop_sepsis)*df_stand[1, 3]
ATE.stand <- (stand_m_trt - stand_m_ctrl)$avg_dead

print("ATE standartized estimator:"); print(ATE.stand)

# CI - ATE Estimator
n00 = df_stand$n[1]
n01 = df_stand$n[2]
n10 = df_stand$n[3]
n11 = df_stand$n[4]
m00 = df_stand$avg_dead[1]
m01 = df_stand$avg_dead[2]
m10 = df_stand$avg_dead[3]
m11 = df_stand$avg_dead[4]
p1 = prop_sepsis
v.m00 = df_stand$adj_var[1]
v.m01 = df_stand$adj_var[2]
v.m10 = df_stand$adj_var[3]
v.m11 = df_stand$adj_var[4]

ATE_stand_asymp_var = (((m11-m10+m00-m01)*p1*(1-p1)/N) + (((1-p1)^2)*(v.m10+v.m00)) + ((p1^2)*(v.m11+v.m01)))

ATE_stand_asymp_sd = sqrt(ATE_stand_asymp_var)
CI_stand = ATE.stand + c(-1, 1)*Z*ATE_stand_asymp_sd
print("CI in 95% for ATE standartized estimator:"); print(CI_stand)

# Confidence interval using bootstrap
B <- 1000
stand_vec <- rep(NA, B)

boot_it_stand <- function(N){
  index = sample(x = 1:N, size = N, replace = T)
  tuncated_df <- data[index,] %>% 
    group_by(treatment, sepsis) %>% 
    summarise(avg_dead = mean(death), .groups = 'drop') 
  
  prop_sepsis_boot = mean(data[index,]$sepsis)
  stand_m_trt_boot = prop_sepsis_boot*tuncated_df[4, 3] + (1-prop_sepsis_boot)*tuncated_df[3, 3]
  stand_m_ctrl_boot = prop_sepsis_boot*tuncated_df[2, 3] + (1-prop_sepsis_boot)*tuncated_df[1, 3]
  output = stand_m_trt_boot-stand_m_ctrl_boot
  return(output$avg_dead)
}

stand_boot_res = sapply(rep(N, B), boot_it_stand)

# CI - bootstrap method - option 1
sd.boot.ATE.stand = sd(stand_boot_res)
CI_boot_a = ATE.stand + c(-1, 1)*Z*sd.boot.ATE.stand
print("bootstraped CI in 95% for ATE standartized estimator, option 1:"); print(CI_boot_a)


# option 2
lower = quantile(x = stand_boot_res, probs = alpha/2)
upper = quantile(x = stand_boot_res, probs = 1-(alpha/2))
print("bootstraped CI in 95% for ATE standartized estimator, option 2:"); print(c(lower, upper))










alpha = 0.05


df_stand <- data %>% 
  group_by(treatment, sepsis) %>% 
  summarise(avg_dead = mean(death),
            sum_dead = sum(death),
            var_dead = var(death),
            n = n(),
            adj_var = var_dead/n) 

prop_sepsis <- mean(data$sepsis)

stand_m_trt <- prop_sepsis*df_stand[4, 3] + (1 - prop_sepsis)*df_stand[3, 3]
stand_m_ctrl <- prop_sepsis*df_stand[2, 3] + (1 - prop_sepsis)*df_stand[1, 3]
ATE.stand <- (stand_m_trt - stand_m_ctrl)$avg_dead

print("ATE standartized estimator:"); print(ATE.stand)

# CI - ATE Estimator
n00 = df_stand$n[1]
n01 = df_stand$n[2]
n10 = df_stand$n[3]
n11 = df_stand$n[4]
m00 = df_stand$avg_dead[1]
m01 = df_stand$avg_dead[2]
m10 = df_stand$avg_dead[3]
m11 = df_stand$avg_dead[4]
p1 = prop_sepsis
v.m00 = df_stand$adj_var[1]
v.m01 = df_stand$adj_var[2]
v.m10 = df_stand$adj_var[3]
v.m11 = df_stand$adj_var[4]

ATE_stand_asymp_var = (((m11-m10+m00-m01)*p1*(1-p1)/N) + (((1-p1)^2)*(v.m10+v.m00)) + ((p1^2)*(v.m11+v.m01)))

ATE_stand_asymp_sd = sqrt(ATE_stand_asymp_var)
CI_stand = ATE.stand + c(-1, 1)*Z*ATE_stand_asymp_sd
print("CI in 95% for ATE standartized estimator:"); print(CI_stand)

# Confidence interval using bootstrap
B <- 1000
stand_vec <- rep(NA, B)

boot_it_stand <- function(N){
  index = sample(x = 1:N, size = N, replace = T)
  tuncated_df <- data[index,] %>% 
    group_by(treatment, sepsis) %>% 
    summarise(avg_dead = mean(death), .groups = 'drop') 
  
  prop_sepsis_boot = mean(data[index,]$sepsis)
  stand_m_trt_boot = prop_sepsis_boot*tuncated_df[4, 3] + (1-prop_sepsis_boot)*tuncated_df[3, 3]
  stand_m_ctrl_boot = prop_sepsis_boot*tuncated_df[2, 3] + (1-prop_sepsis_boot)*tuncated_df[1, 3]
  output = stand_m_trt_boot-stand_m_ctrl_boot
  return(output$avg_dead)
}

stand_boot_res = sapply(rep(N, B), boot_it_stand)

# CI - bootstrap method - option 1
sd.boot.ATE.stand = sd(stand_boot_res)
CI_boot_a = ATE.stand + c(-1, 1)*Z*sd.boot.ATE.stand
print("bootstraped CI in 95% for ATE standartized estimator, option 1:"); print(CI_boot_a)


# option 2
lower = quantile(x = stand_boot_res, probs = alpha/2)
upper = quantile(x = stand_boot_res, probs = 1-(alpha/2))
print("bootstraped CI in 95% for ATE standartized estimator, option 2:"); print(c(lower, upper))

library(dplyr)
library(ggplot2)
library(ggcharts)

df_hist = tibble(res = c(stand_boot_res, boot_res),
                 type = c(rep('Standardized', length(stand_boot_res)),
                          rep('Naive', length(boot_res))))

df_mu = tibble(res = c(ATE.stand, naive_ATE),
               type = c('Standardized', 'Naive'))
   
df_hist %>%
ggplot(aes(x=res, fill = type)) + 
  geom_density(alpha = 0.5) +
  theme_hermit() + 
  labs(title = 'ATE Estimator Distribution',
       subtitle = 'by Calculation Method',
       fill = '',
       color = '') + 
  theme(legend.position = 'top') + 
  geom_vline(data=df_mu, aes(xintercept=res, col=type)) + 
  scale_fill_manual(values = c('#F6D55C', '#ED553B', '#3CAEA3', '#20639B')) +
  scale_color_manual(values = c('#F6D55C', '#ED553B', '#3CAEA3', '#20639B')) 
  

# loading packages
library(dplyr)
library(ggplot2)
library(reshape2)

data <- read.csv('data/rhc.csv', header = T, sep = ';') %>% 
  mutate(age = as.numeric(gsub(',', '', age)))

N = data %>% nrow()

# treatment + female + ARF + age + sepsis

logistic_model = glm(
  formula = death ~ treatment + female + ARF + age + sepsis,
  data = data, family = 'binomial'
  )

all_trt = all_untrt = data
all_trt$treatment <- 1
all_untrt$treatment <- 0

pred.trt = predict(object = logistic_model, newdata = all_trt, type = 'response')
pred.untrt = predict(object = logistic_model, newdata = all_untrt, type = 'response')

# ATE
ate.logistic_model_stand_diff <- mean(pred.trt) - mean(pred.untrt)
ate.logistic_model_stand_diff

# Ratio
ate.logistic_model_stand_ratio <- mean(pred.trt)/mean(pred.untrt)
ate.logistic_model_stand_ratio

# Odds Ratio
ate.logistic_model_stand_OR <- (mean(pred.trt)/(1-mean(pred.trt)))/(mean(pred.untrt)/(1-mean(pred.untrt)))
ate.logistic_model_stand_OR

data %>% head(10)

data 
