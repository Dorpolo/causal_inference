# loading packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggcharts)

### Device Type

df_device_type = read.csv(file = '~/Desktop/udger_data/device_type.csv') %>%
  mutate(device_type= ifelse(device_type == 'Tablet',
                             paste0(device_type, '-', library_type), device_type),
         device_type_adj = ifelse(X. > 0.01, device_type, ''),
         label_pct = paste0(device_type_adj, '\n(', X., '%)'))

dev_ordered = (df_device_type %>% arrange(X.))$device_type


df_device_type %>% 
  select(library=library_type, device_type, `%`= X., label_pct) %>%
  mutate(device_type = factor(device_type, levels = dev_ordered)) %>%
  ggplot(aes(x=library, y=`%`, fill=device_type)) + 
  geom_bar(stat='identity', alpha = 0.5) +
  theme_ggcharts() +
  geom_text(aes(label=label_pct), position = position_stack(0.5), vjust=0.5, size = 4) +
  labs(title = 'Device Type',
       subtitle = 'Udger vs. BitWalker',
       fill = '',
       color = '') + 
  theme(legend.position = 'right')

### Device Type - Diff

df_device_type_diff = read.csv(file = '~/Desktop/udger_data/device_type_diff.csv') %>%
  select(bitwalker = bitwalker_device_type,
         udger = udger_device_type, 
         `%` = X..within.device_type) %>%
         mutate(label_pct = ifelse(`%`>0.1, paste0(udger,' (%',`%`,')'), ''))

dev_dif_ordered = unique((df_device_type_diff %>% arrange(`%`))$udger)

df_device_type_diff %>% group_by(bitwalker) %>%
  mutate(udger = factor(udger, levels = dev_dif_ordered)) %>%
  ggplot(aes(x=bitwalker, y=`%`, fill=udger)) + 
  geom_bar(stat='identity', alpha = 0.5) +
  # theme_hermit() + 
  theme_ggcharts() +
  coord_flip() +
  geom_text(aes(label=label_pct), 
            position = position_stack(0.5), 
            vjust=0.5, 
            size = 3, col='black') +
  labs(title = 'Device Type',
       subtitle = 'BitWalker by Udger Values',
       fill = 'Udger') + 
  theme(legend.position = 'top')
