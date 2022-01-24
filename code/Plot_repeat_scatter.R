# Title     :Plot repeat
# Objective :
# Created by: yiquan
# Created on: 12/5/21
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
require(cowplot)
library(readxl)
library(reshape2)
library(writexl)
HA_df <- read_excel('result/HA_EggPassaging_all.xlsx') %>%
  select(Mutation,Group,Freq)%>%
  filter(Group!='NA')
#long to wide
HA_df_wide <- dcast(HA_df, Mutation~Group,value.var = "Freq")
HA_df_wide[is.na(HA_df_wide)] <- 0
print (paste("Pearson Cor:", cor(HA_df_wide$SwitGPE3rep3,HA_df_wide$SwitGPE3rep1),sep=' '))
write_xlsx(HA_df_wide,'result/HA_EggPassaging_all_wide.xlsx')