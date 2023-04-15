library(phyloseq)
library(tidyverse)
library(here)
library(readxl)

metadata <- readxl::read_xlsx("C:/Users/armetcal/Downloads/Other People's Work/Paula Data/MBI metagenomics part 2/Metagenomics 2_metadata.xlsx", sheet = 'Sheet1') %>% 
  mutate(Diet = ifelse(str_detect(`Sample ID`,'CON'),'CON','LM'),
         Day = ifelse(str_detect(`Sample ID`,'D0'),'Day 0','Day 28')) %>% 
  mutate(Group = paste(Diet,Day,sep=' ')) %>% 
  mutate(Group = factor(Group, levels = c('CON Day 0','LM Day 0','CON Day 28','LM Day 28'))) %>% 
  rename(Sample = `Sample ID`) %>% 
  select(Sample, Diet, Day, Group, everything())

# FUNCTIONAL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load excel data, extract table
df = readxl::read_xlsx("C:/Users/armetcal/Downloads/Other People's Work/Paula Data/MBI metagenomics part 2/Functnional excel files/SUPER-FOCUS/output_all_levels_and_function.xlsx")
names(df) = df[4,]
df = df[5:nrow(df),]

# Select non-normalized data
df = df %>% select(-contains('%'))

names(df)[5:ncol(df)] = sapply(names(df)[5:ncol(df)], function(x) {
  w = which(metadata$MDB_ID == str_sub(x,end=10))
  return(metadata$Sample[w])
})

# Convert to numeric
df = df %>% mutate_at(vars(contains('LM'),contains('CON')),as.numeric)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(df, 'Data Files/functional_abundances2023.csv', row.names = F)
write.csv(metadata, 'Data Files/metadata2023.csv', row.names = F)
