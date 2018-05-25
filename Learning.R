#set working directory
setwd("/Users/sgamoney/Google Drive/mouse_microbiome/")

#load necessary packages -- "library" or "require"
require(tidyverse)
require(data.table)

#load data
abun_raw <- fread("2016_10_13_run2_otus.csv")
abun <- abun_raw
abun[345,4719] <- 0
abun <- abun %>%
  select(-Run, -Day, -Reads) %>%
  rename(Cage = `Cage#`) %>%
  gather(otuid, abun, -SampleID, -Cage, - Abx, -Dose, -Hour, -Cohoused, -Pair) %>%
  group_by(SampleID) %>%
  mutate(abun = abun/sum(abun))

#demo plots
tempz <- abun %>%
  group_by(otuid, Hour) %>%
  summarise(abun = sum(abun)) %>%
  group_by(otuid) %>%
  mutate(total = sum(abun)) %>%
  filter(total > 7)

ggplot(tempz, aes(x = Hour, y = abun, color = otuid)) +
  geom_point() +
  stat_smooth() +
  facet_wrap(~otuid, scale = 'free_y')
