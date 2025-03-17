data = read.csv("anderson_2016_mt/data/quadrat_data/species_list.csv")
data = data[,c(2:6)]
head(data)
data <- data %>%
  mutate(species = paste(species, X, sep = " "))
data = data[,c(1,3:6)]
head(data)

write.csv(data, "anderson_2016_mt/data/quadrat_data/species_list.csv", row.names = F)
