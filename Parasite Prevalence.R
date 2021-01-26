library(readxl)
library(ggplot2)
library(dplyr)
paraprev <- read_excel("Desktop/Masters /Writing/Figure.xlsx", 
                       sheet = "fig1")
View(paraprev)

Parasite_Prevalence = paraprev %>%
  arrange(Prevalence) %>%  
  mutate(Parasite = factor(Parasite, levels=c("C. biliophilus", "Unidentified Trematode", "Strongyloides sp.", 
                              "Pinworm (Enterobius sp.)", "Unknown"))) %>%
  ggplot(aes(x = Parasite, y = Prevalence)) + 
      geom_col() +
      #scale_fill_viridis(discrete = T) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      geom_text(aes(label= c("4%", "10%", "19%", "29%", "43%")), vjust = -0.4) 
      
Parasite_Prevalence
