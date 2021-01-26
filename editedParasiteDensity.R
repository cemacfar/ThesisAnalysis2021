library(readxl)
library(ggplot2)
library(dplyr)
library(wesanderson)
paralocaldens <- read_excel("Desktop/Masters /Writing/Figure.xlsx", 
                            sheet = "fig2edit")
View(paralocaldens)

Parasite_Density = paralocaldens %>%
  arrange(AveDens) %>%  
  mutate(Location = factor(Location, levels=c("Santa Rosa (n = 13)", "Fragmented (n = 8)"))) %>%
  ggplot(aes(x = Location, y = AveDens)) + 
                      geom_col(aes(fill = Location)) +
                      ylab("Density (eggs/mL)") +
  scale_fill_manual(values = wes_palette("Darjeeling2")) +
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin=AveDens-SE, ymax=AveDens+SE), width=.05); Parasite_Density


