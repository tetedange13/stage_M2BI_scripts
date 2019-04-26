#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(reshape2)

# test  <- data.frame(person=c("A", "B", "C", "D", "E"), 
#                     value1=c(100,150,120,80,150),     
#                     value2=c(25,30,45,30,30) , 
#                     value3=c(100,120,150,150,200))
# print(test)
#melted <- melt(test, "person")
#print(melted$variable == 'value1')

# melted$cat <- ''
# melted[melted$variable == 'value1',]$cat <- "first"
# melted[melted$variable != 'value1',]$cat <- "second"

df2 <- read.csv("my_CSVs/test_gglot2.csv", header=T, sep=';')
melted <- melt(df2, "run")

cond_raw = (melted$variable == 'unmap' | melted$variable == 'prec')
cond_lenFilt = (melted$variable == 'unmap_lenFilt' | 
                 melted$variable == 'prec_lenFilt')
melted$cat <- ''
melted[cond_raw,]$cat <- "1_Raw"
melted[cond_lenFilt,]$cat <- "2_LenFilt"
# melted <- ddply(melted, c("run", "cat"),
#                    transform, label_ypos=value-cumsum(value))
                #transform, label_ypos=cumsum(value))
grouped <- melted %>%
            group_by(run, cat) %>%
            summarise(total = sum(value))
melted <- melted %>% arrange(run, cat)

melted['label_ypos'] <- rep(grouped$total, each=2) - melted$value

X11()
ggplot(melted, aes(x = cat, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'stack') +
  geom_text(aes(y=label_ypos, label=value), vjust=-0.5, 
            color="white", size=3.5) +
  facet_grid(~ run) +
  ggtitle("Mapping (p08N25) of complete to Zymo") + 
  theme(plot.title = element_text(hjust = 0.5))

while(identical(readLines(), character(0))){Sys.sleep(1)}
