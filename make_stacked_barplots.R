#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Input dat file NOT specified !", call.=FALSE)
} else if (length(args) == 1) {
  # default output file
  data_to_load <- args[1]
}

library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

# test  <- data.frame(person=c("A", "B", "C", "D", "E"), 
#                     value1=c(100,150,120,80,150),     
#                     value2=c(25,30,45,30,30) , 
#                     value3=c(100,120,150,150,200))
#melted <- melt(test, "person")
#print(melted$variable == 'value1')
# melted$cat <- ''
# melted[melted$variable == 'value1',]$cat <- "first"
# melted[melted$variable != 'value1',]$cat <- "second"

df2 <- read.csv(data_to_load, header=T, sep=';')
melted <- melt(df2, "run")
# print(melted)

test <- str_locate(melted$variable, '_')[, "start"]
vect_debut <- c()
vect_end <- c() 
for (i in 1:length(melted$variable)) {
    cond_name <- toString(melted[i, "variable"])
    vect_debut <- c(vect_debut, 
                 substr(cond_name, start=1, 
                        stop=test[i]-1))
    vect_end <- c(vect_end, 
                 substr(cond_name, start=test[i]+1, 
                        stop=nchar(cond_name)))
}

melted$cat <- factor(vect_end)
melted$variable <- factor(vect_debut)
# print(melted)


# Need to be descendant with 'variable' col, to have 'sens' BEFORE 'FDR':
melted <- melted %>% arrange(run, cat, desc(variable))
# print(melted)

grouped <- melted %>%
            group_by(run, cat) %>%
            do(total=cumsum(.$value)-0.5*.$value)
# print.data.frame(grouped)

melted['label_ypos'] <- unlist(grouped$total)

print(melted)
filename = basename(data_to_load)
title_plot = substr(filename, 1, regexpr("\\.", filename)[1]-1)

X11()
ggplot(melted, aes(x=cat, y=value, fill=variable)) +
  geom_bar(stat='identity', position='stack') +
  theme(axis.text.x=element_text(angle=90, hjust=1, face="bold")) +
  geom_text(aes(y=label_ypos, label=value), 
            vjust=0.5, color="white", size=3.5) +
  facet_grid(~ run) +
  labs(title=paste(title_plot, "- lvl SPECIES"), x="\nDifferent conditions", 
       y="Recall + FDR", fill=" Metrics") + # 'fill' = legend title
  theme(plot.title=element_text(hjust = 0.5), 
        panel.background=element_rect(fill="grey"),
        panel.grid.major=element_blank()) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1))

while(identical(readLines(), character(0))){Sys.sleep(1)}
