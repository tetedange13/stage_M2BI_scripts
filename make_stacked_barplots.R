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
# exit()

# cond_raw <- (melted$variable == 'FDR' | melted$variable == 'sens')
# cond_lenFilt <- (melted$variable == 'FDR_lenFilt' | 
#                  melted$variable == 'sens_lenFilt')
# cond_topOne <- (melted$variable == 'FDR_topOne' | 
#                  melted$variable == 'sens_topOne')

# melted$cat <- ''
# melted[cond_raw, ]$cat <- "1_Centri_Raw"
# melted[cond_lenFilt, ]$cat <- "2_Centri_Filt"
# # melted[cond_topOne, ]$cat <- "3_TopOne"
# melted[cond_lenFilt, ]$variable <- melted[cond_raw, ]$variable

# melted <- ddply(melted, c("run", "cat"),
#                    transform, label_ypos=value-cumsum(value))
                #transform, label_ypos=cumsum(value))

# Need to be descendant with 'variable' col, to have 'sens' BEFORE 'FDR':
melted <- melted %>% arrange(run, cat, desc(variable))
# print(melted)

# select_groups <- function(dd, gr, ...) dd[sort(unlist(attr(dd, "groups")$.rows[ gr ])), ]
grouped <- melted %>%
            group_by(run, cat) %>%
            # do(total=cumsum(.$value)-0.5*.$value)
            do(total=cumsum(.$value)-0.5*.$value)
            # 
# print.data.frame(grouped)

# melted['label_ypos'] <- (rep(grouped$total, each=2) - melted$value)/2
melted['label_ypos'] <- rep(unlist(grouped$total))

print(melted)
filename = basename(data_to_load)
title_plot = substr(filename, 1, regexpr("\\.", filename)[1]-1)

X11()
ggplot(melted, aes(x = cat, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'stack') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(y=label_ypos, label=value), vjust=0.5, 
            color="white", size=3.5) +
  facet_grid(~ run) +
  ggtitle(paste(title_plot, "- lvl SPECIES")) + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "grey"),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(face="bold")) +
  scale_y_continuous(name="recall + FDR", 
                     breaks=c(0, 0.25, 0.5, 0.75, 1))

while(identical(readLines(), character(0))){Sys.sleep(1)}
