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

test_func <- function(df_twoLines, foo_arg) {
    my_var <- -1 * diff(as.matrix(df_twoLines))
    #my_var <- colSums(df_twoLines)
    #print(sign(my_var)*my_var)
    #stat_test <- fisher.test(df_twoLines, simulate.p.value=T);print(stat_test)

    tot_nb_reads <- sum(df_twoLines[1, ]) # Both supposed to have same nb of reads
    my_vect <- c()


    for (i in 1:ncol(df_twoLines)) {
        tmp_var <- df_twoLines[, i] ; tmp_var2 <- rep(tot_nb_reads, 2)
        #stat_test <- fisher.test(cbind(tmp_var, tmp_var2))
        #print(colnames(df_twoLines[i])); print(stat_test)
        stat_test <- chisq.test(df_twoLines[, i])
        print(stat_test)
        val_to_add <- round(stat_test$p.value, 5)-0.5
        my_vect <- c(my_vect, val_to_add)
        #my_vect <- c(my_vect, exp(1-stat_test$p.value))
        #print(stat_test)
    }

    #print(my_vect)
    #fisher.test(df_twoLines, conf.int=F)
    #return(data.frame(t(my_var)))

    give_sign <- sign(my_var)
    give_sign[1, ] <- as.factor(rep(1, ncol(df_twoLines))) # Sign is not taken
    to_return <- data.frame(give_sign*my_vect)
    #print(to_return)
    return(to_return)
}

#data <- read.csv(data_to_load, header=T, sep=";", row.names="Zymo_sp", 
data <- read.csv(data_to_load, header=T, sep=";", 
                 comment.char="#")
print(data)

grouped <- data %>%
            group_by(run) %>%
            #summarise(felix=sum(.data))
            group_map(test_func)
writeLines("\n\n")
print.data.frame(grouped)

melted <- melt(grouped, "run")

#exit()
print(melted)

filename <- basename(data_to_load)
title_plot <- substr(filename, 1, regexpr("\\.", filename)[1]-1)

# Generate plot:
p <- ggplot(data=melted, aes(x=run, y=variable)) +
  geom_tile(aes(fill=value)) +
  labs(title=title_plot,
             x="", 
             y="Species in MC", fill=" p-values scale") + # 'fill' = legend title
  facet_grid(cols=vars(run), scales="free_x", space="free_x") +
  theme(axis.text.x=element_text(angle=45, hjust=1, face="bold", 
                                       size=12)) +
  theme(plot.title=element_text(hjust = 0.5, face="bold"), 
        strip.background=element_rect(fill="grey"),
        strip.text=element_text(face="bold")) +
  scale_fill_gradient2(low="white", mid="#E69F00", high="#56B4E9") #, mid="white")

ggsave(filename="my_plot.png", plot = p)
quit()

#print(paste("MIN_VAL_DF:", data[which(data == min(data), arr.ind = TRUE)]))
#print(paste("MAX_VAL_DF:", data[which(data == max(data), arr.ind = TRUE)]))


# Convert into differential proportion relative to ref_prop:
#oldLot <- c(0.0, 0.157, 0.104, 0.100, 0.188, 0.159, 0.046, 0.113, 0.133)
#newLot <- c(0.0, 0.174, 0.099, 0.101, 0.184, 0.141, 0.042, 0.104, 0.155) # NEW LOT


