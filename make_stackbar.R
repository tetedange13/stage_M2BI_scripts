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



df2 <- read.csv(data_to_load, header=T, sep=';', comment.char="#")

# Check which mode (between 'reproduction of Cusco res' and 'normal') we are:
ncol_less_one <- ncol(df2) - 1
reprod_cusco <- F
if (ncol_less_one%%9 == 0) {
    reprod_cusco <- T
}

if (!reprod_cusco) {
    stopifnot(ncol(df2)%%2 == 1)
}
print(df2)
writeLines("\n")

# Remove 'topOne' cols for report:
rm_topOne = T
if (rm_topOne) {
    df2$sens_2.topOne <- NULL ; df2$FDR_2.topOne <- NULL
    # df2$sens_4.toNCBI16S <- NULL ; df2$FDR_4.toNCBI16S <- NULL
    writeLines("Warning: removed the 'topOne' column from the input df\n")
}

melted <- melt(df2, "run")


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


# Colors palette colorblind-friendly:
my_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
# my_palette <- rev(my_palette) # Reverse order, to have pink and orange 1st


if (!reprod_cusco) {
    # Safeguards:
    stopifnot(length(levels(melted$variable)) == 2)
    stopifnot(length(levels(melted$cat)) == (ncol(df2)-1)/2)
} else {
    # To reproduce data plots from Cusco_2018:
    palette_cusco <- c("#A6CEE3", "#569DA3", "#51AE41", "#F68A89", "#F06C44", 
                       "#F6860F", "#B295C8", "#C7B69A", "#B15928")
    oldLot <- c(0.0, 0.157, 0.104, 0.100, 0.188, 0.159, 0.046, 0.113, 0.133)
    newLot <- c(0.0, 0.174, 0.099, 0.101, 0.184, 0.141, 0.042, 0.104, 0.155) # NEW LOT
    tmp_df <- cbind.data.frame(rep("Mock_db", 9), unique(melted$variable), 
                               oldLot, rep("ref", 9))
    colnames(tmp_df) <- colnames(melted)
    melted <- rbind(tmp_df, melted)
    tmp_df2 <- cbind.data.frame(rep("Mock_db", 9), unique(melted$variable), 
                               newLot, rep("ref_new", 9))
    colnames(tmp_df2) <- colnames(melted)
    melted <- rbind(melted, tmp_df2)

    my_palette <- palette_cusco
}



# Need to be descendant with 'variable' col, to have 'sens' BEFORE 'FDR':
melted <- melted %>% arrange(run, cat, desc(variable))

grouped <- melted %>%
            group_by(run, cat) %>%
            do(total=cumsum(.$value)-0.5*.$value)

melted['label_ypos'] <- unlist(grouped$total)

filename <- basename(data_to_load)
title_plot <- substr(filename, 1, regexpr("\\.", filename)[1]-1)


print(melted)
p <- ggplot(melted, aes(x=cat, y=value, fill=variable)) +
        geom_bar(stat='identity', position='stack') +
        facet_grid(cols=vars(run), scales="free_x", space="free_x") + # No x-tick for missing + bars same widt
        scale_fill_manual(values=my_palette[1:length(levels(melted$variable))]) + 
        geom_text(aes(y=label_ypos, label=round(value, 3), fontface="bold"), 
                  vjust=0.5, color="white", size=3.5) +
        theme(axis.text.x=element_text(angle=45, hjust=1, face="bold", 
                                       size=12)) +
        # labs(title=paste(title_plot, "- lvl SPECIES"), 
        labs(title=title_plot,
             x="\nDifferent conditions", 
             y="Sensitivity + FDR", fill=" Metrics") + # 'fill' = legend title
        theme(plot.title=element_text(hjust = 0.5, face="bold"), 
              panel.background=element_rect(fill="light grey"),
              panel.grid.major=element_blank(),
              strip.background=element_rect(fill="grey"),
              strip.text=element_text(face="bold")) +
        scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1, 1.25))


X11()
print(p)

while(identical(readLines(), character(0))){Sys.sleep(1)}
