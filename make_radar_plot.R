#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Input dat file NOT specified !", call.=FALSE)
} else if (length(args) == 1) {
  # default output file
  data_to_load <- args[1]
}

library(fmsb)


data <- read.csv(data_to_load, header=T, sep=";", row.names="Zymo_sp", 
                 comment.char="#")

print(paste("MIN_VAL_DF:", data[which(data == min(data), arr.ind = TRUE)]))
print(paste("MAX_VAL_DF:", data[which(data == max(data), arr.ind = TRUE)]))

data_old <- data
print(data_old)
writeLines("\n\n")

# Convert into differential proportion relative to ref_prop:
data <- data.frame(sweep(data.matrix(data[-c(1), ]), MARGIN=2, 
                         STATS=as.numeric(data[1, ])) * 100)
data <- rbind(expected=rep(0, 9), data)

min_val <- -15
max_val <- 15
nb_var <- ncol(data)
data <- rbind(rep(max_val, nb_var), rep(min_val, nb_var), data)
print(data)
writeLines("\n\n")



# Plot 2: Same plot with custom features
colors_border=c( rgb(0.6,0.6,0.6,0.4), rgb(0.2,0.5,0.5,0.9), 
                 rgb(0.8,0.2,0.5,0.9), rgb(0.7,0.5,0.1,0.9))
colors_in=c( rgb(0.6,0.6,0.6,0.4), rgb(0.2,0.5,0.5,0.4), 
             rgb(0.8,0.2,0.5,0.4), rgb(0.7,0.5,0.1,0.4))

filename <- basename(data_to_load)
title_plot <- substr(filename, 1, regexpr("\\.", filename)[1]-1)

X11()
par(mar=c(5.1, 4.1, 6.1, 2.1), xpd=TRUE)
radarchart( data, axistype=1, 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="dark grey", cglwd=0.8,
            caxislabels=seq(min_val,max_val,length.out=5), 
            #custom labels
            vlcex=0.8,
            title=paste(title_plot, "\n")
)
legend(x=-1, y=1.6, legend=rownames(data[-c(1,2), ]), bty="n", pch=20, 
       col=colors_in, text.col="dark grey", cex=0.9, pt.cex=3, ncol=2)

while(identical(readLines(), character(0))){Sys.sleep(1)}
