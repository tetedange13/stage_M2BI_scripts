#!/usr/bin/env Rscript

library(fmsb)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Input dat file NOT specified !", call.=FALSE)
} else if (length(args) == 1) {
  # default output file
  data_to_load <- args[1]
}

# data=as.data.frame(matrix( sample( 0:20 , 15 , replace=F) , ncol=5))
data <- read.csv(data_to_load, header=T, sep=";", row.names="Zymo_sp")
#colnames(data)=c("math" , "english" , "biology" , "music" , "R-coding" )
#rownames(data) <- data["Zymo_sp"]
min_val <- -10
max_val <- 10
nb_var <- ncol(data)
data <- rbind(rep(max_val, nb_var) , rep(min_val, nb_var) , data)
#print(data)

# Plot 2: Same plot with custom features
colors_border=c( rgb(0.6,0.6,0.6,0.4), rgb(0.2,0.5,0.5,0.9), 
                 rgb(0.8,0.2,0.5,0.9), rgb(0.7,0.5,0.1,0.9))
colors_in=c( rgb(0.6,0.6,0.6,0.4), rgb(0.2,0.5,0.5,0.4), 
             rgb(0.8,0.2,0.5,0.4), rgb(0.7,0.5,0.1,0.4))

test = X11()
par(mar=c(5.1, 4.1, 6.1, 2.1), xpd=TRUE)
radarchart( data , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="dark grey", 
            caxislabels=seq(min_val,max_val,5), cglwd=0.8,
            #custom labels
            vlcex=0.8,
            title="Cusco2018 - Minimap2 - toZymo - COMPLETE RAW (Majo)"
)
legend(x=-1, y=1.6, legend=rownames(data[-c(1,2), ]), bty="n", pch=20, 
       col=colors_in, text.col="dark grey", cex=0.9, pt.cex=3, ncol=2)

while(identical(readLines(), character(0))){Sys.sleep(1)}
