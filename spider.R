# Library
library(fmsb)

# Create data: note in High school for Jonathan:
data <- as.data.frame(matrix(c(4.7,  3.2, 3.4, 2.6, 4.1, 4), ncol=6))
colnames(data) <- c("math" , "english" , "biology" , "music" , "R-coding", "data-viz")

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
data <- rbind(rep(5,10) , rep(0,10) , data)

# Check your data, it has to look like this!
# head(data)
png("/home/colin/Documents/aberdeen/codes/Specialist_generalist/oupidou.png", height = 1000, width = 1000, bg = "transparent")
# Custom the radarChart !
radarchart( data  , axistype=1 , 
            
            #custom polygon
            pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5) , plwd=8 , 
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,5,1), cglwd=2,
            
            #custom labels
            vlcex=0.00001
)

dev.off()