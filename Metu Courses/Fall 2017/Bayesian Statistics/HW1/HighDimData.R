
##########################################################
#	                                                   #
#     2-D Random Data by Anil Aksu                       #
#   It is developed to show some graphical methods in R  #
#     						               #
#									   #
##########################################################	

## libraries 
library(gplots)
library(rgl)
library(car)
library(RColorBrewer)
library(plotly)
library(scatterplot3d) 

## random data generation
size <- 100             #length of random number vectors
set.seed(1) 
x <- runif(size, 5.0, 7.5)          # generate samples from uniform distribution (0.0, 1.0)
y <- runif(size, 5.0, 7.5)
z <- runif(size, 5.0, 7.5)
df <- matrix( runif(size,1., 10), 10, 10) 

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

x11()
#scatter3d(x,y,z, main="PDF Scatterplot Example")
 plot_ly(mtcars, x = ~wt, y = ~hp, z = ~qsec, color = ~am, colors = c('#BF382A', '#0C4B8E'))
x11()
heatmap.2(
  df
  , Colv = T
  , Rowv = T
  , trace="none"
  , col = colorRampPalette(c('red', 'yellow'))(12)
  , labRow=F
  , labCol=F
  , dendrogram="none"
  , margins = c(5,5)
  , main = "heatmap representation"
)

# Let's use the car dataset proposed by R
data=mtcars
 
# Make the plot
my_colors <- brewer.pal(nlevels(as.factor(data$cyl)), "Set2")

x11()
scatterplotMatrix(~mpg+hp+drat|gear, data=data , reg.line="" , smoother="", col=my_colors , smoother.args=list(col="grey") , cex=1.5 , pch=c(15,16,17) , main="Scatter plot with Three Gear Options")
x11()
pdf("C:/Users/anil.aksu/Dropbox/Metu Phd Work/Courses/Fall 2017/Bayesian Statistics/HW1/scatter3D.pdf")
scatterplot3d(x,y,z, pch=16, highlight.3d=TRUE, main="3D Scatterplot")
dev.off()
