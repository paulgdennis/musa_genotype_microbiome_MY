# Functions

# Paul Dennis 

#install.packages("vegan")
#install.packages("sciplot")
#install.packages("tidyr")
#install.packages("emmeans")
#install.packages("multcomp")
#install.packages("lsmeans")
#install.packages("summarytools")


library(vegan)
library(sciplot)
library(tidyr)
library(emmeans)
library(multcomp)
library(lsmeans)
library(summarytools)
library(ggnested)
library(RColorBrewer)
#library(inlmisc)

## Axis percent 
# Returns the percentage varation on ordination axes

axis.percent <- function(ordination){
  round((100*eigenvals(ordination)[1:2]/ordination$tot.chi[[1]]),digits=2)
}

# END

# Add transparency from https://stackoverflow.com/questions/12995683/any-way-to-make-plot-points-in-scatterplot-more-transparent-in-r
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

# END

