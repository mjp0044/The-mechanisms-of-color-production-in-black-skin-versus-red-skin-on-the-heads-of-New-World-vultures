#Vulture skin and blood
#load preferred packages
library(Hmisc)
library(tidyverse)
library(lme4)
library(cowplot)
library(lmerTest)
library(agricolae)
library(MASS)
library(emmeans)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(car)
library(sf)
library(ggjoy)
library(MuMIn)
library(plotrix)
library(HH)
library(rstatix)

#optional plotting packages
library(betareg)
library(gifski)
library(Rcpp)
library(gganimate)
#Set global theme for plots
theme_set(theme_cowplot())


####Begin plotting of HPLC chromatograms####

datum = read.csv(file = "hplc.chromatograms.vulture.csv")

#Trim first 10 minutes off to remove solvent front peak and reduce white space of graph
#trim off last 7 minutes to remove equilibration period and reduce white space of graph
#Untrimmed data in original data sheet imported above (HPLC_chromatograms.csv)
datumtrim = subset(datum, time >=10 & time <=25)

#You will notice numbers after the Y variable call for each graph. This is just an adjustment to 
#all of the values of the intensity column for that standard to correct for the natural baseline drift
#of the HPLC system



#Sample R23 with standards

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,2,3,4,5)]
#Adjust y values for baseline drift during hplc
dd_sub$lutein = dd_sub$lutein-150
dd_sub$zeaxanthin = dd_sub$zeaxanthin-60
dd_sub$bcarotene = dd_sub$bcarotene

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR23standards <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line:",values=c("goldenrod3","blue2", "green3", "black"),
                      labels=c("Lutein","Zeaxanthin", "Beta-carotene", "R23"))+
  scale_y_continuous(breaks=seq(0,4000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("Relative intensity (mV)")+ xlab("Time (min)")+
  annotate("text", x = 20, y = 100, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 100, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 100, label = "5", size = 4, color ="black")+
  annotate("text", x = 17, y = 2000, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.8, y = 3000, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 4000, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 4000, label = "Sample R23
Hatch Year Male", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
#View plot
pR23standards

jpeg(filename = "Sample R23 with standards.jpg", width = 6, height = 6, units = "in", res = 500)
pR23standards
dev.off()


#Sample R23 

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,5)]


##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR23 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black"),
                      labels=c("Skin"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("Relative intensity (mV)")+ xlab("")+
  annotate("text", x = 20, y = 200, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 200, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 200, label = "5", size = 4, color ="black")+
  annotate("text", x = 16.9, y = 200, label = "1",size = 4, color ="black")+
  annotate("text", x = 18, y = 500, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 2600, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R23
HY Male", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR23




#Sample R20

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,6,16)]
#Adjust y values for baseline drift during hplc
dd_sub$R20 = dd_sub$R20-80

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR20 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("")+
  annotate("text", x = 20, y =200, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 200, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 300, label = "5", size = 4, color ="black")+
  annotate("text", x = 17, y = 300, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.5, y = 900, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 2000, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R20
ASY Male", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR20


#Sample R21

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,7,17)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR21 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("")+
  annotate("text", x = 20, y =100, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 100, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 200, label = "5", size = 4, color ="black")+
  annotate("text", x = 17, y = 600, label = "1",size = 4, color ="black")+
  annotate("text", x = 18.1, y = 700, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 1000, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R21
HY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR21


#Sample R24

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,8,19)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR24 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("")+
  annotate("text", x = 20, y =300, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 400, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 600, label = "5", size = 4, color ="black")+
  annotate("text", x = 16.9, y = 800, label = "1",size = 4, color ="black")+
  annotate("text", x = 18.1, y = 900, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 5700, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R24
ASY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR24


#Sample R25

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,9,20)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR25 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("")+
  annotate("text", x = 20, y =300, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 400, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 800, label = "5", size = 4, color ="black")+
  annotate("text", x = 16.9, y = 1300, label = "1",size = 4, color ="black")+
  annotate("text", x = 18.1, y = 1500, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 4700, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R25
ASY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR25



#Sample R59

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,10,21)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR59 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("Relative intensity (mV)")+ xlab("Time (min)")+
  annotate("text", x = 20, y =100, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 100, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 100, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 500, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 500, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 900, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R59
SY Male", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR59


#Sample R09

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,11,22)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR09 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (min)")+
  annotate("text", x = 20, y =100, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 100, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 100, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 300, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 300, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 700, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R09
ASY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR09


#Sample R08

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,12,23)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR08 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (min)")+
  annotate("text", x = 20, y =200, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 200, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 200, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 1800, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 1900, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 4200, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R08
ASY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR08


#Sample R10

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,13,24)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR10 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line:",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (min)")+
  annotate("text", x = 20, y =300, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 300, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 300, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 1800, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 1900, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 5200, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R10
ASY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = c(-0.6,-0.18), legend.direction = "horizontal", 
        legend.text = element_text(size=16), legend.title = element_text(size =16))+
  theme(plot.margin = unit(c(0.3,0,0.65,0), "cm"))
#View plot
pR10


#Sample R11

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,14,25)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR11 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (min)")+
  annotate("text", x = 20, y =400, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 400, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 400, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 3500, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 3600, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 5800, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R11
ASY Male", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR11


#Sample R12

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,15)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR12 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black"),
                      labels=c("Skin"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (min)")+
  annotate("text", x = 20, y =100, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 100, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 100, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 200, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 200, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 1300, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R12
HY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR12



#Sample R22

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,18)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR22 <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("red1"),
                      labels=c("Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (min)")+
  annotate("text", x = 20, y =100, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 100, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 100, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 1900, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 2000, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 700, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R22
HY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR22


#Make combined figure for all chromatograms
figure1 <- ggarrange(pR23, pR21, pR20, pR24, pR25, pR22, 
                     pR59, pR09, pR08, pR10, pR11, pR12,
                     labels = c("Turkey", "", "", "", "", "", "Black","","","", "", ""), 
                     font.label = list(size=12), hjust = -0.1,vjust = 1, 
                     ncol=6, nrow=2, common.legend = FALSE)

jpeg(file="Summary_same_Yaxis.jpg", units="in", width=23, height=8, res=500)
figure1
dev.off()



#### Plot with individually scaled y axis ####

#Sample R23 

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,5)]


##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR23ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black"),
                      labels=c("Skin"))+
  scale_y_continuous(breaks=seq(0,3000,600))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("Relative Intensity (mV)")+ xlab("")+
  annotate("text", x = 20, y = 100, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 100, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 200, label = "5", size = 4, color ="black")+
  annotate("text", x = 17, y = 200, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 500, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 2500, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 3000, label = "Sample R23
HY Male", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR23ind




#Sample R20

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,6,16)]
#Adjust y values for baseline drift during hplc
dd_sub$R20 = dd_sub$R20-80

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR20ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,3000,600))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("")+
  annotate("text", x = 20, y =200, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 200, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 300, label = "5", size = 4, color ="black")+
  annotate("text", x = 17, y = 300, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.5, y = 900, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 1900, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 3000, label = "Sample R20
ASY Male", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR20ind


#Sample R21

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,7,17)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR21ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,1000,200))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("")+
  annotate("text", x = 20, y =100, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 100, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 200, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.2, y = 600, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 700, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 900, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 1000, label = "Sample R21
HY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR21ind


#Sample R24

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,8,19)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR24ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("")+
  annotate("text", x = 20, y =300, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 400, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 600, label = "5", size = 4, color ="black")+
  annotate("text", x = 17, y = 800, label = "1",size = 4, color ="black")+
  annotate("text", x = 18, y = 900, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 5700, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R24
ASY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR24ind


#Sample R25

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,9,20)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR25ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,5000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("")+
  annotate("text", x = 20, y =300, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 400, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 800, label = "5", size = 4, color ="black")+
  annotate("text", x = 17, y = 1300, label = "1",size = 4, color ="black")+
  annotate("text", x = 18, y = 1500, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 4600, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 5000, label = "Sample R25
ASY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR25ind



#Sample R59

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,10,21)]
#Adjust y values for baseline drift during hplc
dd_sub$R59 = dd_sub$R59+40


##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR59ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(-100,600,100))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("Relative Intensity (mV)")+ xlab("Time (min)")+
  annotate("text", x = 20, y =30, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 30, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 30, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 500, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 500, label = "2", size = 4, color ="black")+
  annotate("text", x = 24.5, y = 600, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 600, label = "Sample R59
SY Male", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR59ind


#Sample R09

  ##Subset the necessary columns (2-4 are standards while first and last column numbers are 
  #feather samples changed as needed)
  dd_sub = datumtrim[,c(1,11,22)]
  #Adjust y values for baseline drift during hplc
  dd_sub$R09 = dd_sub$R09-20
  
  ##Then rearrange your data frame
  dd = melt(dd_sub, id=c("time"))
  
  #Make plot
  pR09ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
    scale_colour_manual(name="Line",values=c("black", "red1"),
                        labels=c("Skin", "Blood"))+
    scale_y_continuous(breaks=seq(-160,400,80))+
    scale_x_continuous(breaks=seq(10,25,3))+
    ylab("")+ xlab("Time (min)")+
    annotate("text", x = 20, y =-30, label = "3", size = 4, color ="black")+
    annotate("text", x = 21.5, y = -30, label = "4", size = 4, color ="black")+
    annotate("text", x = 23, y = -30, label = "5", size = 4, color ="black")+
    annotate("text", x = 17.1, y = 200, label = "1",size = 4, color ="black")+
    annotate("text", x = 17.9, y = 200, label = "2", size = 4, color ="black")+
    annotate("text", x = 24, y = 400, label = "6", size = 4, color ="black")+
    annotate("text", x = 12, y = 400, label = "Sample R09
  ASY Female", size = 3.5, color ="black")+
    guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
    theme(legend.position = "none", legend.direction = "horizontal")+
    theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
  #View plot
  pR09ind


#Sample R08

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,12,23)]
#Adjust y values for baseline drift during hplc
dd_sub$R08 = dd_sub$R08-30
##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR08ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,4000,800))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (min)")+
  annotate("text", x = 20, y =200, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 200, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 200, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 1800, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 1900, label = "2", size = 4, color ="black")+
  annotate("text", x = 24.5, y = 4000, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 4000, label = "Sample R08
ASY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR08ind


#Sample R10

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,13,24)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR10ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line:",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,5000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (Min)")+
  annotate("text", x = 20, y =300, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 300, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 300, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 1800, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 1900, label = "2", size = 4, color ="black")+
  annotate("text", x = 24.5, y = 5000, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 5000, label = "Sample R10
ASY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = c(-0.6,-0.18), legend.direction = "horizontal", 
        legend.text = element_text(size=16), legend.title = element_text(size =16))+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR10ind


#Sample R11

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,14,25)]
#Adjust y values for baseline drift during hplc
dd_sub$R11 = dd_sub$R11+40
##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR11ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black", "red1"),
                      labels=c("Skin", "Blood"))+
  scale_y_continuous(breaks=seq(0,6000,1000))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (min)")+
  annotate("text", x = 20, y =400, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 400, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 400, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 3500, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 3600, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 5800, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 6000, label = "Sample R11
ASY Male", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR11ind


#Sample R12

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,15)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR12ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("black"),
                      labels=c("Skin"))+
  scale_y_continuous(breaks=seq(0,1000,200))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("Time (min)")+
  annotate("text", x = 20, y =-10, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = -10, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = -10, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 50, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 50, label = "2", size = 4, color ="black")+
  annotate("text", x = 24.5, y = 1000, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 1000, label = "Sample R12
HY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.3,0,0.7,0), "cm"))
#View plot
pR12ind



#Sample R22

##Subset the necessary columns (2-4 are standards while first and last column numbers are 
#feather samples changed as needed)
dd_sub = datumtrim[,c(1,18)]
#Adjust y values for baseline drift during hplc

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))

#Make plot
pR22ind <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="Line",values=c("red1"),
                      labels=c("Blood"))+
  scale_y_continuous(breaks=seq(0,2000,400))+
  scale_x_continuous(breaks=seq(10,25,3))+
  ylab("")+ xlab("")+
  annotate("text", x = 20, y =100, label = "3", size = 4, color ="black")+
  annotate("text", x = 21.5, y = 100, label = "4", size = 4, color ="black")+
  annotate("text", x = 23, y = 100, label = "5", size = 4, color ="black")+
  annotate("text", x = 17.1, y = 1900, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.9, y = 1900, label = "2", size = 4, color ="black")+
  annotate("text", x = 24, y = 600, label = "6", size = 4, color ="black")+
  annotate("text", x = 12, y = 2000, label = "Sample R22
HY Female", size = 3.5, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 1, ncol = 7))+
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
#View plot
pR22ind


#Make combined figure for all chromatograms
figure2 <- ggarrange(pR23ind, pR21ind, pR20ind, pR24ind, pR25ind, pR22ind, 
                     pR59ind, pR09ind, pR08ind, pR10ind, pR11ind, pR12ind,
                     labels = c("Turkey", "", "", "", "", "", "Black","","","", "", ""), 
                     font.label = list(size=12), hjust = -0.1,vjust = 1, 
                     ncol=6, nrow=2, common.legend = FALSE)

jpeg(file="Summary_scaled_Yaxis.jpg", units="in", width=23, height=8, res=500)
figure2
dev.off()



#### Stats and concentration comparisons ####

#Read in data and specify cells marked with "trace", "ND", and "NA" to be read in as "empties"
datum2 = read.csv(file = "HPLC.data.csv", na.strings = c("Trace", "ND", "NA"))
#Check column codings
str(datum2)

#combine tissue and species categories into label for plotting and analysis
datum2$tissue_and_species = paste(datum2$Species,datum2$Tissue)

#combine tissue and age categories into label for plotting and analysis
datum2$tissue_and_age = paste(datum2$Age.Broad,datum2$Tissue)

#combine tissue and sex categories into label for plotting and analysis
datum2$tissue_and_sex = paste(datum2$Sex,datum2$Tissue)

#combine species and age categories into label for plotting and analysis
datum2$species_and_age = paste(datum2$Age.Broad,datum2$Species)

#combine species and sex categories into label for plotting and analysis
datum2$species_and_sex = paste(datum2$Sex,datum2$Species)

#Melt data frame into lengthwise order
datum2.molten <- melt(datum2[,c(4,5,7,8,10, 23,24,25,26,27)], 
                      id = c("Sample.ID", "Tissue", "Age.Broad", "Species", "Sex"))

#Relevel treatment variable to put in the desired order for plotting
datum2.molten$variable<- factor(datum2.molten$variable,
                         levels = c("Total.st", "Zeaxanthin.st", "Bcarotene.st",
                                    "Lutein.st", "minor.peaks.st"))
#check order with table command
table(datum2.molten$variable)

#Rename carotenoid variable levels
datum2.molten %>% 
  mutate(variable = fct_recode(variable, "Total carotenoids" = "Total.st" ,
                          "Zeaxanthin" = "Zeaxanthin.st", "Beta Carotene" = "Bcarotene.st", "Lutein" = "Lutein.st",
                          "Minor carotenoids" = "minor.peaks.st")) -> datum2.molten

datum2$Total.st.tf = log(datum2$Total.st)
datum2$Zeaxanthin.st.tf = log(datum2$Zeaxanthin.st)
datum2$Bcarotene.st.tf = log(datum2$Bcarotene.st)
datum2$Lutein.st.tf = log(datum2$Lutein.st)
datum2$minor.peaks.st.tf = log(datum2$minor.peaks.st)

####Analysis of overall effect of each variable controlling for other effects of all variables####

#Overall effect of age, sex, tissue, and species on minor carotenoid concentrations
mod.overall.aov = aov(Total.st.tf~Tissue*Species*Age.Broad*Sex, data = datum2)
anova_summary(mod.overall.aov)
summary(mod.overall.aov)

#Overall effect of age, sex, tissue, and species on zeaxanthin proportion
mod.overall.zea = aov(Zeaxanthin.st.tf~Tissue*Species*Age.Broad*Sex, data = datum2)
anova_summary(mod.overall.zea)
summary(mod.overall.zea)

#Overall effect of age, sex, tissue, and species on bcarotene proportion
mod.overall.bcaro = aov(Bcarotene.st.tf~Tissue*Species*Age.Broad*Sex, data = datum2)
anova_summary(mod.overall.bcaro)
summary(mod.overall.bcaro)

#Overall effect of age, sex, tissue, and species on lutein proportion
mod.overall.lutein = aov(Lutein.st.tf~Tissue*Species*Age.Broad*Sex, data = datum2)
anova_summary(mod.overall.lutein)
summary(mod.overall.lutein)

#Overall effect of age, sex, tissue, and species on minor carotenoid proportion
mod.overall.minor = aov(minor.peaks.st.tf~Tissue*Species*Age.Broad*Sex, data = datum2)
anova_summary(mod.overall.minor)
summary(mod.overall.minor)

#There seems to be an overall effect of species on the total carotenoid content and all individual carotenoids except lutein.
#Sex is not significant for total carotenoids or any specific carotenoid
#Effect of age seems to be important for total carotenoid content, driven by effect on zeaxanthin 
#There is a significant interaction between tissue and species for total carotenoids, zeaxanthin, and bcarotene
#There's an interaction between age, tissue, and species in total carotenoid content, driven by effect on zeaxanthin
#All variables seem to have a significant effect on minor carotenoid concentrations, 
#but this is probably due to noise from low concentrations

#Need to include random effect of tissue for age and sex analysis and random effect of age for species analysis
#Block age and sex by species since it has a large effect

####Analysis by species####


#Run linear model for total carotenoids
mod.1 = lmer(Total.st.tf~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.1)
#Rerun without normalization to get means
mod.1 = lmer(Total.st~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.1)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.1, pairwise ~ tissue_and_species) 

#Run linear model for Zeaxanthin
mod.2 = lmer(Zeaxanthin.st.tf~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.2)
#Rerun without normalization to get means
mod.2 = lmer(Zeaxanthin.st~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.2)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.2, pairwise ~ tissue_and_species) 

#Run linear model for beta carotene
mod.3 = lmer(Bcarotene.st.tf~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.3)
#Rerun without normalization to get means
mod.3 = lmer(Bcarotene.st~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.3)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.3, pairwise ~ tissue_and_species) 

#Run linear model for lutein
mod.4 = lmer(Lutein.st.tf~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.4)
#Rerun without normalization to get means
mod.4 = lmer(Lutein.st~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.4)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.4, pairwise ~ tissue_and_species) 

#Run linear model for minor carotenoids
mod.5 = lmer(minor.peaks.st.tf~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.5)
#Rerun without normalization to get means
mod.5 = lmer(minor.peaks.st~tissue_and_species + (1|Age.Broad), data = datum2)
summary(mod.5)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.5, pairwise ~ tissue_and_species) 


#build plot
species <-ggplot(data=datum2.molten, aes(x=Species, y=value, fill=Tissue)) + 
  facet_wrap(~variable,  ncol=5, strip.position = "bottom")+
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+
  scale_fill_manual(values = c("red", "grey")) +
  labs(x="", y=bquote('Carotenoid concentration ('~mu*'g'~'g'^-1*')'), fill="Tissue:") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 16), axis.text.x = element_text(),
        axis.text.y = element_text(size=14),
        legend.position = "bottom", strip.placement = 'outside')
species

#output figure
jpeg(file="Carotenoid concentrations per tissue and species.jpg", units="in", width=12, height=8, res=300)
species
dev.off()


####Analysis by age####

#Run linear model for total carotenoids
mod.6 = lmer(Total.st.tf~species_and_age + (1|Tissue), data = datum2)
summary(mod.6)
#Rerun without normalization to get means
mod.6 = lmer(Total.st~species_and_age + (1|Tissue), data = datum2)
summary(mod.6)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.6, pairwise ~ species_and_age) 

#Run linear model for Zeaxanthin
mod.7 = lmer(Zeaxanthin.st.tf~species_and_age + (1|Tissue), data = datum2)
summary(mod.7)
#Rerun without normalization to get means
mod.7 = lmer(Zeaxanthin.st~species_and_age + (1|Tissue), data = datum2)
summary(mod.7)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.7, pairwise ~ species_and_age) 

#Run linear model for beta carotene
mod.8 = lmer(Bcarotene.st.tf~species_and_age + (1|Tissue), data = datum2)
summary(mod.8)
#Rerun without normalization to get means
mod.8 = lmer(Bcarotene.st~species_and_age + (1|Tissue), data = datum2)
summary(mod.8)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.8, pairwise ~ species_and_age) 

#Run linear model for lutein
mod.9 = lmer(Lutein.st.tf~species_and_age + (1|Tissue), data = datum2)
summary(mod.9)
#Rerun without normalization to get means
mod.9 = lmer(Lutein.st~species_and_age + (1|Tissue), data = datum2)
summary(mod.9)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.9, pairwise ~ species_and_age) 

#Run linear model for minor carotenoids
mod.10 = lmer(minor.peaks.st.tf~species_and_age + (1|Tissue), data = datum2)
summary(mod.10)
#Rerun without normalization to get means
mod.10 = lmer(minor.peaks.st~species_and_age + (1|Tissue), data = datum2)
summary(mod.10)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.10, pairwise ~ species_and_age) 


#build plot
age <-ggplot(data=datum2.molten, aes(x=Species, y=value, fill=Age.Broad)) + 
  facet_wrap(~variable,  ncol=5, strip.position = "bottom")+
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+
  scale_fill_manual(values = c("blue", "lightblue")) +
  labs(x="", y=bquote('Carotenoid concentration ('~mu*'g'~'g'^-1*')'), fill="Age:") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 16), axis.text.x = element_text(),
        axis.text.y = element_text(size=14),
        legend.position = "bottom", strip.placement = 'outside')
age

#output figure
jpeg(file="Carotenoid concentrations per species and age.jpg", units="in", width=12, height=8, res=300)
age
dev.off()




####Analysis by sex####

#Run linear model for total carotenoids
mod.11 = lmer(Total.st.tf~species_and_sex + (1|Tissue), data = datum2)
summary(mod.11)
#Rerun without normalization to get means
mod.11 = lmer(Total.st~species_and_sex + (1|Tissue), data = datum2)
summary(mod.11)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.11, pairwise ~ species_and_sex) 

#Run linear model for Zeaxanthin
mod.12 = lmer(Zeaxanthin.st.tf~species_and_sex + (1|Tissue), data = datum2)
summary(mod.12)
#Rerun without normalization to get means
mod.12 = lmer(Zeaxanthin.st~species_and_sex + (1|Tissue), data = datum2)
summary(mod.12)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.12, pairwise ~ species_and_sex) 

#Run linear model for beta carotene
mod.13 = lmer(Bcarotene.st.tf~species_and_sex + (1|Tissue), data = datum2)
summary(mod.13)
#Rerun without normalization to get means
mod.13 = lmer(Bcarotene.st~species_and_sex + (1|Tissue), data = datum2)
summary(mod.13)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.13, pairwise ~ species_and_sex) 

#Run linear model for lutein
mod.14 = lmer(Lutein.st.tf~species_and_sex + (1|Tissue), data = datum2)
summary(mod.14)
#Rerun without normalization to get means
mod.14 = lmer(Lutein.st~species_and_sex + (1|Tissue), data = datum2)
summary(mod.14)
#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.14, pairwise ~ species_and_sex) 

#Run linear model for minor carotenoids
mod.15 = lmer(minor.peaks.st.tf~species_and_sex + (1|Tissue), data = datum2)
summary(mod.15)

#Rerun without normalization to get means
mod.15 = lmer(minor.peaks~species_and_sex + (1|Tissue), data = datum2)
summary(mod.15)

#Look at pairwise comparisons from model with tukeys correction
emmeans(mod.15, pairwise ~ species_and_sex) 


#build plot
sex <-ggplot(data=datum2.molten, aes(x=Species, y=value, fill=Sex)) + 
  facet_wrap(~variable,  ncol=5, strip.position = "bottom")+
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+
  scale_fill_manual(values = c("darkgreen", "chartreuse")) +
  labs(x="", y=bquote('Carotenoid concentration ('~mu*'g'~'g'^-1*')'), fill="Sex:") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 16), axis.text.x = element_text(),
        axis.text.y = element_text(size=14),
        legend.position = "bottom", strip.placement = 'outside')
sex

#output figure
jpeg(file="Carotenoid concentrations per species and sex.jpg", units="in", width=12, height=8, res=300)
sex
dev.off()




####Run overall without blocking by tissue or random effects####
mod.16 = lmer(Total.st.tf ~ Age.Broad+Species+Sex+(1|Tissue), data = datum2)
summary(mod.16)

#retrieve means
emmeans(mod.16, pairwise~Age.Broad)
emmeans(mod.16, pairwise~Species)
emmeans(mod.16, pairwise~Sex)


#Does total blood concentrations correlate with total skin concentrations?

#Remove birds with only one tissue
datum3 = datum2[-c(1,10,16),]
#Cast data to reorganize it widthwise
datum4 =dcast(datum3, Sample.ID+Species~Tissue, value.var  = 'Total.st')

#Run model 
mod.17 = lm(Skin~Blood, data = subset(datum4, Species == "Black"))
summary(mod.17)
#Get R-squared value for model fit
r.squaredGLMM(mod.17)

#Run model 
mod.17 = lm(Skin~Blood, data = subset(datum4, Species == "Black"))
summary(mod.17)
#Get R-squared value for model fit
r.squaredGLMM(mod.17)

#Run model overall
mod.17 = lmer(Skin~Blood + (1|Species), data = datum4)
summary(mod.17)
#Get R-squared value for model fit
r.squaredGLMM(mod.17)

#Run model for black
mod.17b = lm(Skin~Blood, data = subset(datum4, Species == "Black"))
summary(mod.17b)
#Get R-squared value for model fit
r.squaredGLMM(mod.17b)

#Run model for turkey
mod.17c = lm(Skin~Blood, data = subset(datum4, Species == "Turkey"))
summary(mod.17c)
#Get R-squared value for model fit
r.squaredGLMM(mod.17c)

#plot scatterplot
bloodvsskin <-ggplot(datum4, aes(x=Blood, y=Skin, col=Species)) + 
  geom_point(size =5)+
  stat_smooth(method = "lm")+
  stat_smooth(method="lm", col="black")+
  scale_color_manual(values = c("blue", "red"))+
  #scale_y_continuous(name = "", breaks = 0:1, labels=c("Not Expressed", "Expressed"))+
  labs(x=bquote('Blood Carotenoid concentration ('~mu*'g'~'g'^-1*')'), 
       y=bquote('Skin Carotenoid concentration ('~mu*'g'~'g'^-1*')'))+
  coord_cartesian(ylim = c(-1, 2))+
  annotate("text", x = 0.6, y = 1.8, label = expression(italic(beta)*"=0.38"~italic(p)*"=0.490"), size = 4, color ="black")+
  annotate("text", 0.6, 1.75, label =  bquote('R'^2*'=0.018'), size=3.5, color ="black")+
  annotate("text", x = 0.6, y = 1.6, label = expression(italic(beta)*"=0.21"~italic(p)*"=0.017"), size = 4, color ="blue")+
  annotate("text", 0.6, 1.55, label =  bquote('R'^2*'=0.855'), size=3.5, color ="blue")+
  annotate("text", x = 0.6, y = 2, label = expression(italic(beta)*"=1.03"~italic(p)*"=0.636"), size = 4, color ="red")+
  annotate("text", 0.6, 1.95, label =  bquote('R'^2*'=0.092'), size=3.5, color ="red")
bloodvsskin

#output figure
jpeg(file="Blood vs skin scatter.jpg", units="in", width=6, height=6, res=300)
bloodvsskin
dev.off()


#### Proportions analysis####

#Make percentage columns for plotting
datum2$zea.perc=datum2$zea.prop*100
datum2$bcaro.perc=datum2$bcaro.prop*100
datum2$lutein.perc=datum2$lutein.prop*100
datum2$acrypto.perc=datum2$acrypto.prop*100
datum2$bcrypto.perc=datum2$bcrypto.prop*100
datum2$echine.perc=datum2$echine.prop*100


#Broad divisions

##Species donut plot
#Melt data frame into lengthwise order
datum2.molten.species <- melt(datum2[,c(10, 43,44,45,46,47,48)], id = c("Species"))
#Replace NA's with 0 since only trace or undetectable concentrations
datum2.molten.species[is.na(datum2.molten.species)] = 0

datum2.molten.species %>%
  group_by(Species, variable) %>%
  summarize(mean.conc = mean(value)) -> datum2.molten.species2

#Round for plotting
datum2.molten.species2$mean.conc=round(datum2.molten.species2$mean.conc, digits=0)

#Make plot
species.donut <- ggplot(datum2.molten.species2, aes(x = Species, y = mean.conc, fill = variable, color=Species)) +
  geom_col(size=1.1)+ 
  scale_x_discrete(limits = c(" ", "Black","Turkey"))+
  geom_text(aes(label = mean.conc),
            position = position_stack(vjust = 0.5), col="black", size=4, vjust=0) +
  scale_fill_manual(values = c("Grey96", "Grey85", "Grey75", "Grey65", "Grey55", "Grey45"), labels = c("Zeaxanthin", "B-carotene", "Lutein", "A-crytpoxanthin", "B-cryptoxanthin", "Echinenone"))+
  scale_color_manual(values = c("Grey", "Red"), labels = c("Black Vulture", "Turkey Vulture"))+
  coord_polar("y")+
  theme_void()+
  theme(legend.position = "left")+
  labs(fill='Carotenoid', col='Species')+
  ggtitle ("Species")
species.donut


##Age donut plot
#Melt data frame into lengthwise order
datum2.molten.age <- melt(datum2[,c(4, 43,44,45,46,47,48)], id = c("Age.Broad"))
#Replace NA's with 0 since only trace or undetectable concentrations
datum2.molten.age[is.na(datum2.molten.age)] = 0

datum2.molten.age %>%
  group_by(Age.Broad, variable) %>%
  summarize(mean.conc = mean(value)) -> datum2.molten.age2

#Round for plotting
datum2.molten.age2$mean.conc=round(datum2.molten.age2$mean.conc, digits=0)

#Make plot
age.donut <- ggplot(datum2.molten.age2, aes(x = Age.Broad, y = mean.conc, fill = variable, color=Age.Broad)) +
  geom_col(size=1.1)+ 
  scale_x_discrete(limits = c(" ", "ASY","HY/SY"))+
  geom_text(aes(label = mean.conc),
            position = position_stack(vjust = 0.5), col="black", size=4, vjust=0) +
  scale_fill_manual(values = c("Grey96", "Grey85", "Grey75", "Grey65", "Grey55", "Grey45"), labels = c("Zeaxanthin", "B-carotene", "Lutein", "A-crytpoxanthin", "B-cryptoxanthin", "Echinenone"))+
  scale_color_manual(values = c("Blue", "Skyblue"), labels = c("After second year", "Hatch year/Second year"))+
  coord_polar("y")+
  theme_void()+
  theme(legend.position = "left")+
  labs(fill='Carotenoid', col='Age group')+
  ggtitle ("Age Group")
age.donut

##Sex donut plot
#Melt data frame into lengthwise order
datum2.molten.sex <- melt(datum2[,c(5, 43,44,45,46,47,48)], id = c("Sex"))
#Replace NA's with 0 since only trace or undetectable concentrations
datum2.molten.sex[is.na(datum2.molten.sex)] = 0

datum2.molten.sex %>%
  group_by(Sex, variable) %>%
  summarize(mean.conc = mean(value)) -> datum2.molten.sex2

#Round for plotting
datum2.molten.sex2$mean.conc=round(datum2.molten.sex2$mean.conc, digits=0)

#Make plot
sex.donut <- ggplot(datum2.molten.sex2, aes(x = Sex, y = mean.conc, fill = variable, color=Sex)) +
  geom_col(size=1.1)+ 
  scale_x_discrete(limits = c(" ", "Female","Male"))+
  geom_text(aes(label = mean.conc),
            position = position_stack(vjust = 0.5), col="black", size=4, vjust=0) +
  scale_fill_manual(values = c("Grey96", "Grey85", "Grey75", "Grey65", "Grey55", "Grey45"), labels = c("Zeaxanthin", "B-carotene", "Lutein", "A-crytpoxanthin", "B-cryptoxanthin", "Echinenone"))+
  scale_color_manual(values = c("Darkgreen", "Chartreuse"), labels = c("Female", "Male"))+
  coord_polar("y")+
  theme_void()+
  theme(legend.position = "left")+
  labs(fill='Carotenoid', col='Sex')+
  ggtitle ("Sex")
sex.donut

##Tissue donut plot
#Melt data frame into lengthwise order
datum2.molten.tissue <- melt(datum2[,c(8, 43,44,45,46,47,48)], id = c("Tissue"))
#Replace NA's with 0 since only trace or undetectable concentrations
datum2.molten.tissue[is.na(datum2.molten.tissue)] = 0

datum2.molten.tissue %>%
  group_by(Tissue, variable) %>%
  summarize(mean.conc = mean(value)) -> datum2.molten.tissue2

#Round for plotting
datum2.molten.tissue2$mean.conc=round(datum2.molten.tissue2$mean.conc, digits=0)

#Make plot
tissue.donut <- ggplot(datum2.molten.tissue2, aes(x = Tissue, y = mean.conc, fill = variable, color=Tissue)) +
  geom_col(size=1.1)+ 
  scale_x_discrete(limits = c(" ", "Blood","Skin"))+
  geom_text(aes(label = mean.conc),
            position = position_stack(vjust = 0.5), col="black", size=4, vjust=0) +
  scale_fill_manual(values = c("Grey96", "Grey85", "Grey75", "Grey65", "Grey55", "Grey45"), labels = c("Zeaxanthin", "B-carotene", "Lutein", "A-crytpoxanthin", "B-cryptoxanthin", "Echinenone"))+
  scale_color_manual(values = c("Red", "Darkgrey"), labels = c("Blood", "Skin"))+
  coord_polar("y")+
  theme_void()+
  theme(legend.position = "left")+
  labs(fill='Carotenoid', col='Tissue Type')+
  ggtitle ("Tissue Type")
tissue.donut

donut.figure <- ggarrange(species.donut, age.donut, sex.donut, tissue.donut, ncol=2, nrow=2)

jpeg(file = "overall effect of species tissue age and sex on proportion of each carotenoids.jpg", 
     units = "in", width = 9, height = 9, res = 300)
donut.figure
dev.off()


##Donuts using just turkey vultures for tissue, sex, and age donuts

#I dont want extra clutter with new data frames so just rewrite the old ones
#If you want the previous donut plot, just rerun the code blocks above.

##Species donut plot
#Melt data frame into lengthwise order
datum2.molten.species <- melt(datum2[,c(10, 43,44,45,46,47,48)], id = c("Species"))
#Replace NA's with 0 since only trace or undetectable concentrations
datum2.molten.species[is.na(datum2.molten.species)] = 0

datum2.molten.species %>%
  group_by(Species, variable) %>%
  summarize(mean.conc = mean(value)) -> datum2.molten.species2

#Round for plotting
datum2.molten.species2$mean.conc=round(datum2.molten.species2$mean.conc, digits=0)

#Make plot
species.donut <- ggplot(datum2.molten.species2, aes(x = Species, y = mean.conc, fill = variable, color=Species)) +
  geom_col(size=1.1)+ 
  scale_x_discrete(limits = c(" ", "Black","Turkey"))+
  geom_text(aes(label = mean.conc),
            position = position_stack(vjust = 0.5), col="black", size=4, vjust=0) +
  scale_fill_manual(values = c("Grey96", "Grey85", "Grey75", "Grey65", "Grey55", "Grey45"), labels = c("Zeaxanthin", "B-carotene", "Lutein", "A-crytpoxanthin", "B-cryptoxanthin", "Echinenone"))+
  scale_color_manual(values = c("Grey", "Red"), labels = c("Black Vulture", "Turkey Vulture"))+
  coord_polar("y")+
  theme_void()+
  theme(legend.position = "left")+
  labs(fill='Carotenoid', col='Species')+
  ggtitle ("Species")
species.donut


##Age donut plot
#Melt data frame into lengthwise order
datum2.molten.age <- melt(subset(datum2, Species == "Turkey")[,c(4, 43,44,45,46,47,48)], id = c("Age.Broad"))
#Replace NA's with 0 since only trace or undetectable concentrations
datum2.molten.age[is.na(datum2.molten.age)] = 0

datum2.molten.age %>%
  group_by(Age.Broad, variable) %>%
  summarize(mean.conc = mean(value)) -> datum2.molten.age2

#Round for plotting
datum2.molten.age2$mean.conc=round(datum2.molten.age2$mean.conc, digits=0)

#Make plot
age.donut <- ggplot(datum2.molten.age2, aes(x = Age.Broad, y = mean.conc, fill = variable, color=Age.Broad)) +
  geom_col(size=1.1)+ 
  scale_x_discrete(limits = c(" ", "ASY","HY/SY"))+
  geom_text(aes(label = mean.conc),
            position = position_stack(vjust = 0.5), col="black", size=4, vjust=0) +
  scale_fill_manual(values = c("Grey96", "Grey85", "Grey75", "Grey65", "Grey55", "Grey45"), labels = c("Zeaxanthin", "B-carotene", "Lutein", "A-crytpoxanthin", "B-cryptoxanthin", "Echinenone"))+
  scale_color_manual(values = c("Blue", "Skyblue"), labels = c("After second year", "Hatch year/Second year"))+
  coord_polar("y")+
  theme_void()+
  theme(legend.position = "left")+
  labs(fill='Carotenoid', col='Age group')+
  ggtitle ("Age Group (Turkey Vultures Only)")+
  theme(plot.title = element_text(hjust = 0)) 
age.donut

##Sex donut plot
#Melt data frame into lengthwise order
datum2.molten.sex <- melt(subset(datum2, Species == "Turkey")[,c(5, 43,44,45,46,47,48)], id = c("Sex"))
#Replace NA's with 0 since only trace or undetectable concentrations
datum2.molten.sex[is.na(datum2.molten.sex)] = 0

datum2.molten.sex %>%
  group_by(Sex, variable) %>%
  summarize(mean.conc = mean(value)) -> datum2.molten.sex2

#Round for plotting
datum2.molten.sex2$mean.conc=round(datum2.molten.sex2$mean.conc, digits=0)

#Make plot
sex.donut <- ggplot(datum2.molten.sex2, aes(x = Sex, y = mean.conc, fill = variable, color=Sex)) +
  geom_col(size=1.1)+ 
  scale_x_discrete(limits = c(" ", "Female","Male"))+
  geom_text(aes(label = mean.conc),
            position = position_stack(vjust = 0.5), col="black", size=4, vjust=0) +
  scale_fill_manual(values = c("Grey96", "Grey85", "Grey75", "Grey65", "Grey55", "Grey45"), labels = c("Zeaxanthin", "B-carotene", "Lutein", "A-crytpoxanthin", "B-cryptoxanthin", "Echinenone"))+
  scale_color_manual(values = c("Darkgreen", "Chartreuse"), labels = c("Female", "Male"))+
  coord_polar("y")+
  theme_void()+
  theme(legend.position = "left")+
  labs(fill='Carotenoid', col='Sex')+
  ggtitle ("Sex (Turkey Vultures Only)")+
  theme(plot.title = element_text(hjust = 0)) 
sex.donut

##Tissue donut plot
#Melt data frame into lengthwise order
datum2.molten.tissue <- melt(subset(datum2, Species == "Turkey")[,c(8, 43,44,45,46,47,48)], id = c("Tissue"))
#Replace NA's with 0 since only trace or undetectable concentrations
datum2.molten.tissue[is.na(datum2.molten.tissue)] = 0

datum2.molten.tissue %>%
  group_by(Tissue, variable) %>%
  summarize(mean.conc = mean(value)) -> datum2.molten.tissue2

#Round for plotting
datum2.molten.tissue2$mean.conc=round(datum2.molten.tissue2$mean.conc, digits=0)

#Make plot
tissue.donut <- ggplot(datum2.molten.tissue2, aes(x = Tissue, y = mean.conc, fill = variable, color=Tissue)) +
  geom_col(size=1.1)+ 
  scale_x_discrete(limits = c(" ", "Blood","Skin"))+
  geom_text(aes(label = mean.conc),
            position = position_stack(vjust = 0.5), col="black", size=4, vjust=0) +
  scale_fill_manual(values = c("Grey96", "Grey85", "Grey75", "Grey65", "Grey55", "Grey45"), labels = c("Zeaxanthin", "B-carotene", "Lutein", "A-crytpoxanthin", "B-cryptoxanthin", "Echinenone"))+
  scale_color_manual(values = c("Red", "Darkgrey"), labels = c("Blood", "Skin"))+
  coord_polar("y")+
  theme_void()+
  theme(legend.position = "left")+
  labs(fill='Carotenoid', col='Tissue Type')+
  ggtitle ("Tissue Type (Turkey Vultures Only)")+
  theme(plot.title = element_text(hjust = 0)) 
tissue.donut

donut.figure2 <- ggarrange(species.donut, age.donut, sex.donut, tissue.donut, ncol=2, nrow=2)

jpeg(file = "overall effect on proportions by species tissue age and sex JUST TURKEY.jpg", 
     units = "in", width = 9, height = 9, res = 300)
donut.figure2
dev.off()





####Make animated concentration comparison across all individuals####

#Melt data frame into lengthwise order
datum2.molten.anim <- melt(datum2[,c(4,5,7,8,10, 25,26,27,28,29,30)], id = c("Sample.ID", "Tissue", "Age.Broad", "Species", "Sex"))

#Recode carotenoid names for plotting
datum2.molten.anim %>% 
  mutate(variable = fct_recode(variable,  "Beta-carotene" = "Bcarotene.st", "Zeaxanthin" = "Zeaxanthin.st",
                               "Lutein" = "Lutein.st" , "Alpha-cryptoxanthin" = "Acrypto.st",
                               "Beta-cryptoxanthin" = "Bcrypto.st", "Echinenone" = "Echinenone.st")) -> datum2.molten.anim

#combine tissue and species categories into label for plotting
datum2.molten.anim$group = paste(datum2.molten.anim$Sample.ID, datum2.molten.anim$Species, datum2.molten.anim$Age.Broad,
                                 datum2.molten.anim$Sex)

#Turn NAs to zeros for smooth transitions
datum2.molten.anim[is.na(datum2.molten.anim)] = 0

#By tissue
carotenoid.anim <- ggplot(datum2.molten.anim, aes(x=variable, y=value, fill=variable)) + 
  facet_wrap(~Tissue,  ncol=2, strip.position = "bottom")+
  geom_bar(stat='identity') +
  #geom_text(aes(x=6, y=1.3 , label = as.factor(group)), alpha = 0.9,  col = "black", size = 6) +
  scale_fill_brewer(palette="Dark2")+
  xlab("")+ ylab("Carotenoids (micrograms per g tissue)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust =1),
        legend.position = "")+
  # gganimate specific bits:
  transition_states(group, transition_length = 1, state_length = 1) +
  labs(title = '{closest_state}')+
  ease_aes('linear')

anim <- animate(carotenoid.anim, nframes = 1000, fps = 20, height = 461, width = 644)

anim_save("carotenoids across individuals tissues side by side.gif", anim)

#By individual
carotenoid.anim <- ggplot(datum2.molten.anim, aes(x=variable, y=value, fill=variable)) + 
  facet_wrap(~group,  ncol=6, nrow=2, strip.position = "bottom")+
  geom_bar(stat='identity') +
  #geom_text(aes(x=6, y=1.3 , label = as.factor(group)), alpha = 0.9,  col = "black", size = 6) +
  scale_fill_brewer(palette="Dark2")+
  xlab("")+ ylab("Carotenoids (micrograms per g tissue)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16), axis.text.x = element_text(angle = 65, hjust =1),
        legend.position = "")+
  # gganimate specific bits:
  transition_states(Tissue, transition_length = 1, state_length = 1) +
  labs(title = '{closest_state}')+
  ease_aes('linear')

anim <- animate(carotenoid.anim, nframes = 200, fps = 10, height = 650, width = 1000)

anim_save("carotenoids across individuals tissues transition.gif", anim)

