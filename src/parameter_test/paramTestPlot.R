rm(list=ls())
library(ggplot2)
library(reshape2)
args<-commandArgs(trailingOnly = T)
dataDenovo<-read.table(args[1])

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

dataDenovo$V7<-as.numeric(sub("%","",dataDenovo$V7))
dataDenovo$V6<-dataDenovo$V6/dataDenovo$V5*100
dataDenovo$V9<-dataDenovo$V3-dataDenovo$V4
dataDenovo<-dataDenovo[-c(5)]
colnames(dataDenovo)<-c("minDepth","clusterID","number of clusters","number of unassembled cluster","percentage multimapped","mapping percentage","average depth","number of assembled clusters")


melteddataDenovo<-melt(dataDenovo,id.vars = c("minDepth","clusterID"))


ggplot(melteddataDenovo,aes(x=minDepth,y=value,col=as.factor(clusterID)))+geom_line(size=0.5)+facet_wrap(.~variable,scales = "free")+theme_light()+xlab("minimal depth")+geom_point(size=2)+
  scale_colour_manual("Clustering ID %",values = colorBlindBlack8)+
  ggsave(args[2],height=9,width=9,dpi="retina")

