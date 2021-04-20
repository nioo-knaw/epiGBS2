rm(list=ls())
library(ggplot2)
library(reshape2)
args<-commandArgs()
data<-read.table(args[1])

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

data$V7<-as.numeric(sub("%","",data$V7))
data$V6<-data$V6/data$V5*100
data$V9<-data$V3-data$V4
data<-data[-c(5)]
colnames(data)<-c("minDepth","clusterID","number of clusters","number of unassembled cluster","percentage multimapped","mapping percentage","average depth","number of assembled clusters")


meltedData<-melt(data,id.vars = c("minDepth","clusterID"))


ggplot(meltedData,aes(x=minDepth,y=value,col=as.factor(clusterID)))+geom_line(size=0.5)+facet_wrap(.~variable,scales = "free")+theme_light()+xlab("minimal depth")+geom_point(size=2)+
  scale_colour_manual("Clustering ID %",values = colorBlindBlack8)+
  ggsave(args[2],height=9,width=9,dpi="retina")

