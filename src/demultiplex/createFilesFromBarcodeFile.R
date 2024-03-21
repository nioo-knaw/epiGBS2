rm(list=ls())

args = commandArgs(trailingOnly=TRUE)


barcodes<-read.table(args[1],h=T,sep="\t")
barcodes$run<-sub("_R1.fq.gz","",barcodes$rawR1)

barcodeStacks<-data.frame(sample=c(paste0(barcodes$Sample,"-Watson"),paste0(barcodes$Sample,"-Crick")),
                          barcode1=c(paste0(barcodes$Barcode_R1,"T"),paste0(barcodes$Barcode_R1,"C")),
                          barcode2=c(paste0(barcodes$Barcode_R2,"C"),paste0(barcodes$Barcode_R2,"T")))
write.table(barcodeStacks,paste0(args[2],"/barcode_stacks.tsv"),row.names = F,col.names = F,quote=F,sep="\t")

for(run in barcodes$run){
  barcodesSub<-barcodes[barcodes$run==run,]
barcodeStacks<-data.frame(sample=c(paste0(barcodes$Sample,"-Watson"),paste0(barcodes$Sample,"-Crick")),
                          barcode1=c(paste0(barcodes$Barcode_R1,"T"),paste0(barcodes$Barcode_R1,"C")),
                          barcode2=c(paste0(barcodes$Barcode_R2,"C"),paste0(barcodes$Barcode_R2,"T")))
  write.table(barcodeStacks,paste0(args[2],"/barcode_stacks",run,".tsv"),row.names = F,col.names = F,quote=F,sep="\t")
              
              
}