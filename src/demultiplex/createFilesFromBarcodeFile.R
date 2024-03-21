rm(list=ls())

args = commandArgs(trailingOnly=TRUE)


barcodes<-read.table(args[1],h=T,sep="\t")
barcodes$run<-sub("_R1.fq.gz","",barcodes$rawR1)

barcodeStacks<-data.frame(sample=c(paste0(barcodes$sample,"-Watson"),paste0(barcodes$sample,"-Crick")),
                          barcode1=c(paste0(barcodes$barcode1,"T"),paste0(barcodes$barcode1,"C")),
                          barcode2=c(paste0(barcodes$barcode2,"C"),paste0(barcodes$barcode2,"T")))
write.table(barcodeStacks,paste0(args[2],"/barcodeStacks.tsv"),row.names = F,col.names = F,quote=F,sep="\t")

for(run in barcodes$run){
  barcodesSub<-barcodes[barcodes$run==run,]
  barcodeStacks<-data.frame(barcode1=paste0(barcodesSub$barcode1,"C"),
                            barcode2=paste0(barcodesSub$barcode2,"C"),
                            sample=barcodesSub$sample)
  write.table(barcodeStacks,paste0(args[2],"/barcodeStacks",run,".tsv"),row.names = F,col.names = F,quote=F,sep="\t")
              
              
}