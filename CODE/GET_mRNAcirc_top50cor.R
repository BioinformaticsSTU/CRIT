
cal_mRNA_circ_cor <- function(cancertype){

  #load packages
  library(xlsx)
  library(magrittr)
  library(data.table)
  library(dplyr)

  #set working path
  #setwd("./MiOncoCirc")

  #read mRNA expression file
  FPKM <- fread("./input/MiOncoCirc/fpkm_matrix.csv", sep=",", header=TRUE,data.table = F,check.names=FALSE) %>%
          .[,-1]
  #53986row,878col,1col:ENSG,878col:gene_symbol

  #delete NULL expression row
  valuesum <- rowSums(FPKM[,-c(1,878)])
  FPKM <- FPKM[-which(valuesum == 0),]

  #Keep the first duplicate
  FPKM <- FPKM[!duplicated(FPKM[,878]),]
  rownames(FPKM)<-FPKM$V2

  #read circRNA expression file
  reads <- fread("./input/MiOncoCirc/v0.1.release.txt", header=TRUE,data.table = F,check.names=FALSE) %>% 
           .[,c("reads","symbol","sample")]
  #10492676row,7col,1~3col:chr,4col:reads,5col:symbol,6col:sample,7col:release

  #read matadata of sample
  v00meta <- read.xlsx("./input/MiOncoCirc/meta_update.xlsx", sheetIndex ="v0.0", check.names = FALSE, header=TRUE) %>% 
             .[,c("ID","Analysis Cohort")]
  v01meta <- read.xlsx("./input/MiOncoCirc/meta_update.xlsx", sheetIndex ="v0.1_additional", check.names = FALSE, header=TRUE) %>% 
             .[,c("ID","Analysis Cohort")]
  meta    <- rbind(v00meta,v01meta)


  #spearson correlation of mRNA and circRNA per cancer type
  for(i in 1:length(cancertype)){
    type_sample<-meta[which(meta$"Analysis Cohort" == cancertype[i]),1]
    inter_sample<-intersect(intersect(type_sample,names(FPKM)),reads$sample)

    mRNAexp<-FPKM[,inter_sample] %>% .[,order(colnames(.))]

    circ_list<-reads[which(reads$sample %in% inter_sample),]
    circ_symbol<-unique(circ_list$symbol)
    circ_sample<-unique(circ_list$sample)
    circexp<-matrix(NA, ncol=length(circ_sample), nrow=length(circ_symbol))
    rownames(circexp)<-circ_symbol
    colnames(circexp)<-circ_sample
      for(j in 1:length(circ_sample)){
         for(n in 1:length(circ_symbol)){
	        #for duplicate circRNA expression,get sum expression
            circexp[n,j]<-sum(circ_list[intersect(which(circ_list$sample == circ_sample[j]),which(circ_list$symbol == circ_symbol[n])),1])
         }
      }
	
    #delete gene of 0 expression > 1/3 samples 
    symbol_count0<-(apply(circexp, 1, function(x){length(which(x == 0))}))
    circexp<-circexp[-which(symbol_count0>length(circ_sample)/3),] %>% .[,order(colnames(.))]
   
    #spearman correlation of mRNAexp and circRNA
    mat_mRNAcirc_cor<-cor(t(mRNAexp),t(circexp), method =  "spearman")#row:gene,col:circ
   
    #checkpoint
    #print(dim(mat_mRNAcirc_cor))
   
    #output files
    #write.table(mat_mRNAcirc_cor,paste("./output/",cancertype[i],"_mat_mRNAcirc_cor.txt",sep=""),sep="\t",quote=F)
   
    symbol_top50_meancor<-apply(mat_mRNAcirc_cor, 1, function(x){mean(abs(as.vector(x[order(desc(abs(x)))][1:50])))} )
   
    #output files
    write.table(symbol_top50_meancor,paste("./output/",cancertype[i],"_symbol_top50_meancor.txt",sep=""),sep="\t",quote=F,col.names=F)
  }  
}
