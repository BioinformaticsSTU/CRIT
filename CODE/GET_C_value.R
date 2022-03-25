
cal_C_value <- function(RBP_list){

  #library packages
  library(magrittr)
  library(readr)
  library(stringr)

  #get uniprot_id
  RBP_uniprot_list <- na.omit(unique(RBP_list$uniprotswissprot))

  #get uniprot id and score>=30
  f <- function(x, pos) subset(x,x$uniprot_accession %in%  RBP_uniprot_list & x$catrapid_score_max_normalised >= 30)
  res<-read_delim_chunked(file="./input/RNAct/catrapid_human_basic_normalised.txt",
                          delim="\t",DataFrameCallback$new(f), chunk_size = 10000,col_names = TRUE)
  #output 
  #write_excel_csv(res,"./output/RBP_scores.csv")

  #get enst list
  enst_list<-unique(res$ensembl_transcript_id)

  #get enst str location information
  #read annotation file
  gencode <- fread("./input/RNAct/gencode.v29lift37.annotation.gtf",skip="__auto__",header=F,data.table=F,sep="\t") %>%
                   .[.$V3 == "transcript",]
  index <- str_split_fixed(gencode[,9], "\\s+",n=30)

  #keep enst location 
  trans_pos <- cbind(substr(index[,4], 2,18), gencode[,c(1,4,5,7)])
  colnames(trans_pos) <- c("enst","chr","start","end","strand")

  #get circRNA str location information
  circ_pos <- fread("./input/RNAct/hsa_hg19_circRNA.txt",,header=T,data.table=F,sep="\t") %>%
              .[,c("circRNA ID","# chrom","start","end","strand","best transcript","gene symbol")]

  #get enst list loc
  enst_pos <- cbind(substr(trans_pos$enst,1,15),trans_pos[,2:5]) %>% .[which(.[,1] %in% substr(enst_list,1,15)),]
  colnames(enst_pos)[1]<- "enst"

  #get list of circRNA matched in RNA 
  get_circRNA<-function(x){
               
			   chr <- as.character(x[2])
			   strand <-as.character(x[5])
			   start <- as.numeric(x[3])
			   end <- as.numeric(x[4])
			   
			   index <- subset(circ_pos,circ_pos$'# chrom' == chr & circ_pos$strand == strand )
			   index <- index[which(index$start >= start  & index$end <= end),]
			   circRNA_list <- paste(index$"circRNA ID"[1:length(index$"circRNA ID")],collapse=",")            
  }


  #get circRNA list
  circRNA_list <- apply(enst_pos,1,get_circRNA)

  enst_circ <-cbind(enst_pos,circRNA_list)

  uniprot_list <- as.data.frame(unique(res$uniprot_accession))

  #get number of unique circ list
  count_uniprot_circ <- function(x){

                    index <- subset(res,res$uniprot_accession == x)                
                    index <- substr(unique(index$ensembl_transcript_id),1,15)
                    length(setdiff(unique(unlist(strsplit(subset(enst_circ,enst_circ$enst %in% index)$circRNA_list ,","))),"NA"))
  }

  uniprot_circ_counts<-apply(uniprot_list,1,count_uniprot_circ)

  uniprot_circ <-cbind(uniprot_list,uniprot_circ_counts)
  colnames(uniprot_circ) <- c("uniprotswissprot","C_value")
  C_value_matrix <- merge(RBP_list,uniprot_circ,by="uniprotswissprot",all.x=T)

  #output
  write.table(C_value_matrix,paste("./output/", "C_value_matrix.txt", sep=""), sep="\t", quote=FALSE,col.names=T,row.names=F)
}