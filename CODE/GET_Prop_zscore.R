
cal_zscore_p <- function(path){
  
  #library packages
  library('data.table')
  library('magrittr')

  #TCGA input data
  FILE_mutation <- "data_mutations_extended.txt"
  FILE_CNV      <- "data_CNA.txt"
  FILE_mRNA     <- "data_mRNA_median_Zscores.txt"
  FILE_mRNA_2   <- "data_RNA_Seq_v2_mRNA_median_Zscores.txt"

  #choose one mRNA input data 
  files=list.files(path)
    if(!(FILE_mRNA %in% files)){
        file_mRNA=FILE_mRNA_2 
	  }else{
	    file_mRNA=FILE_mRNA 
	}

  #read mutation data 
  mutation_list <- fread(paste(path,FILE_mutation,sep="/"),data.table=F)   %>% 
                   .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")]

  #read mRNA data
  mat_mRNA <- fread(paste(path,file_mRNA,sep="/"), header=TRUE,data.table=F) %>%         #read files
              na.omit(.) %>%                                                             #delete NA rows
              .[,-2] %>%                                                                 #delete id column
              .[!duplicated(.[,'Hugo_Symbol']),]                                         #Keep the first duplicate
  rownames(mat_mRNA) <- mat_mRNA$Hugo_Symbol
  mat_mRNA <- mat_mRNA[,-1]

  #read CNV data
  mat_CNV<-fread(paste(path,FILE_CNV,sep="/"), sep="\t", header=TRUE,data.table=F) %>%   #read files
           na.omit(.) %>%                                                                #delete NA rows
           .[,-2] %>%                                                                    #delete id column
           .[!duplicated(.[,'Hugo_Symbol']),]                                            #Keep the first duplicate
  rownames(mat_CNV) <- mat_CNV$Hugo_Symbol
  mat_CNV<- mat_CNV[,-1]


  #get all sample and all symbol 
  sample_all <- union(union(colnames(mat_mRNA),colnames(mat_CNV)),mutation_list$Tumor_Sample_Barcode)
  symbol_all <- union(union(rownames(mat_mRNA),rownames(mat_CNV)),mutation_list$Hugo_Symbol)


  #get NULL bin matrix
  tmp <- as.data.frame(cbind(symbol_all,matrix(0, ncol=length(sample_all), nrow=length(symbol_all))))
  rownames(tmp) <- symbol_all
  colnames(tmp) <- c("symbol",sample_all)


  #get binary mRNA matrix
  bin_mRNA<- function(mRNA){
    index <- colnames(mat_mRNA)[union(which(mat_mRNA[mRNA[1],] >= 2),which(mat_mRNA[mRNA[1],] <= -2))]
	loc1 <- which(names(mRNA) %in% index )
	
    mRNA[loc1] <- 1                                                                      #convert expression_z_score>=2 or <=-2 to 1
    return(mRNA)
  }

  bin_mRNA_matrix<-t(apply(tmp,1,bin_mRNA))

  #get binary CNV matrix
  bin_CNV<- function(cnv){

    index <- colnames(mat_CNV)[union(which(mat_CNV[cnv[1],] == 2),which(mat_CNV[cnv[1],] == -2))]
	loc1 <- which(names(cnv) %in% index )                                                #convert cnv ==2 or ==-2 to 1         

	cnv[loc1] <- 1	
	return(cnv)
  }


  bin_CNV_matrix<-t(apply(tmp,1,bin_CNV))

  #get binary mutation matrix
  bin_mutation<- function(mutation){
    
	index <- intersect(which(mutation_list$Hugo_Symbol %in% mutation[1]),which(mutation_list$Tumor_Sample_Barcode %in% names(mutation)))    
    loc1 <- which(names(mutation) %in% mutation_list$Tumor_Sample_Barcode[index])
    
    mutation[loc1] <- 1                                                                  #convert any mutation to 1
    return(mutation)
  }

  bin_mutation_matrix<-t(apply(tmp,1,bin_mutation))

  mat_altered <- t(apply(bin_CNV_matrix[,-1],1,as.numeric)+apply(bin_mRNA_matrix[,-1],1,as.numeric)+apply(bin_mutation_matrix[,-1],1,as.numeric))
  colnames(mat_altered) <- sample_all
  
  #output
  write.table(mat_altered, paste("./output/",cancer, "_mat_altered.txt", sep=""), sep="\t", quote=FALSE,col.names=T,row.names=T)

  #get proportion of gene alter
  prop_alter<- function(altered){
    
	index <- length(which(altered >0))/length(altered)	
	return(index)
}

  proportion_alter<-apply(mat_altered,1,prop_alter)

  proportion_scores <- as.matrix(proportion_alter, ncol=1)
  #output
  #write.table(proportion_scores, paste("./output/",cancer, "_proportion_scores.txt", sep=""), sep="\t", quote=FALSE,col.names=F)

  proportion_std     <- sd(proportion_alter)*sqrt((length(proportion_alter)-1)/(length(proportion_alter)))
  proportion_mean    <- mean(proportion_alter)
  proportion_Z_Score <- (proportion_alter-proportion_mean)/proportion_std
  proportion_Z_Score <- as.matrix(proportion_Z_Score, ncol=1)
	
  #output	
  write.table(proportion_Z_Score, paste("./output/",cancer, "_proportion_Z_Score.txt", sep=""), sep="\t", quote=FALSE,col.names=F)
}