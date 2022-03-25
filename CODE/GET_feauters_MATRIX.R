#library packages
library(data.table)


#setwd("D:/project/CRIT_github/DATA")

path_list <- list.files('./TCGA')

RBP_list <- fread("./output/RBP_symbol_id_uniprot.txt")

Nor <- function(x){
	   min_x <- min(x,na.rm=T)
	   max_x <- max(x,na.rm=T)
	   nor_x <- (x-min_x)/(max_x-min_x)
	   return(nor_x)
}


#get Zp and Zs feautures matrix
MATRIX_Zp <-NULL
MATRIX_Zs <-NULL

for(i in 1:length(path_list)){
  
  cancer <- path_list[i]
  
  Zp <- read.table(paste("./output/",cancer, "_proportion_Z_Score.txt", sep=""),  sep="\t", header=F,row.names=1)
  Zs <- read.table(paste("./output/",cancer, "_sur_p_value_Z_Score.txt", sep=""), sep="\t", header=F,row.names=1)
  
  index_Zp <- Nor(Zp[RBP_list$RBP,])
  index_Zs <- Nor(Zs[RBP_list$RBP,])

  MATRIX_Zp <- cbind(MATRIX_Zp,index_Zp)
  colnames(MATRIX_Zp)[i] <- paste(cancer,"_Zp",sep="")
	
  MATRIX_Zs <- cbind(MATRIX_Zs,index_Zs)
  colnames(MATRIX_Zs)[i] <- paste(cancer,"_Zs",sep="")
}

rownames(MATRIX_Zp) <- RBP_list$RBP
rownames(MATRIX_Zs) <- RBP_list$RBP
#output
#write.table(MATRIX_Zp , paste("./output/", "MATRIX_Zp.txt", sep=""), sep="\t", quote=F)
#write.table(MATRIX_Zs , paste("./output/", "MATRIX_Zs.txt", sep=""), sep="\t", quote=F)


#get Zg feautures matrix
Zg <- read.table(paste("./output/","GO_sim_Z_Score.txt", sep=""), header=F,row.names=1)
Zg_value  <- Nor(Zg[as.character(RBP_list$entrezgene_id),])
MATRIX_Zg <- cbind(RBP_list$RBP,Zg_value)
#output
#write.table(MATRIX_Zg , paste("./output/", "MATRIX_Zg.txt", sep=""), sep="\t", quote=F,row.names=F)


#get cor feautures matrix
MATRIX_cor<-NULL
cancertype <- c("BLCA","BRCA","COLO","ESCA","GBM","KDNY","AML","LUNG","OV","PRAD","STAD")

for(i in 1:length(cancertype)){
  cor <- read.table(paste("./output/",cancertype[i],"_symbol_top50_meancor.txt",sep=""))
  
  index_cor <- Nor(cor[RBP_list$RBP,])
  MATRIX_cor <- cbind(MATRIX_cor ,index_cor )

  colnames(MATRIX_cor)[i] <- paste(cancertype[i],"_cor",sep="")	
}
rownames(MATRIX_cor) <- RBP_list$RBP
#output
#write.table(MATRIX_cor, paste("./output/", "MATRIX_cor.txt", sep=""), sep="\t", quote=F)


#get C value feautures matrix	
MATRIX_C<-NULL

C_value <- read.table(paste("./output/","C_value_matrix.txt", sep=""),header=T)
rownames(C_value) <- C_value$RBP

MATRIX_C <- as.matrix(Nor(C_value[RBP_list$RBP,"C_value"]),ncol=1)
rownames(MATRIX_C) <- RBP_list$RBP
colnames(MATRIX_C) <- "c_value"
#output
#write.table(MATRIX_C, paste("./output/", "MATRIX_C.txt", sep=""), sep="\t", quote=F)
	
	
MATRIX<-cbind(MATRIX_Zp, MATRIX_Zs, MATRIX_Zg[,2], MATRIX_cor, MATRIX_C)
colnames(MATRIX)[39] <- "Zg"
MATRIX<-na.omit(MATRIX)

#output 
write.table(MATRIX , paste("./output/", "MATRIX.txt", sep=""), sep="\t", quote=F,row.names=T,col.names=T)