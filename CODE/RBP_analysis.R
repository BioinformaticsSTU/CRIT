library(magrittr)

RBP_paper<-read.xlsx("./input/gold_regulators.xlsx",header=F,sheet=2) %>% .$V1
RBP_Li<-read.csv("./input/gold_regulators.xlsx",header=F,sheet=1) %>% .$V1

RBP_verified_all<-union(RBP_paper,RBP_Li)


sapply(n_run,function(run){

      sapply(select_cluster_num,function(num){

	    result <- read.table(paste("./output/","NMFgroup_",select_cluster_num ,"cluster" ,run,"run.txt",sep=""))

	    group <- result$V2
	    names(group) <- result$V1
	    group_names <- names(table(group))

	    MATRIX_pvalue <- t(sapply(group_names,function(g){	
	
          RBP_group <- names(group)[which(group== g)]
          RBP_group_verified <- intersect(RBP_verified_all, RBP_group)		
          RBP_all_verified<-intersect(RBP_verified_all,names(group))
          RBP_allgroup<-names(group)
		
          len_group<-length(RBP_group)
          len_group_verified<-length(RBP_group_verified)
          len_all_verified<-length(RBP_all_verified)
          len_allgroup<-length(RBP_allgroup)
		
	    	mat=matrix(c(len_group_verified, 
		                 len_group-len_group_verified, 
			    		 len_all_verified-len_group_verified, 
		    			 len_allgroup-len_all_verified-len_group+len_group_verified) ,ncol=2, nrow=2)
		    #print(mat)
	    	p_value<-fisher.test(mat)$p.value	
		
	    	index<-c(paste(sort(RBP_group),collapse=","),len_group,paste(sort(RBP_group_verified),collapse=","),len_group_verified,p_value)		
    	}
		))
	
	    rownames(MATRIX_pvalue) <- group_names
	    colnames(MATRIX_pvalue) <-c("Total gene","Num of genes","Included in gold RBPs","Num of golden RBPs","Fisher_p_value")
        #output	
	    write.table(MATRIX_pvalue,paste("./output/","NMF_result",select_cluster_num,"cluster",run,"_run.txt",sep=""),row.names=F,quote=F,col.names=T,sep="\t")
})
})


options(scipen=200)

RBP_list <- read.table("./output/NMF_result36cluster50_run.txt",header=T,sep="\t",check.names=F)   %>% 
            .[which(.$Fisher_p_value == min(.$Fisher_p_value)),"Total gene"]      %>%
	      strsplit(.,split=",") %>% unlist()

RBP_verified<-intersect(RBP_verified_all,RBP_list)
RBP_unknown<-setdiff(RBP_list,RBP_verified)


####
pro_physical<-read.table("./input/Harmonizome/physical/gene_attribute_edges.txt",header=T,skip=1)

verified_physical<-pro_physical[union(which(pro_physical[,1] %in% RBP_verified),which(pro_physical[,4] %in% RBP_verified)),]
#81040

verified_all_physical<-pro_physical[union(which(pro_physical[,1] %in% RBP_verified_all),which(pro_physical[,4] %in% RBP_verified_all)),]
#461608

all_pro<-union(pro_physical[,1],pro_physical[,4])


p_val<-NULL

all_verified_RBPs<-NULL

allnums<-NULL



for(i in 1:length(RBP_unknown)){

  index1<-union(pro_physical[which(pro_physical[,1] %in% RBP_unknown[i]),4],pro_physical[which(pro_physical[,4] %in% RBP_unknown[i]),1])
  #all protein
  index2<-length(index1)
  #number

  index3<-union(verified_all_physical[which(verified_all_physical[,1] %in% RBP_unknown[i]),4],verified_all_physical[which(verified_all_physical[,4] %in% RBP_unknown[i]),1])
  #all protein in gold
  index4<-length(index3)
  #number

  allnums<-c(allnums,index4)
  index10<-paste(sort(index3),collapse=",")
  #gold list 
  all_verified_RBPs<-c(all_verified_RBPs,index10)
  

  pval<-fisher.test(matrix(c(index4,length(RBP_verified_all)-index4,index2-index4,length(all_pro)-index2-length(RBP_verified_all)+index4),nrow=2,ncol=2))$p.value
  p_val<-c(p_val,pval)

}
 

FDR <- p.adjust(p_val, method="fdr")#fdr 

result<-cbind(RBP_unknown,all_verified_RBPs,allnums,p_val,FDR)
colnames(result) <- c("Candidate regulators","Physical interaction gold regulators","Num of gold regulators","Fisher_p_val","FDR")

write.table(result,"./output/physical_result.txt",col.names=T,row.names=F,sep="\t")






pro_coexp<-read.table("./input/Harmonizome/coexp/gene_attribute_edges.txt",header=T,skip=1)

verified_coexp<-pro_coexp[union(which(pro_coexp[,1] %in% RBP_verified),which(pro_coexp[,4] %in% RBP_verified)),]
#548

verified_all_coexp<-pro_coexp[union(which(pro_coexp[,1] %in% RBP_verified_all),which(pro_coexp[,4] %in% RBP_verified_all)),]
#2924

all_pro<-union(pro_coexp[,1],pro_coexp[,4])
#4870

dim(pro_coexp)
#41327


p_val<-NULL

all_verified_RBPs<-NULL

allnums<-NULL



for(i in 1:length(RBP_unknown)){

  index1<-setdiff(union(pro_coexp[which(pro_coexp[,1] ==RBP_unknown[i]),4],pro_coexp[which(pro_coexp[,4] ==RBP_unknown[i]),1]),RBP_unknown[i])
  #all coexp protein
  index2<-length(index1)
  #number

  index3<-union(verified_all_coexp[which(verified_all_coexp[,1] == RBP_unknown[i]),4],verified_all_coexp[which(verified_all_coexp[,4] == RBP_unknown[i]),1])
  #all coexp gold 
  index4<-length(index3)
  #number

  allnums<-c(allnums,index4)
  index10<-paste(sort(index3),collapse=",")
  #coexp list
  all_verified_RBPs<-c(all_verified_RBPs,index10)

  
  pval<-fisher.test(matrix(c(index4,length(RBP_verified_all)-index4,index2-index4,length(all_pro)-index2-length(RBP_verified_all)+index4),nrow=2,ncol=2))$p.value
  p_val<-c(p_val,pval)

}
 

FDR <- p.adjust(p_val, method="fdr")#fdr 

result<-cbind(RBP_unknown,all_verified_RBPs,allnums,p_val,FDR)
colnames(result) <- c("Candidate regulators","coexp gold regulators","Num of gold regulators","Fisher_p_val","FDR")

write.table(result,"./output/coexp_result.txt",col.names=T,row.names=F,sep="\t")




pro_tf<-read.table("./input/Harmonizome/TF/gene_attribute_edges.txt",header=T,skip=1)


#pro_tf, first col is gene, fourth col is tf

verified_tf<-pro_tf[which(pro_tf[,1] %in% RBP_verified),]
#201

verified_all_tf<-pro_tf[which(pro_tf[,1] %in% RBP_verified_all),]
#1788

all_pro<-unique(pro_tf[,1])
#13216
len1<-length(which(RBP_verified_all %in% all_pro))

#110
len2<-length(which(RBP_verified %in% all_pro))
#12


dim(pro_tf)
#100560


p_val<-NULL

all_verified_RBPs<-NULL

allnums<-NULL

cotf<-NULL

for(i in 1:length(RBP_unknown)){

  index1<-setdiff(pro_tf[which(pro_tf[,1] == RBP_unknown[i]),4],RBP_unknown[i])
  #all tf
  index2<-unique(pro_tf[which(pro_tf[,4] %in% index1),1])
  #length(index2)
  #number of gene

  index3<-index2[which(index2 %in% RBP_verified_all)]
  #all gene in gold
  index4<-length(index3)
  #number
  
  allnums<-c(allnums,index4)
  index10<-paste(sort(index3),collapse=",")
  #tf list
  all_verified_RBPs<-c(all_verified_RBPs,index10)

  
  pval<-fisher.test(matrix(c(index4,len1-index4,length(index2)-index4,length(all_pro)-length(index2)-len1+index4),nrow=2,ncol=2))$p.value
  p_val<-c(p_val,pval)

  cotf <- c(cotf,paste(sort(index1),collapse=","))


}
 

FDR <- p.adjust(p_val, method="fdr")#fdr

result<-cbind(RBP_unknown,cotf,all_verified_RBPs,allnums,p_val,FDR)
colnames(result) <- c("Candidate regulators","co tf","Gold regulators","Num of gold regulators","Fisher_p_val","FDR")

write.table(result,"./output/tf_result.txt",col.names=T,row.names=F,sep="\t")

