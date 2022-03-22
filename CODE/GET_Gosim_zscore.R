#library packages
library(org.Hs.eg.db)
library(clusterProfiler)
library(GOSemSim)
library(GOSim)

#get all ENSG info
index <- keys(org.Hs.eg.db,keytype = "ENSEMBL")

#convert ENSG to symbol and ENTREZID
index_list <- bitr(index,fromType = 'ENSEMBL',
                   toType = c('ENTREZID','SYMBOL'),
                   OrgDb='org.Hs.eg.db',drop = F) 

#get unique ENTREZID list
ID <- as.matrix(as.numeric(unique(index_list$ENTREZID),ncol=1))

#set database
hsGO <- godata('org.Hs.eg.db', ont="BP")

cal_go <- function(id){
  tryCatch({                                                                             #skip error
    if(require(annotate)){
	  setOntology("BP")
	  GO_list <- getGOInfo(id)
    }
	
    GO_names <- as.character(unlist(GO_list[1,]))
    GO_sim   <- mgoSim(GO_names,"GO:0000398" ,semData=hsGO , measure="Lin", combine="BMA")
    return(GO_sim)
  },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
}

ID_sim <- apply(ID,1,cal_go)

sim_matrix <- cbind(ID,data.frame(matrix(t(sapply(ID_sim,c)),ncol=1)))
colnames(sim_matrix) <- c("ENTREZID","GOsim")

del_loc <- union(union(which(sim_matrix$GOsim == "NULL"),which(sim_matrix$GOsim == "")),which(is.na(sim_matrix$GOsim)))
loc     <- setdiff(1:dim(sim_matrix)[1], del_loc) 


#sim_matrix$GOsim == "NULL" : #Error in getGOInfo(3) : No GO information available for these genes!
#sim_matrix$GOsim == ""     : #named list()
#is.na(sim_matrix$GOsim)    : #have GO terms but no sim value


sim <- as.matrix(unlist(sim_matrix$GOsim[loc]), ncol=1)
#output
#write.table(cbind(sim_matrix$ENTREZID[loc],sim),paste("./output/", "GO_sim.txt", sep=""), sep="\t", quote=FALSE,col.names=F,row.names=F)

sim_std     <- sd(sim)*sqrt((length(sim)-1)/(length(sim)))
sim_mean    <- mean(sim)
sim_Z_Score <- (sim-sim_mean)/sim_std

sim_Z_Score_min <- min(sim_Z_Score,na.rm = T)

sim_Z_Score_matrix <- as.matrix(sim_Z_Score, ncol=1)       %>% 
                      cbind(sim_matrix$ENTREZID[loc],.)    %>% 
                      rbind(.,cbind(sim_matrix$ENTREZID[del_loc],rep(sim_Z_Score_min,times=length(del_loc))))

#output
write.table(sim_Z_Score_matrix,paste("./output/", "GO_sim_Z_Score.txt", sep=""), sep="\t", quote=FALSE,col.names=F,row.names=F)
	   






