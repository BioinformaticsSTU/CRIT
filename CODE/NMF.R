setwd("D:/project/github_CRIT/DATA")

#library packages
library(NMF)
library(ggplot2)
library(ComplexHeatmap)

#read features matrix
MATRIX <- read.table(paste("./output/", "MATRIX.txt", sep=""), sep="\t",header=T, row.names=1)

#set run number
runnum <- c(10,30,50,80,100)
ranks  <- 2:100

list_dis_co <- lapply(runnum, function(num){

  nmf.input <- t(MATRIX)

  estim <- lapply(ranks, function(r){

	fit <- nmf(nmf.input, rank = r, nrun = num)

      group <- predict(fit)
      names(group) <- rownames(MATRIX)
      write.table(group,paste("./output/","NMFgroup_",r,"cluster",num,"run",".txt",sep=""),row.names=T,quote=F,col.names=F,sep="\t")
			
      consensus_matrix <- consensus(fit)
      write.table(consensus_matrix,paste("./output/","consensusmatrix_",r,"cluster",num,"run",".txt",sep=""),row.names=T,quote=F,col.names=T,sep="\t")
      pdf(paste("./output/consensusmatrix_",r,"cluster",num,"run",".pdf",sep=""),width=10,height=10)
	pheatmap(consensus_matrix,show_rownames=F,show_colnames=F)
      dev.off()		    

      print(paste("cluster",r,"run",num,dispersion(fit),sep="_"))
      print(paste("cluster",r,"run",num,cophcor(fit),sep="_"))
      list( dispersion = dispersion(fit),cophenetic = cophcor(fit))

  })


index_dispersion <- sapply(estim, '[[', 'dispersion')
index_cophenetic <- sapply(estim, '[[', 'cophenetic')

list(index_dispersion,index_cophenetic)

})



MATRIX_dispersion <- sapply(list_dis_co, '[[', 1)
MATRIX_cophenetic <- sapply(list_dis_co, '[[', 2)

colnames(MATRIX_dispersion) <- paste("runnum_",runnum,sep="")
rownames(MATRIX_dispersion) <- paste("rank_",ranks,sep="")
write.table(MATRIX_dispersion , paste("./output/", "MATRIX_dispersion.txt", sep=""), sep="\t", quote=F)	

colnames(MATRIX_cophenetic) <- paste("runnum_",runnum,sep="")
rownames(MATRIX_cophenetic) <- paste("rank_",ranks,sep="")
write.table(MATRIX_cophenetic, paste("./output/", "MATRIX_cophenetic.txt", sep=""), sep="\t", quote=F)	

mean_dispersionvalue<-rowMeans(MATRIX_dispersion)
max_dispersionvalue<-apply (MATRIX_dispersion,1,max)
min_dispersionvalue<-apply(MATRIX_dispersion,1,min)

mean_copheneticvalue<-rowMeans(MATRIX_cophenetic)
max_copheneticvalue<-apply (MATRIX_cophenetic,1,max)
min_copheneticvalue<-apply(MATRIX_cophenetic,1,min)

MATRIX_dispersionplot<-data.frame(cbind(mean_dispersionvalue,ranks))
colnames(MATRIX_dispersionplot)<-c("Dispersion_coefficient","number_of_clusters")

MATRIX_copheneticplot<-data.frame(cbind(mean_copheneticvalue,ranks))
colnames(MATRIX_copheneticplot)<-c("Cophenetic_coefficient","number_of_clusters")

pdf("./output/dis_co_coefficient.pdf",width=30,height=15)


ggplot(MATRIX_copheneticplot, aes(x=number_of_clusters, y=Cophenetic_coefficient)) + 
  geom_line(color="#67A2D3",size=2) +
  geom_errorbar(aes(ymin=min_copheneticvalue, ymax=max_copheneticvalue), width=.5, position=position_dodge(0.05),size=1)+
  geom_point(data = MATRIX_copheneticplot,aes(x=number_of_clusters, y=Cophenetic_coefficient),size=3,shape = 21, colour = "black", fill = "#67A2D3",stroke = 2)+
    theme_bw()+
      theme(axis.title.x=element_text(face="italic",size=16),   ##name of x.axis 
            axis.text.x=element_text(face="bold",size=16),
		axis.title.y=element_text(face="italic",size=16), 
	      axis.text.y=element_text(face="bold",size=16) ) + 
  geom_line(data = MATRIX_dispersionplot,  aes(x=number_of_clusters, y=Dispersion_coefficient),color = "#AB2C24",size=2)+
  geom_errorbar(aes(ymin=min_dispersionvalue, ymax=max_dispersionvalue), width=.5, position=position_dodge(0.05),size=1)+
  geom_point(data = MATRIX_dispersionplot,  aes(x=number_of_clusters, y=Dispersion_coefficient),size=3,shape = 21, colour = "black", fill = "#AB2C24",stroke = 2)


dev.off()
