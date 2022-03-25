
cal_zscore_s <- function(path){
  
  #library packages
  library(survival)
  
  #TCGA input data
  FILE_clinical <- "data_clinical_patient.txt"
  
  #read clinical data 
  mat_cli <- fread(paste(path,FILE_clinical,sep="/"),data.table=F,skip=4) %>% 
             .[,c("PATIENT_ID","OS_MONTHS","OS_STATUS")] 
  mat_cli$OS_MONTHS <- as.numeric(mat_cli$OS_MONTHS) 
  mat_cli <- na.omit(mat_cli)

  rownames(mat_cli) <- mat_cli$PATIENT_ID 
   
  #TCGA input altered data   
  file_altered <- paste("./output/",cancer, "_mat_altered.txt", sep="")

  mat_altered <- read.table(file_altered, header=T,sep = "\t",check.names=F)
  colnames(mat_altered) <- substr(colnames(mat_altered),1,12)                            #sample name convert

  #get intersect sample 
  sample_intersect <- intersect(colnames(mat_altered),mat_cli$PATIENT_ID)

  mat_cli <- mat_cli[sample_intersect,]
  mat_altered <- mat_altered[,sample_intersect]

  get_sur_p <- function(alter){

        alt <- as.numeric(as.logical(alter))
            if(sum(alt) == 0){
                p_value <- NA
            }else{
                index1 <- cbind(mat_cli,alt)
                
                loc1 <- which(index1$OS_STATUS == "DECEASED")
                index1[loc1,"OS_STATUS"]  <- 1
                index1[setdiff(1:dim(index1)[1],loc1),"OS_STATUS"] <- 0
		
                surv_diff <- survdiff(Surv(OS_MONTHS, as.numeric(OS_STATUS)) ~ alt, data =index1)
                p_value <- 1-pchisq(surv_diff$chisq, length(surv_diff$n)-1)
            }
    }    

  sur_p_value <- apply(mat_altered,1,get_sur_p)


  sur_p_value_matrix <- as.matrix(sur_p_value, ncol=1)
  #output
  #write.table(sur_p_value_matrix, paste("./output/",cancer, "_sur_p_value.txt", sep=""), sep="\t", quote=FALSE,col.names=F)


  sur_p_value_Z_Score <- qnorm(1-sur_p_value)

  del_loc <- union(which(sur_p_value_Z_Score == "-Inf"),which(sur_p_value_Z_Score == Inf))
  loc <- setdiff (1:length(sur_p_value_Z_Score), del_loc)  

  sur_p_value_Z_Score_min <- min(sur_p_value_Z_Score[loc],na.rm = T)
  sur_p_value_Z_Score_max <- max(sur_p_value_Z_Score[loc],na.rm = T)

  sur_p_value_Z_Score[union(which(sur_p_value_Z_Score == -Inf),which(is.na(sur_p_value_Z_Score)))] <- sur_p_value_Z_Score_min
  sur_p_value_Z_Score[which(sur_p_value_Z_Score ==  Inf)] <- sur_p_value_Z_Score_max

  sur_p_value_Z_Score_matrix <- as.matrix(sur_p_value_Z_Score, ncol=1)

  #output	
  write.table(sur_p_value_Z_Score_matrix, paste("./output/",cancer, "_sur_p_value_Z_Score.txt", sep=""), sep="\t", quote=FALSE,col.names=F)
}

