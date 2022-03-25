######################################## SETTINGS AND COMMANDS ########################################
## 1. collection of features

# prepare: set directory and files
# Working directory:

DATApath <- "./DATA"
setwd(DATApath)

# 1.1 calculate feature of Zp and Zs for TCGA data

source("../CODE/GET_Prop_zscore.R")			# provided.     
source("../CODE/GET_Sur_zscore.R")			# provided.

path_list <- list.files('./input/TCGA/')

for(i in 1:length(path_list)){

  path <- paste('./input/TCGA/',path_list[i],sep="")
  cancer <- path_list[i]

  cal_zscore_p(path)
  cal_zscore_s(path)
}

# output File
#./output/XX_mat_altered.txt 			    #alternation matrix of genes in XX cancer type of TCGA  
#./output/XX_proportion_Z_Score.txt         #Zp feature of genes in XX cancer type of TCGA
#./output/XX _sur_p_value_Z_Score.txt       #Zs feature of genes in XX cancer type of TCGA

# 1.2 caluate Zg for all genes in ENSEMBL database

source("../CODE/GET_Gosim_zscore.R")	    # provided.

# output File
#./output/GO_sim_Z_Score.txt		#Zs feature of genes in ENSEMBL database.

# 1.3 calculate feature of ρm for MiOncoCirc data

# setting cancertype of MiOncoCirc:

cancertype<-c("BLCA","BRCA","COLO","ESCA","GBM","KDNY","AML","LUNG","OV","PRAD","STAD")

source("../CODE/GET_mRNAcirc_top50cor.R")	# provided.

for(i in 1:length(cancertype)){
  cal_mRNA_circ_cor(cancertype[i])
}

# output File
#./output/XX_symbol_top50_meancor.txt		#ρm feature of genes in XX cancer type of MiOncoCirc

# 1.4 calculate feature of C value for RNAct data

# setting candidate regulator list:

load("../CODE/RBP_list.RData")			    # provided RBP_list
source("../CODE/GET_C_value.R")			    # provided.

cal_C_value(RBP_list)

# output File
#./output/C_value_matrix.txt		        #C value feature of genes in candidate regulators.

# 1.5 integrate features of genes

source("../GET_feauters_MATRIX.R")			# provided.

# output File
#./output/MATRIX.txt                    	#all features of candidate regulators.

#------------------------------------------------------------------------------------------------------
## 2. NMF algorithm

#set run number and ranks number

runnum <- c(10,30,50,80,100)                #set run numbers 
ranks  <- 2:50                              #set ranks

source("../CODE/NMF.R")			    # provided.

# output File
#./output/NMFgroup_rclusternrun.txt         #RBP group in r cluster and n run
#./output/consensusmatrix__rclusternrun.txt #consensus matrix of r cluster and n run
#./output/consensusmatrix__rclusternrun.pdf #heatmap of consensus matrix of r cluster and n run
#./output/MATRIX_dispersion.txt             #dispersion matrix of RBP group
#./output/MATRIX_cophenetic.txt             #cophenetic matrix of RBP group
#./output/dis_co_coefficient.pdf            #plot of cophenetic and dispersion matrix

#------------------------------------------------------------------------------------------------------
## 3. Validate candidates

#set run number and ranks number 

select_cluster_num <- 36                    #set ranks
n_run <- 50                                 #set run numbers

source("../CODE/RBP_analysis.R")			    # provided.

# output File
#./output/NMF_resultrclustern_run.txt       #p value of each cluster of r cluster and n run
#./output/physical_result.txt
#./output/coexp_result.txt
#./output/tf_result.txt

#######################################################################################################