#### CT2_edgeR.R


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### libraries

library("edgeR")
library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(cowplot)
library(stringr)
library(gtable)
library(pheatmap)
library(RColorBrewer)
require(vegan)
library(pvclust)
library(raster)
library("SuperExactTest")
library(cocor)

print (sessionInfo())

# R version 3.5.1 (2018-07-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.15.5

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] cocor_1.1-3          SuperExactTest_1.0.7 raster_3.0-2         sp_1.3-1             pvclust_2.0-0        vegan_2.5-5          permute_0.9-5        RColorBrewer_1.1-2   pheatmap_1.0.12      gtable_0.3.0         stringr_1.4.0        cowplot_1.0.0       
# [13] lattice_0.20-38      ggplot2_3.2.1        gridExtra_2.3        VennDiagram_1.6.20   futile.logger_1.4.3  edgeR_3.18.1         limma_3.38.3        

# loaded via a namespace (and not attached):
 # [1] Rcpp_1.0.2           pillar_1.4.2         compiler_3.5.1       formatR_1.7          futile.options_1.0.1 tools_3.5.1          nlme_3.1-141         tibble_2.1.3         mgcv_1.8-28          pkgconfig_2.0.2      rlang_0.4.0          Matrix_1.2-17       
# [13] parallel_3.5.1       cluster_2.1.0        withr_2.1.2          dplyr_0.8.3          locfit_1.5-9.1       tidyselect_0.2.5     glue_1.3.1           R6_2.4.0             purrr_0.3.2          lambda.r_1.2.3       magrittr_1.5         codetools_0.2-16    
# [25] splines_3.5.1        MASS_7.3-51.4        scales_1.0.0         assertthat_0.2.1     colorspace_1.4-1     stringi_1.4.3        lazyeval_0.2.2       munsell_0.5.0        crayon_1.3.4        
# > 


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter and normalisation (code)

### takes: DGE structure, cpm cutoff (min cpm), number of samples to apply 
#### way to filter more specifically by sample class, so don't keep sex-limited genes.

filt_and_norm_maj <- function(y,cpm_cut,cut_in_Nsams){
	
	cat("\nNumber number of genes / samples in orig data\n")
	print(dim(y)) ### number of genes / samples
	head(cpm(y)) 
	keep <- 
	  rowSums(cpm(y[,1:3])> cpm_cut) >= cut_in_Nsams & 
      rowSums(cpm(y[,4:6])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,7:9])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,10:12])> cpm_cut) >= cut_in_Nsams      
   

	y <- y[keep,]

	y$samples$lib.size <- colSums(y$counts) # if filter recalc lib size
	y <- calcNormFactors(y)
	cat("\nLib norm factors\n")
	print(y$samples)	
	cat("\nNumber number of genes / samples after filtering\n")
	print(dim(y))
	
	return(y)

}



### MDS function
mds_fact_sex_temp <- function(y){
	mds <- plotMDS(y)
	t1 = as.data.frame(mds$cmdscale.out)
	t1$names <- rownames(t1) 
	a1 <- str_split_fixed(t1$names, "_", 5)
	for_mds <- cbind (t1[,], a1 )
	names(for_mds)[names(for_mds)=="1"] <- "pop"
	names(for_mds)[names(for_mds)=="2"] <- "temp"
	names(for_mds)[names(for_mds)=="3"] <- "sex"
	names(for_mds)[names(for_mds)=="4"] <- "rep"
	for_mds$group <- paste(for_mds$sex, for_mds$temp, sep = "_")
	P0 <- ggplot(for_mds, aes(x=V1, V2, shape = group)) + 
			geom_point(size = 4) +
			scale_shape_manual(values=c(1, 16, 2, 17)) +
			theme_bw() 

	print(for_mds)
	
	return(P0)
}


#################################################################################################
### get DE_gene_table (code)

## returns full table and sig DE genes as a vector
get_DE_genes <- function(fita,FDRa){
	TT1 = topTags(fita, n =3000000000)
	TT2 = TT1$table
	temp_sig <- subset(TT2, TT2$FDR <= FDRa)	
	sig_genes = temp_sig$genes
	sig_logFC = temp_sig$logFC
	
	N_sig_genes <- length(sig_genes)
	cat("Number of sig genes: ", N_sig_genes )
	
	r_list <- list("table" = TT2, "S_gene_list" = sig_genes, "S_logFC_list" = sig_logFC )
	return(r_list)
}


################################################
##### DATA

rawdata <- read.csv("data/Gene_expression/readcounts/CT_H2E.counts.csv", check.names=FALSE, stringsAsFactors=FALSE)

head(rawdata)
colnames(rawdata)
 # [1] "Gene_name"           "KL_Temp19_SexF_Rep1" "KL_Temp19_SexF_Rep2" "KL_Temp19_SexF_Rep3" "KL_Temp19_SexM_Rep1" "KL_Temp19_SexM_Rep2" "KL_Temp19_SexM_Rep3" "KL_Temp6_SexF_Rep1"  "KL_Temp6_SexF_Rep2"  "KL_Temp6_SexF_Rep3"  "KL_Temp6_SexM_Rep1" 
# [12] "KL_Temp6_SexM_Rep2"  "KL_Temp6_SexM_Rep3" 
length(colnames(rawdata))
# 13

##### output

dir.create("CT_analysis")
dir.create("CT_analysis/EdgeR_out")
setwd("CT_analysis/EdgeR_out")

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
## into DGE structure - all species together first

y_all <- DGEList(counts=rawdata[,c(
"KL_Temp6_SexF_Rep1",  "KL_Temp6_SexF_Rep2",  "KL_Temp6_SexF_Rep3",
"KL_Temp19_SexF_Rep1", "KL_Temp19_SexF_Rep2", "KL_Temp19_SexF_Rep3",
"KL_Temp6_SexM_Rep1",  "KL_Temp6_SexM_Rep2",  "KL_Temp6_SexM_Rep3",
"KL_Temp19_SexM_Rep1", "KL_Temp19_SexM_Rep2", "KL_Temp19_SexM_Rep3")], genes=rawdata[,1:1])


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter and normalisation (samples)

y_all_F <- filt_and_norm_maj(y_all,0.5,2) ### 0.5 in at least 2/3 libs for each treatment


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
########get cpm vals

y_all_F_logcpm <- as.data.frame(cpm(y_all_F, log = T))

head(y_all_F_logcpm)

y_all_F_logcpm$KL_Temp6_SexF_meanlogcpm  <- (y_all_F_logcpm$KL_Temp6_SexF_Rep1  +  y_all_F_logcpm$KL_Temp6_SexF_Rep2  +  y_all_F_logcpm$KL_Temp6_SexF_Rep3) / 3
y_all_F_logcpm$KL_Temp19_SexF_meanlogcpm <- (y_all_F_logcpm$KL_Temp19_SexF_Rep1 +  y_all_F_logcpm$KL_Temp19_SexF_Rep2 +  y_all_F_logcpm$KL_Temp19_SexF_Rep3) / 3
y_all_F_logcpm$KL_Temp6_SexM_meanlogcpm  <- (y_all_F_logcpm$KL_Temp6_SexM_Rep1  +  y_all_F_logcpm$KL_Temp6_SexM_Rep2  +  y_all_F_logcpm$KL_Temp6_SexM_Rep3) / 3
y_all_F_logcpm$KL_Temp19_SexM_meanlogcpm <- (y_all_F_logcpm$KL_Temp19_SexM_Rep1 +  y_all_F_logcpm$KL_Temp19_SexM_Rep2 +  y_all_F_logcpm$KL_Temp19_SexM_Rep3) / 3


########## output CPM

write.csv(cbind(y_all_F$genes, y_all_F_logcpm), "logCPM_all.csv", row.names = F)


cor_19 <- cor.test(y_all_F_logcpm$KL_Temp19_SexF_meanlogcpm, y_all_F_logcpm$KL_Temp19_SexM_meanlogcpm, method = 'spearman')
cor_6  <- cor.test(y_all_F_logcpm$KL_Temp6_SexF_meanlogcpm, y_all_F_logcpm$KL_Temp6_SexM_meanlogcpm, method = 'spearman')

### test if cors are different

cocor.indep.groups(cor_19$estimate, cor_6$estimate, length(y_all_F_logcpm[,1]), length(y_all_F_logcpm[,1]), alternative = "two.sided",
test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0,
data.name = NULL, var.labels = NULL, return.htest = FALSE)

# fisher1925: Fisher's z (1925)
  # z = -2.4025, p-value = 0.0163
  # Null hypothesis rejected


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
## plot mds (samples)


MDS_all_F_plot <- mds_fact_sex_temp(y_all_F)
MDS_all_F_plot 



############ DE analysis
####### design matrix
####### design matrix: tissues and wb_match_tiss
sex = factor(c(
"S1F", "S1F", "S1F", "S1F", "S1F", "S1F",
"S0M", "S0M", "S0M", "S0M", "S0M", "S0M"
))

temp = factor(c(
"6", "6", "6", "19", "19", "19",
"6", "6", "6", "19", "19", "19"
))


data.frame(Sample=colnames(y_all),sex,temp)


### design matrix
design_all <- model.matrix(~sex * temp )
rownames(design_all) <- colnames(y_all_F)
design_all

## est dis

y_all_F <- estimateDisp(y_all_F, design_all )

#####################################################################################################################################################################
###### fit glm for interaction  


fit_fullmodel_all_F <- glmQLFit(y_all_F, design_all,robust=TRUE)
colnames(fit_fullmodel_all_F)

# [1] "(Intercept)"  "sexS1F"       "temp6"        "sexS1F:temp6"


### interaction 

fit_fullmodel_all_F_inter_c <- glmQLFTest(fit_fullmodel_all_F, coef=4)
TTT_fullmodel_all_F_inter_c <- get_DE_genes(fit_fullmodel_all_F_inter_c, 0.05)

# Number of sig genes:  64


#####################################################################################################################################################################
###### fit glm for contrasts (code) 

make_Group_matrix <- function(y, fact1,fact2){
	a_samp = data.frame(Sample=colnames(y),fact1,fact2)
	Group <- factor(paste(a_samp$fact1,a_samp$fact2, sep="."))
	print(Group)
	cbind(a_samp,Group=Group)
	cat(Group)
	G_design <- model.matrix(~0+Group)
	colnames(G_design) <- levels(Group)
	print(G_design)
	return(G_design)
}


####### fit glm (samples)

design_all_group <- make_Group_matrix(y_all_F,sex,temp)
fit_groupsmodel_all_F <- glmQLFit(y_all_F, design_all_group ,robust=TRUE)  


#################################################################################################
### make contrasts


all.contrasts <- makeContrasts(
female_CT =	S1F.6 - S1F.19,   ### +ve FC = higher 6C
male_CT	  =	S0M.6 - S0M.19 , 
SB_19     =	S1F.19 - S0M.19 ,  ### +ve FC = higher experssion in female
SB_6	      =	S1F.6 - S0M.6,
levels=design_all_group)

female_CT_contrast  <- glmQLFTest(fit_groupsmodel_all_F, contrast=all.contrasts[,"female_CT"])
male_CT_contrast    <- glmQLFTest(fit_groupsmodel_all_F, contrast=all.contrasts[,"male_CT"])
SB_19_contrast      <- glmQLFTest(fit_groupsmodel_all_F, contrast=all.contrasts[,"SB_19"])
SB_6_contrast       <- glmQLFTest(fit_groupsmodel_all_F, contrast=all.contrasts[,"SB_6"])


TTT_female_CT_contrast <- get_DE_genes(female_CT_contrast, 0.05)
# Number of sig genes:  1062

TTT_male_CT_contrast <- get_DE_genes(male_CT_contrast, 0.05)
# Number of sig genes:  1236

TTT_SB_19_contrast <- get_DE_genes(SB_19_contrast, 0.05)
# Number of sig genes:  6747

TTT_SB_6_contrast <- get_DE_genes(SB_6_contrast, 0.05)
# Number of sig genes:  7005


####### export full tables

write.table(TTT_fullmodel_all_F_inter_c$table,     "TTT_CT_sex_temp_inter.csv",       sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_female_CT_contrast$table,          "TTT_female_CT_contrast.csv",      sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_male_CT_contrast$table,            "TTT_male_CT_contrast.csv",        sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_SB_19_contrast$table,              "TTT_SB_19_contrast.csv",          sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_SB_6_contrast$table,               "TTT_SB_6_contrast.csv",           sep = ',', quote = FALSE, row.names = FALSE)










