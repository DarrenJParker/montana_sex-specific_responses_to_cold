#### CT_overlap.R


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


#######################################################
### read in edgeR output

dat_CT_edgeR <- read.csv("CT_analysis/EdgeR_out/TTT_CT2_sig_CTSB_tidied.csv", header = T)
head(dat_CT_edgeR)
all_expressed_genes <- as.character(dat_CT_edgeR$gene_name)
length(all_expressed_genes)

genes_DE_in_males   <- as.character(subset(dat_CT_edgeR, dat_CT_edgeR$CT_M_FDR < 0.05)$gene_name)
genes_DE_in_females <- as.character(subset(dat_CT_edgeR, dat_CT_edgeR$CT_F_FDR < 0.05)$gene_name)

length(genes_DE_in_males)
length(genes_DE_in_females)

######################################################
### Overlap between pop_divergence data (Wiberg et al 2020) and gene expression data

dat_wiberg_2020_CTmin <- read.csv("data/Wiberg_etal_2020/TableS3_bscan_topSNPs_genes_CTmin_DM.csv", header = T)


#### subset by genes that were expressed
dat_wiberg_2020_CTmin$Dmon_ID  <- as.character(dat_wiberg_2020_CTmin$Dmon_ID)
dat_wiberg_2020_CTmin_exp <- subset(dat_wiberg_2020_CTmin, Dmon_ID %in% all_expressed_genes )

length(dat_wiberg_2020_CTmin[,1])
length(dat_wiberg_2020_CTmin_exp[,1])


## test 3-way overlap
res=supertest(list("males_DE" = genes_DE_in_males, "females_DE" = genes_DE_in_females, "wiberg_CTmin" = dat_wiberg_2020_CTmin_exp$Dmon_ID ), n= length(all_expressed_genes))
summary(res)


### plot venn
## 3-way
venn.plot.3 <- venn.diagram(
list("males_DE" = genes_DE_in_males, "females_DE" = genes_DE_in_females, "wiberg_CTmin" = dat_wiberg_2020_CTmin_exp$Dmon_ID ), filename = NULL,
                            cat.col = c("dodgerblue", "firebrick2" , "black"),
                            fill=c("dodgerblue", "firebrick2" , "black"), margin = 0.2)


grid.arrange(gTree(children=venn.plot.3) )

Venn_female_male_wiberg <- arrangeGrob(gTree(children=venn.plot.3),ncol = 1 )
ggsave(file="CT_analysis/EdgeR_out/Venn_female_male_wiberg.pdf", Venn_female_male_wiberg)



### output data
dat_wiberg_2020_CTmin$CT_males_DE   <- ifelse(dat_wiberg_2020_CTmin$Dmon_ID %in% genes_DE_in_males, 1, 0)
dat_wiberg_2020_CTmin$CT_females_DE <- ifelse(dat_wiberg_2020_CTmin$Dmon_ID %in% genes_DE_in_females, 1, 0)
write.csv(dat_wiberg_2020_CTmin, "CT_analysis/EdgeR_out/dat_wiberg_2020_CTmin_DE.csv", row.names = F, quote = F)




######################################################
### Looking at specific gene sets


### innate immune genes 
immune_mel_genes <- read.table("data/Other/innate_immune.txt")
dat_all <-  read.csv("data/Functional_processes/TTT_CT2_sig_CTSB_tidied_with_DAVID.csv", header = T)


#### subset immune genes
head(dat_all)

ss_immune <- subset(dat_all  , Dmel_FBN  %in% immune_mel_genes$V1 )

## innate immune gene DE in females only
subset(ss_immune, ss_immune$CT_F_FDR < 0.05 & ss_immune$CT_M_FDR > 0.05)

## innate immune gene DE in males only
subset(ss_immune, ss_immune$CT_M_FDR < 0.05 & ss_immune$CT_F_FDR > 0.05)



### sex diff genes
tra <- ("FBgn0003741")
fru <- ("FBgn0004652")
dsx <- ("FBgn0000504")

subset(dat_all  , Dmel_FBN  %in% tra) ## tra to lowly expressed in my data
subset(dat_all  , Dmel_FBN  %in% fru)
subset(dat_all  , Dmel_FBN  %in% dsx)


