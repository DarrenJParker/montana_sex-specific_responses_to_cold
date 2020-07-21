library(ggplot2)
library(pheatmap)
library(cowplot)
library(grid)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(cowplot)
library(gtable)
library(RColorBrewer)
library(pvclust)
library("VennDiagram")
library("SuperExactTest")
library(raster)


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
 # [1] raster_3.0-2         sp_1.3-1             SuperExactTest_1.0.7 VennDiagram_1.6.20   futile.logger_1.4.3  pvclust_2.0-0        RColorBrewer_1.1-2   gtable_0.3.0         lattice_0.20-38      gridExtra_2.3        cowplot_1.0.0        pheatmap_1.0.12     
# [13] ggplot2_3.2.1       

# loaded via a namespace (and not attached):
 # [1] Rcpp_1.0.2           magrittr_1.5         tidyselect_0.2.5     munsell_0.5.0        colorspace_1.4-1     R6_2.4.0             rlang_0.4.0          dplyr_0.8.3          lambda.r_1.2.3       withr_2.1.2          lazyeval_0.2.2       assertthat_0.2.1    
# [13] tibble_2.1.3         crayon_1.3.4         formatR_1.7          purrr_0.3.2          codetools_0.2-16     futile.options_1.0.1 glue_1.3.1           compiler_3.5.1       pillar_1.4.2         scales_1.0.0         pkgconfig_2.0.2     


################################################################################################################################################
##### read in data

dat_all = read.csv("data/Functional_processes/TTT_CT2_sig_CTSB_tidied_with_DAVID.csv")
head(dat_all , n = 20)


###################################################################################################################################################
### plot heatmap (logFC)

#### interaction
david_inter_1_df <- subset(dat_all, dat_all$david_I1 == 1)

david_inter_1_df_FC <-
cbind(
david_inter_1_df$CT_F_logFC,
david_inter_1_df$CT_M_logFC)
rownames(david_inter_1_df_FC) <- david_inter_1_df$Dmel_FBsym
colnames(david_inter_1_df_FC ) <- c("Females", "Males")

breaksList = c(-4.25, -3.75, -3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25,  0.75,  1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25) 


png(filename = "heatmap_david_inter_1.png", width = 3, height = 4, units = "in", pointsize = 12, bg = "white", res = 300)
pheatmap(david_inter_1_df_FC , clustering_distance_rows = "euclidean", cluster_cols = FALSE, 
clustering_method = "ward.D2", color = colorRampPalette(c("#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D","#67000D"))(length(breaksList)), breaks = breaksList, show_rownames = T,border_color = NA)
dev.off()


#### immune


david_immune_df  <- rbind(subset(dat_all, dat_all$david_B2 == 1), subset(dat_all, dat_all$david_B12 == 1))
david_immune_df_FC <-
cbind(
david_immune_df$CT_F_logFC,
david_immune_df$CT_M_logFC)
rownames(david_immune_df_FC) <- david_immune_df$Dmel_FBsym
colnames(david_immune_df_FC) <- c("Females", "Males")

png(filename = "heatmap_david_immune.png", width = 3, height = 3, units = "in", pointsize = 12, bg = "white", res = 300)
pheatmap(david_immune_df_FC , clustering_distance_rows = "euclidean", cluster_cols = FALSE, 
clustering_method = "ward.D2", color = colorRampPalette(c("#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D","#67000D"))(length(breaksList)), breaks = breaksList, show_rownames = T,border_color = NA)
dev.off()






	
	
	
	
	
	
	
	

	
	
########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "CT2_plots_heatmaps_GOs.R_sessionInfo.txt")







