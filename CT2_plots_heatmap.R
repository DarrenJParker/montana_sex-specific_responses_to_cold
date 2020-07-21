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

dat_logcpm = read.csv("CT_analysis/EdgeR_out/logCPM_all.csv")
head(dat_logcpm, n = 20)


###################################################################################################################################################
### plot heatmap (logFC)

dat1 <- cbind(
dat_logcpm$KL_Temp6_SexF_Rep1,
dat_logcpm$KL_Temp6_SexF_Rep2 ,
dat_logcpm$KL_Temp6_SexF_Rep3 ,
dat_logcpm$KL_Temp19_SexF_Rep1 ,
dat_logcpm$KL_Temp19_SexF_Rep2 ,
dat_logcpm$KL_Temp19_SexF_Rep3 ,
dat_logcpm$KL_Temp6_SexM_Rep1 ,
dat_logcpm$KL_Temp6_SexM_Rep2 ,
dat_logcpm$KL_Temp6_SexM_Rep3 ,
dat_logcpm$KL_Temp19_SexM_Rep1,
dat_logcpm$KL_Temp19_SexM_Rep2 ,
dat_logcpm$KL_Temp19_SexM_Rep3)

rownames(dat1) <- dat_logcpm$genes
colnames(dat1) <- c("KL_Temp6_SexF_Rep1", "KL_Temp6_SexF_Rep2", "KL_Temp6_SexF_Rep3", "KL_Temp19_SexF_Rep1", "KL_Temp19_SexF_Rep2", "KL_Temp19_SexF_Rep3", "KL_Temp6_SexM_Rep1", "KL_Temp6_SexM_Rep2", "KL_Temp6_SexM_Rep3", "KL_Temp19_SexM_Rep1", "KL_Temp19_SexM_Rep2", "KL_Temp19_SexM_Rep3" )


png(filename = "heatmap_all_logcpm.png", width = 5, height = 12, units = "in", pointsize = 12, bg = "white", res = 300)
pheatmap(dat1, clustering_distance_rows = "euclidean", clustering_distance_cols = "correlation", clustering_method = "ward.D2", show_rownames = F)
dev.off()


### this takes a while to run.
N_boot = 1000

Pv_clust_plot <- pvclust(dat1, method.hclust="ward.D2", method.dist="correlation", nboot= N_boot ) 
plot(Pv_clust_plot)

png(filename = paste("PPv_clust_allsamples_logcpm_n=", N_boot, ".png", sep = ""), width = 5, height = 7, units = "in", pointsize = 12, bg = "white", res = 300)
plot(Pv_clust_plot)
dev.off()


	
	
########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "CT2_plots_heatmaps.R_sessionInfo.txt")







