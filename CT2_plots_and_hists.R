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
# > 


################################################################################################################################################
##### read in data

dat_CT2_all = read.csv("CT_analysis/EdgeR_out/TTT_CT2_sig_CTSB_tidied.csv") ## produced by CT_SB_edgeR_tidier.py
### see CT_analysis/EdgeR_out/README_TTT_CT2_sig_CTSB_tidied.txt for info on columns

head(dat_CT2_all, n = 20)
setwd("CT_analysis/EdgeR_out")


#################################################################################################################################################
###### refine sex bias by FC (not log FC)
	
FC_cutoff = 2

dat_CT2_all$SB_6_sexbias2       <- ifelse(2^(sqrt(dat_CT2_all$SB_6_logFC         * dat_CT2_all$SB_6_logFC))         > FC_cutoff , as.character(dat_CT2_all$SB_6_sexbias)      , "Unbiased")
dat_CT2_all$SB_19_sexbias2      <- ifelse(2^(sqrt(dat_CT2_all$SB_19_logFC        * dat_CT2_all$SB_19_logFC))        > FC_cutoff , as.character(dat_CT2_all$SB_19_sexbias)     , "Unbiased")
dat_CT2_all$vir_sexbias2        <- ifelse(2^(sqrt(dat_CT2_all$vir_log2FoldChange * dat_CT2_all$vir_log2FoldChange)) > FC_cutoff , as.character(dat_CT2_all$vir_sexbias), "Unbiased")

#### SB in both conditions

dat_CT2_all$SB_619_sexbias_j <- paste(dat_CT2_all$SB_6_sexbias, dat_CT2_all$SB_19_sexbias, sep = "_")
dat_CT2_all$SB_619_sexbias   <- ifelse(dat_CT2_all$SB_619_sexbias_j == "female_biased_female_biased", "female_biased", 
					             ifelse(dat_CT2_all$SB_619_sexbias_j == "male_biased_male_biased", "male_biased", "Unbiased"))

dat_CT2_all$SB_619_sexbias2_j <- paste(dat_CT2_all$SB_6_sexbias2, dat_CT2_all$SB_19_sexbias2, sep = "_")
dat_CT2_all$SB_619_sexbias2   <- ifelse(dat_CT2_all$SB_619_sexbias2_j == "female_biased_female_biased", "female_biased", 
					             ifelse(dat_CT2_all$SB_619_sexbias2_j == "male_biased_male_biased", "male_biased", "Unbiased"))
### as ordered factor

dat_CT2_all$SB_619_sexbias_ord           <- ordered(dat_CT2_all$SB_619_sexbias, levels = c("male_biased", "Unbiased", "female_biased"))   
dat_CT2_all$SB_619_sexbias2_ord          <- ordered(dat_CT2_all$SB_619_sexbias2, levels = c("male_biased", "Unbiased", "female_biased"))   

dat_CT2_all$vir_sexbias_ord           <- ordered(dat_CT2_all$vir_sexbias, levels = c("male_biased", "Unbiased", "female_biased"))   
dat_CT2_all$vir_sexbias2_ord          <- ordered(dat_CT2_all$vir_sexbias2, levels = c("male_biased", "Unbiased", "female_biased"))   

### overlapping vir and mont

dat_CT2_all$SB_619_virsexbias2_j <- paste(dat_CT2_all$SB_619_sexbias2, dat_CT2_all$vir_sexbias2, sep = "_")
dat_CT2_all$SB_619_virsexbias2     <- ifelse(dat_CT2_all$SB_619_virsexbias2_j == "female_biased_female_biased", "female_biased", 
					                  ifelse(dat_CT2_all$SB_619_virsexbias2_j == "male_biased_male_biased", "male_biased", "Unbiased"))
					                
					                
dat_CT2_all$SB_619_virsexbias2_ord          <- ordered(dat_CT2_all$SB_619_virsexbias2, levels = c("male_biased", "Unbiased", "female_biased"))   					                


         
##############################################################################################################################################
############### Do males and females use the same genes in response to cold?
###############################################################################################################################################


### raw numbers of DE genes

female_sig_DE <- subset(dat_CT2_all, dat_CT2_all$CT_F_FDR < 0.05 )[,1]
male_sig_DE <- subset(dat_CT2_all, dat_CT2_all$CT_M_FDR < 0.05 )[,1]
inter_sig_DE <- subset(dat_CT2_all, dat_CT2_all$CT_SB_inter_FDR < 0.05 )[,1]


### Venn

three_DE_venn <- function(Female,Male,Inter,title){
	venny.plot <- venn.diagram(
	list("Males" = Male, "Female" = Female, "Interaction" = Inter), filename = NULL,
                            fill = c("dodgerblue", "firebrick2", "darkorange1"),
                            cat.col = c("dodgerblue", "firebrick2", "darkorange1"),
                            margin = 0.6, cat.dist = 0.23, main = title, main.pos = c(0.5,0.8), main.cex = 2, main.fontface = "bold", cat.cex = 2)
	return(venny.plot)
}

DE_Venn_female_male_inter <- three_DE_venn(female_sig_DE, male_sig_DE, inter_sig_DE, "cc")

grid.arrange(gTree(children=DE_Venn_female_male_inter),ncol = 1 )

# ## to output to pdf

DE_Venn_female_male_inter_1 <- arrangeGrob(gTree(children=DE_Venn_female_male_inter),ncol = 1 )
ggsave(file="DE_Venn_female_male_inter_1.pdf", DE_Venn_female_male_inter_1)

### test overlap

head(dat_CT2_all)
FM_DE_lists <- list(female_sig_DE, male_sig_DE )
res=supertest(FM_DE_lists , n= length(dat_CT2_all[,1]))
summary(res)

# A msets object with 2 sets: Set1 Set2 
# Background size: 9338 
# Summary of intersections:
   # Intersections Degree Observed.Overlap Expected.Overlap       FE       P.value                 Elements
# 01          Set2      1             1236               NA       NA            NA augustus_masked-scaf ...
# 10          Set1      1             1062               NA       NA            NA augustus_masked-scaf ...
# 11   Set1 & Set2      2              407         140.5689 2.895378 1.739285e-110 augustus_masked-scaf ...

## plot correlations

### class 

dat_CT2_all$CT_F_sigclass <- ifelse(dat_CT2_all$CT_F_FDR < 0.05, "sig_F", "non-sig_F")
dat_CT2_all$CT_M_sigclass <- ifelse(dat_CT2_all$CT_M_FDR < 0.05, "sig_M", "non-sig_M")
dat_CT2_all$CT_FM_sigclass <- paste(dat_CT2_all$CT_F_sigclass, dat_CT2_all$CT_M_sigclass, sep = "_")
dat_CT2_all$CT_inter_sigclass <- ifelse(dat_CT2_all$CT_SB_inter_FDR < 0.05 , "sig_inter", "non-sig_inter")

max(dat_CT2_all$CT_F_logFC)
max(dat_CT2_all$CT_M_logFC)
min(dat_CT2_all$CT_F_logFC)
min(dat_CT2_all$CT_M_logFC)


head(dat_CT2_all)

#### alter order for plotting

dat_CT2_all_ordered <- rbind(
subset(dat_CT2_all, dat_CT2_all$CT_FM_sigclass == "non-sig_F_non-sig_M"),
subset(dat_CT2_all, dat_CT2_all$CT_FM_sigclass == "sig_F_non-sig_M"),
subset(dat_CT2_all, dat_CT2_all$CT_FM_sigclass == "non-sig_F_sig_M"),
subset(dat_CT2_all, dat_CT2_all$CT_FM_sigclass == "sig_F_sig_M")
)

dat_CT2_allsig_ordered <- rbind(
subset(dat_CT2_all, dat_CT2_all$CT_FM_sigclass == "sig_F_non-sig_M"),
subset(dat_CT2_all, dat_CT2_all$CT_FM_sigclass == "non-sig_F_sig_M"),
subset(dat_CT2_all, dat_CT2_all$CT_FM_sigclass == "sig_F_sig_M")
)

dat_CT2_DEinbothsexes_ordered <- rbind(
subset(dat_CT2_all, dat_CT2_all$CT_FM_sigclass == "sig_F_sig_M")
)

length(dat_CT2_all_ordered[,1])
length(dat_CT2_allsig_ordered[,1])
length(dat_CT2_DEinbothsexes_ordered[,1])
length(dat_CT2_all[,1])


### without reg line and all data
MF_Plot4 <- 
ggplot(dat_CT2_all_ordered, aes(x=CT_M_logFC, CT_F_logFC, col = CT_FM_sigclass)) + 
geom_point(aes(shape = CT_inter_sigclass), size = 1,alpha = 0.5) + 
ylim(-10,10) + xlim(-10,10)  + theme_bw() + 
scale_colour_manual(values = c("grey","dodgerblue4","firebrick3","orange"))  + 
xlab("Log2 fold expression change in females") + ylab("Log2 fold expression change in males")




png(filename = "MF_Plot4.png", width = 6, height = 4, units = "in", bg = "white", res = 300)
MF_Plot4
dev.off()


### corr coffs

cor.test(dat_CT2_all_ordered$CT_M_logFC , dat_CT2_all_ordered$CT_F_logFC, method = "spearman")
# 0.5872869, p-value < 2.2e-16

cor.test(dat_CT2_allsig_ordered$CT_M_logFC , dat_CT2_allsig_ordered$CT_F_logFC, method = "spearman")
# 0.7311189, p-value < 2.2e-16

cor.test(dat_CT2_DEinbothsexes_ordered$CT_M_logFC , dat_CT2_DEinbothsexes_ordered$CT_F_logFC, method = "spearman")
# 0.9178907, p-value < 2.2e-16


##############################################################################################################################################
### interaction genes


dat_CT2_sig_inter <- subset(dat_CT2_all, dat_CT2_all$CT_inter_sigclass == "sig_inter")
length(dat_CT2_sig_inter[,1])

head(dat_CT2_sig_inter)


dat_CT2_sig_inter_l <- as.data.frame(cbind(
c(dat_CT2_sig_inter$gene_name, dat_CT2_sig_inter$gene_name),
c(dat_CT2_sig_inter$CT_F_logFC, dat_CT2_sig_inter$CT_M_logFC),
c(rep("F", length(dat_CT2_sig_inter[,1])), rep("M", length(dat_CT2_sig_inter[,1]))),
c(dat_CT2_sig_inter$CT_FM_sigclass, dat_CT2_sig_inter$CT_FM_sigclass)
))

colnames(dat_CT2_sig_inter_l) <- c("gene_name", "logFC", "sex", "CT_FM_sigclass")
dat_CT2_sig_inter_l$logFC <- as.numeric(as.character(dat_CT2_sig_inter_l$logFC))

dat_CT2_sig_inter_l$CT_FM_sigclass2 <- ifelse(dat_CT2_sig_inter_l$CT_FM_sigclass == "sig_F_sig_M", "Both", "Single")

### plotting subsets so the Both or on top
Inter_P1 <- ggplot()  + geom_hline(yintercept=0) +  geom_line(data=subset(dat_CT2_sig_inter_l, dat_CT2_sig_inter_l$CT_FM_sigclass2 != "Both"), aes( x=sex, y=logFC, group=gene_name, col=CT_FM_sigclass2))  + theme_bw() + 
           geom_line(data=subset(dat_CT2_sig_inter_l, dat_CT2_sig_inter_l$CT_FM_sigclass2 == "Both"), aes( x=sex, y=logFC, group=gene_name, col=CT_FM_sigclass2))

Inter_P2 <- Inter_P1 + scale_color_manual(values=c('#E69F00', '#999999')) + ylim(c(-9,9)) + scale_x_discrete(expand = c(0, .2))


pdf("Inter_P2.pdf", width = 6, height = 8)
Inter_P2
dev.off()
getwd() ## where has my plot gone....?




##############################################################################################################################################
############### Shifts in SB genes
################################################################################################################################################


##### plot boxplots

plot_hist_and_box = function(df, sb_type){
	
	##########################################################################################################################################
	### boxplots
	
	xlim_val_low  = -4
	xlim_val_high = 4
	
	print(sb_type)
	female_box_p1 <- ggplot(eval(parse(text=paste(df, sep=""))), aes(eval(parse(text=paste(df,"$",sb_type, sep=""))), CT_F_logFC)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(eval(parse(text=paste(df,"$",sb_type, sep=""))))), position=position_dodge(0.65), width = 0.45, outlier.size = 1, lwd = 0.5, fatten = 1) +
		ylab ("log2 fold change in expression ") +
		xlab ("") + 
		scale_fill_manual(values=c("royalblue2", "grey", "firebrick2")) + 
		#geom_hline(yintercept = median(subset(dat_CT2_all,  dat_CT2_all$SB_overall_sexbias2_ord == "Unbiased")$CT_M_logFC)) + 
		labs(fill='') +
		ggtitle("") + coord_flip(ylim=c(xlim_val_low, xlim_val_high))	

	male_box_p1 <- ggplot(eval(parse(text=paste(df, sep=""))), aes(eval(parse(text=paste(df,"$",sb_type, sep=""))), CT_M_logFC)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(eval(parse(text=paste(df,"$",sb_type, sep=""))))), position=position_dodge(0.65), width = 0.45, outlier.size = 1, lwd = 0.5, fatten = 1) +
		ylab ("log2 fold change in expression ") +
		xlab ("") + 
		scale_fill_manual(values=c("royalblue2", "grey", "firebrick2")) + 
		#geom_hline(yintercept = median(subset(dat_CT2_all,  dat_CT2_all$SB_overall_sexbias2_ord == "Unbiased")$CT_M_logFC)) + 
		labs(fill='') +
		ggtitle("") + coord_flip(ylim=c(xlim_val_low, xlim_val_high))	

	####################################################################################################################################
	### hist

	### ALL
	#FB <- subset(dat_CT2_all, dat_CT2_all$SB_overall_sexbias2 == "female_biased")
	FB_all <- subset(eval(parse(text=paste(df, sep=""))), eval(parse(text=paste(df,"$",sb_type, sep=""))) == "female_biased")
	MB_all <- subset(eval(parse(text=paste(df, sep=""))), eval(parse(text=paste(df,"$",sb_type, sep=""))) == "male_biased")	
	UB_all <- subset(eval(parse(text=paste(df, sep=""))), eval(parse(text=paste(df,"$",sb_type, sep=""))) == "Unbiased")
	


	######
	
	p0_female_all_hist = ggplot(eval(parse(text=paste(df, sep=""))),aes(x=CT_F_logFC)) + 
		theme_classic() +	
		geom_histogram(aes(y = ..density..), data= MB_all, color="darkblue", fill="blue", binwidth=0.1 ,alpha = 0.2) +
		geom_histogram(aes(y = ..density..), data= FB_all,color="darkred", fill="red", binwidth=0.1,alpha = 0.2) +
		geom_density(data=MB_all, color="darkblue", size = 1) + 
		geom_density(data=FB_all, color="darkred", size = 1) + 
		coord_cartesian(xlim=c(xlim_val_low, xlim_val_high))  + 
		scale_y_continuous(expand = c(0,0.05)) +  #### makes it start at y = 0
		geom_vline(xintercept = 0,  linetype = "longdash") +
		geom_hline(yintercept = 0) 

	p0_male_all_hist = ggplot(eval(parse(text=paste(df, sep=""))),aes(x=CT_M_logFC)) + 
		theme_classic() +	
		geom_histogram(aes(y = ..density..), data= MB_all, color="darkblue", fill="blue", binwidth=0.1 ,alpha = 0.2) +
		geom_histogram(aes(y = ..density..), data= FB_all,color="darkred", fill="red", binwidth=0.1,alpha = 0.2) +
		geom_density(data=MB_all, color="darkblue", size = 1) + 
		geom_density(data=FB_all, color="darkred", size = 1) + 
		coord_cartesian(xlim=c(xlim_val_low, xlim_val_high))  + 
		scale_y_continuous(expand = c(0,0.05)) +  #### makes it start at y = 0
		geom_vline(xintercept = 0,  linetype = "longdash") +
		geom_hline(yintercept = 0) 


	######################################################################################################
	###### wilcox
	

	FB_all_female_wil <- wilcox.test(FB_all$CT_F_logFC, UB_all$CT_F_logFC, paired = FALSE)
	MB_all_female_wil <- wilcox.test(MB_all$CT_F_logFC, UB_all$CT_F_logFC, paired = FALSE)
	FB_all_male_wil   <- wilcox.test(FB_all$CT_M_logFC, UB_all$CT_M_logFC, paired = FALSE)
	MB_all_male_wil   <- wilcox.test(MB_all$CT_M_logFC, UB_all$CT_M_logFC, paired = FALSE)

	wilcox_all_out <- as.data.frame(
	cbind(
	c(FB_all_female_wil$statistic, MB_all_female_wil$statistic, FB_all_male_wil$statistic, MB_all_male_wil$statistic),
	c(FB_all_female_wil$p.value, MB_all_female_wil$p.value, FB_all_male_wil$p.value, MB_all_male_wil$p.value),
	c(median(FB_all$CT_F_logFC), median(MB_all$CT_F_logFC), median(FB_all$CT_M_logFC), median(MB_all$CT_M_logFC)),
	c(median(UB_all$CT_F_logFC), median(UB_all$CT_F_logFC), median(UB_all$CT_M_logFC), median(UB_all$CT_M_logFC))	
	)
	)
	
	colnames(wilcox_all_out) <- c("W","p","med_SB","med_UB")
	rownames(wilcox_all_out) <- c("Female_FB-UB","Female_MB-UB","Male_FB-UB","Male_MB-UB")

	wilcox_all_out$FDR <- p.adjust(wilcox_all_out$p, method = "BH")
	
	
	outlist <- list("female_box_p1" = female_box_p1,         "male_box_p1" = male_box_p1,
	                "p0_female_all_hist" = p0_female_all_hist, "p0_male_all_hist" = p0_male_all_hist,
	                "wilcox_all_out" =  wilcox_all_out,
	                "FB_all" = FB_all, "MB_all" = MB_all, "UB_all" = UB_all)
	return(outlist)	

}


##### SB 

ALL_female_SB619_2_box <- plot_hist_and_box("dat_CT2_all", "SB_619_sexbias2_ord")$female_box_p1
ALL_male_SB619_2_box   <- plot_hist_and_box("dat_CT2_all", "SB_619_sexbias2_ord")$male_box_p1

ALL_female_SB619_2_hist <- plot_hist_and_box("dat_CT2_all", "SB_619_sexbias2_ord")$p0_female_all_hist
ALL_male_SB619_2_hist   <- plot_hist_and_box("dat_CT2_all", "SB_619_sexbias2_ord")$p0_male_all_hist

plot_grid(ALL_female_SB619_2_hist, ALL_male_SB619_2_hist, ALL_female_SB619_2_box, ALL_male_SB619_2_box, ncol = 2, nrow = 2)

pdf("ALL_SB619_2_hist.pdf", width = 10, height = 4)
plot_grid(ALL_female_SB619_2_hist, ALL_male_SB619_2_hist, ncol = 2, nrow = 1)
dev.off()
getwd() ## where has my plot gone....?

pdf("ALL_SB619_2_box.pdf", width = 14, height = 1.4)
plot_grid(ALL_female_SB619_2_box, ALL_male_SB619_2_box, ncol = 2, nrow = 1)
dev.off()
getwd() ## where has my plot gone....?





###### wilcox

wilcox_619_SB2_all <- plot_hist_and_box("dat_CT2_all", "SB_619_sexbias2_ord")$wilcox_all_out
write.csv(wilcox_619_SB2_all, "wilcox_619_SB2_all.csv")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### SHIFT FROM UB

FB_all_619_SB2 <- plot_hist_and_box("dat_CT2_all", "SB_619_sexbias2_ord")$FB_all 
MB_all_619_SB2 <- plot_hist_and_box("dat_CT2_all", "SB_619_sexbias2_ord")$MB_all 

### females 
median(FB_all_619_SB2$CT_F_logFC) 
median(MB_all_619_SB2$CT_F_logFC)

### make rel to unbiased med
median(FB_all_619_SB2$CT_F_logFC - -0.08923119) 
median(MB_all_619_SB2$CT_F_logFC - -0.08923119)

FB_F_UBadj <- FB_all_619_SB2$CT_F_logFC - -0.08923119
MB_F_UBadj <- MB_all_619_SB2$CT_F_logFC - -0.08923119

c(
median(FB_all_619_SB2$CT_F_logFC - -0.08923119),
wilcox.test(FB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$estimate,
wilcox.test(FB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[1],
wilcox.test(FB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[2])

c(
median(MB_all_619_SB2$CT_F_logFC - -0.08923119),
wilcox.test(MB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$estimate,
wilcox.test(MB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[1],
wilcox.test(MB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[2])

### males 

median(FB_all_619_SB2$CT_M_logFC) 
median(MB_all_619_SB2$CT_M_logFC)

### make rel to unbiased med
median(FB_all_619_SB2$CT_M_logFC - -0.02756493) 
median(MB_all_619_SB2$CT_M_logFC - -0.02756493)

FB_M_UBadj <- FB_all_619_SB2$CT_M_logFC - -0.02756493
MB_M_UBadj <- MB_all_619_SB2$CT_M_logFC - -0.02756493

c(
median(FB_all_619_SB2$CT_M_logFC - -0.02756493),
wilcox.test(FB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$estimate,
wilcox.test(FB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[1],
wilcox.test(FB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[2])

c(
median(MB_all_619_SB2$CT_M_logFC - -0.02756493),
wilcox.test(MB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$estimate,
wilcox.test(MB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[1],
wilcox.test(MB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[2])


######## join together

FB_all_619_SB2_CIs_UBadj <- as.data.frame(rbind(
c(
median(FB_all_619_SB2$CT_F_logFC - -0.08923119),
wilcox.test(FB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$estimate,
wilcox.test(FB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[1],
wilcox.test(FB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[2]),

c(
median(MB_all_619_SB2$CT_F_logFC - -0.08923119),
wilcox.test(MB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$estimate,
wilcox.test(MB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[1],
wilcox.test(MB_F_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[2]),

c(
median(FB_all_619_SB2$CT_M_logFC - -0.02756493),
wilcox.test(FB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$estimate,
wilcox.test(FB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[1],
wilcox.test(FB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[2]),

c(
median(MB_all_619_SB2$CT_M_logFC - -0.02756493),
wilcox.test(MB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$estimate,
wilcox.test(MB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[1],
wilcox.test(MB_M_UBadj, paired = FALSE, conf.int=TRUE,conf.level=0.95, correct=TRUE)$conf.int[2])))


rownames(FB_all_619_SB2_CIs_UBadj ) <- c("FB_females", "MB_females", "FB_males", "MB_males")
colnames(FB_all_619_SB2_CIs_UBadj ) <- c("med", "pseudomed", "lower", "upper")

FB_all_619_SB2_CIs_UBadj $class <- rownames(FB_all_619_SB2_CIs_UBadj )
FB_all_619_SB2_CIs_UBadj $class_ord <- ordered(FB_all_619_SB2_CIs_UBadj $class, levels = c("MB_males", "FB_males",   "MB_females" , "FB_females"))
FB_all_619_SB2_CIs_UBadj $bias      <- c("MB", "FB",   "MB" , "FB")

FB_all_619_SB2_CIs_P1 <- ggplot() + 
geom_errorbar(data=FB_all_619_SB2_CIs_UBadj , mapping=aes(x=class_ord, ymin=upper, ymax=lower,  color= bias), width=0.2, size=1) +  theme_bw() + 
geom_point(data=FB_all_619_SB2_CIs_UBadj , mapping=aes(x = class_ord, y=pseudomed), size=4, shape=21, fill="white") + coord_flip() + ylim(c(-0.1, 0.1)) + geom_hline(yintercept=0) 


write.csv(FB_all_619_SB2_CIs_UBadj , "FB_all_619_SB2_CIs_UBadj .csv")


pdf("FB_all_619_SB2_CIs_P1.pdf", width = 10, height = 3)
FB_all_619_SB2_CIs_P1
dev.off()
getwd() ## where has my plot gone....?







##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### SB using genes SB in vir and mont

head(dat_CT2_all)
ALL_female_vir_mont619_SB_2_box <- plot_hist_and_box("dat_CT2_all", "SB_619_virsexbias2_ord")$female_box_p1
ALL_male_vir_mont619_SB_2_box   <- plot_hist_and_box("dat_CT2_all", "SB_619_virsexbias2_ord")$male_box_p1

ALL_female_vir_mont619_SB_2_hist <- plot_hist_and_box("dat_CT2_all", "SB_619_virsexbias2_ord")$p0_female_all_hist
ALL_male_vir_mont619_SB_2_hist   <- plot_hist_and_box("dat_CT2_all", "SB_619_virsexbias2_ord")$p0_male_all_hist


pdf("ALL_vir_mont619_SB_2_hist.pdf", width = 10, height = 4)
plot_grid(ALL_female_vir_mont619_SB_2_hist , ALL_male_vir_mont619_SB_2_hist, ncol = 2, nrow = 1)
dev.off()
getwd() ## where has my plot gone....?

pdf("ALL_vir_mont619_SB_2_box.pdf", width = 14, height = 1.4)
plot_grid(ALL_female_vir_mont619_SB_2_box , ALL_male_vir_mont619_SB_2_box , ncol = 2, nrow = 1)
dev.off()
getwd() ## where has my plot gone....?

wilcox_vir_mont619_SB_2_all <- plot_hist_and_box("dat_CT2_all", "SB_619_virsexbias2_ord")$wilcox_all_out

write.csv(wilcox_vir_mont619_SB_2_all, "wilcox_vir_mont619_SB_2_all.csv")




################################################################################################################################################
############### Vir and mont SB genes
################################################################################################################################################

vir_MB_2_genes <- subset(dat_CT2_all, dat_CT2_all$vir_sexbias2_ord == "male_biased")[,1]
SB619_MB_2_genes <- subset(dat_CT2_all, dat_CT2_all$SB_619_sexbias2_ord == "male_biased")[,1]

vir_FB_2_genes <- subset(dat_CT2_all, dat_CT2_all$vir_sexbias2_ord == "female_biased")[,1]
SB619_FB_2_genes <- subset(dat_CT2_all, dat_CT2_all$SB_619_sexbias2_ord == "female_biased")[,1]

length(vir_MB_2_genes )


male_DE_venn <- function(DE_vect,inter_vect,title){

	venny.plot <- venn.diagram(
	list("DE" = DE_vect, "Interaction" = inter_vect), filename = NULL,
                            cat.col = c( "blue",   "green"),
                            fill=c("blue",   "green"), margin = 0.2, main = title, main.pos = c(0.5,0.8), cat.dist = 0.03, main.cex = 2, main.fontface = "bold", cat.cex = 0)
	return(venny.plot)
}


female_DE_venn <- function(DE_vect,inter_vect,title){

	venny.plot <- venn.diagram(
	list("DE" = DE_vect, "Interaction" = inter_vect), filename = NULL,
                            cat.col = c( "red",   "yellow"),
                            fill=c("red",   "yellow"), margin = 0.2, main = title, main.pos = c(0.5,0.8), cat.dist = 0.03, main.cex = 2, main.fontface = "bold", cat.cex = 0)
	return(venny.plot)
}


MB_vir_SB619_venn <- male_DE_venn(vir_MB_2_genes, SB619_MB_2_genes, "Male-biased genes")
grid.arrange(gTree(children=MB_vir_SB619_venn),ncol = 1 )



FB_vir_SB619_venn <- female_DE_venn(vir_FB_2_genes, SB619_FB_2_genes, "Female-biased genes")
grid.arrange(gTree(children=FB_vir_SB619_venn),ncol = 1 )



### lots of overlap

# ## to output to pdf

MB_vir_SB619_venn_out <- arrangeGrob(gTree(children=MB_vir_SB619_venn),ncol = 1 )
ggsave(file="MB_vir_SB619_venn_out.pdf", MB_vir_SB619_venn_out, width = 7, height = 7)

FB_vir_SB619_venn_out <- arrangeGrob(gTree(children=FB_vir_SB619_venn),ncol = 1 )
ggsave(file="FB_vir_SB619_venn_out.pdf", FB_vir_SB619_venn_out,  width = 7, height = 7)




########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "CT2_plots_and_hists.R_sessionInfo.txt")









