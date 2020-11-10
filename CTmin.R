## CTmin

library(ggplot2)
library(fitdistrplus)
library(car)
library(lme4)
print (sessionInfo())

# R version 3.5.1 (2018-07-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.15.7

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] lme4_1.1-21         Matrix_1.2-17       car_3.0-3           carData_3.0-2       fitdistrplus_1.0-14 npsurv_0.4-0        lsei_1.2-0          survival_3.1-7      MASS_7.3-51.4       ggplot2_3.3.2      

# loaded via a namespace (and not attached):
 # [1] zip_2.0.3         Rcpp_1.0.2        nloptr_1.2.1      cellranger_1.1.0  pillar_1.4.2      compiler_3.5.1    tools_3.5.1       forcats_0.4.0     boot_1.3-23       nlme_3.1-141      lifecycle_0.2.0   tibble_2.1.3      gtable_0.3.0     
# [14] lattice_0.20-38   pkgconfig_2.0.2   rlang_0.4.8       openxlsx_4.1.0.1  curl_4.0          haven_2.1.1       rio_0.5.16        withr_2.1.2       dplyr_1.0.2       generics_0.0.2    vctrs_0.3.4       hms_0.5.1         grid_3.5.1       
# [27] tidyselect_1.1.0  glue_1.4.2        data.table_1.12.2 R6_2.4.0          readxl_1.3.1      foreign_0.8-72    minqa_1.2.4       purrr_0.3.2       magrittr_1.5      scales_1.0.0      splines_3.5.1     abind_1.4-5       colorspace_1.4-1 
# [40] munsell_0.5.0     crayon_1.3.4     


##### data
## see  data/phenotypic/README_Ctmin_data.txt for explination of columns
dat1 <- read.table("data/phenotypic/Ctmin_data.csv", sep = ",", header = TRUE)
head(dat1)

dat1$group <- paste(dat1$Sex, dat1$Cold_treatment, sep = "_")

### add 5 to all values to avoid logging neg numbers
min(dat1$Ctmin)
dat1$Ctmin_5 <- dat1$Ctmin + 5
min(dat1$Ctmin_5)

################
# plot

P1 <- ggplot(dat1, aes(group, Ctmin)) + 
	theme_bw() +
	geom_boxplot(aes(fill = factor(Sex)), outlier.shape = NA) +
	geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), aes(fill = Sex), stroke = 1, dotsize=0.5 , binwidth = 0.1) + 
	ylab ("CT min") +
	xlab ("Treatment group") + 
	scale_fill_manual(values=c("firebrick2", "royalblue2"))


######### anova

### non-norm dist 

hist(dat1$Ctmin)
shapiro.test(dat1$Ctmin) # W = 0.94699, p-value = 4.247e-05
m1 = glm(dat1$Ctmin ~ dat1$Sex * dat1$CA_treatment + dat1$weight_mg)	
summary(m1)
plot(m1)

### looks log-norm

fit_n   <- fitdist(as.numeric(na.omit(dat1$Ctmin_5)), "norm")
fit_ln  <- fitdist(as.numeric(na.omit(dat1$Ctmin_5)), "lnorm")
fit_g   <- fitdist(as.numeric(na.omit(dat1$Ctmin_5)), "gamma")

denscomp(list(fit_n, fit_ln, fit_g), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n, fit_ln, fit_g), legendtext = c("norm", "lnorm", "gamma"))

shapiro.test(log(dat1$Ctmin_5)) # W = 0.99052, p-value = 0.4812

m2 = glm(log(dat1$Ctmin_5) ~ dat1$Sex * dat1$Cold_treatment + dat1$weight_mg)	## I am happy to not reduce the model, but will do so just to check
m3 = glm(log(dat1$Ctmin_5) ~ dat1$Sex + dat1$Cold_treatment + dat1$weight_mg)	
m4 = glm(log(dat1$Ctmin_5) ~ dat1$Sex + dat1$Cold_treatment)	
m5 = glm(log(dat1$Ctmin_5) ~ dat1$Cold_treatment)	

Anova(m2, type = 3, test.statistic = "F")

# # Analysis of Deviance Table (Type III tests)

# Response: log(dat1$Ctmin_5)
# Error estimate based on Pearson residuals 

                            # Sum Sq  Df F values  Pr(>F)   
# dat1$Sex                    0.0038   1   0.0391 0.84360   
# dat1$CA_treatment           1.0815   1  11.1732 0.00108 **
# dat1$weight_mg              0.0838   1   0.8657 0.35386   
# dat1$Sex:dat1$CA_treatment  0.1345   1   1.3900 0.24053   
# Residuals                  12.7769 132       

Anova(m3, type = 3, test.statistic = "F")
Anova(m4, type = 3, test.statistic = "F")
Anova(m5, type = 3, test.statistic = "F")


### doesnt change results - acclimation has an effect, but no effect of sex or weight 


####################
### GLMM to account for batch - similar results as from the glm's above - batch has a significant effect

mix_2 = lmer(log(dat1$Ctmin_5) ~ dat1$Sex * dat1$Cold_treatment + dat1$weight_mg +(1|dat1$Batch)) 
Anova(mix_2  ,type="III") ## sig of fixed effects with Wald

# Analysis of Deviance Table (Type III Wald chisquare tests)

# Response: log(dat1$Ctmin_5)
                             # Chisq Df Pr(>Chisq)    
# (Intercept)                58.0928  1    2.5e-14 ***
# dat1$Sex                    0.6822  1    0.40883    
# dat1$CA_treatment           4.9315  1    0.02637 *  
# dat1$weight_mg              0.2171  1    0.64130    
# dat1$Sex:dat1$CA_treatment  0.0067  1    0.93485    

### sig of batch by LRT
LL_mix1 = logLik(m2)
LL_mix2 = logLik(mix_2)
LRT_1_2 = 2 * (as.numeric(LL_mix1) - as.numeric(LL_mix2))
pchisq(as.numeric(LRT_1_2), df=1, lower.tail=F) ## 0.0002363689








