# montana_sex-specific_responses_to_cold

This is the repository for the collected scripts used in the study:

> Parker, D. J., Envall, T., Ritchie, M. G., and Kankare, M. (2020) Sex-specific responses to cold in Drosophila montana.

## Data

* Phenotype (CTmin) data: 
    * `data/phenotype/Ctmin_data.csv`
    * `data/phenotype/README_Ctmin_data.txt`

* Gene_expression data: 
    * `data/Gene_expression/readcounts/CT_H2E.counts.csv` - readcounts
    * `data/Gene_expression/gff/D.mont_freeze_v1.4.gff.txt` - D. montana annotation info
    * `data/Gene_expression/SAGD_SBgenes/SAGD_00026.csv` - D. virilis sex-bias info from `Shi M-W, Zhang N-A, Shi C-P, Liu C-J, Luo Z-H, Wang D-Y, et al. (2019). SAGD: a comprehensive sex-associated gene database from transcriptomes. Nucleic Acids Res 47: D835–D840`

* Gene expression output files provided for convenience
    * `CT_analysis/EdgeR_out/logCPM_all.csv` - log counts per million for all genes
    * `CT_analysis/EdgeR_out/TTT_CT2_sig_CTSB_tidied.csv` - EdgeR output: expression differences between sexes and temperatures 

* Functional processes
     * `data/Functional_processes/TTT_CT2_sig_CTSB_tidied_with_DAVID.csv` - Gene expression output file, with D. melanogaster orthologs from [FlyMine](https://www.flymine.org/flymine/begin.do), and output from [DAVID](https://david.ncifcrf.gov/) attached.
     * `data/DAVID_genelists` contains D. melanogaster gene IDs used for DAVID analyses:
         * `both_dmel_genelist_DU.txt` - genes DE in response to cold in both males and females 
         * `CT_female_only_dmel_genelist_DU.txt` - genes DE in response to cold in females only
         * `CT_male_only_dmel_genelist_DU.t` - genes DE in response to cold in males only
         * `interaction_dmel_genelist_DU.txt` - Genes showing a significant sex by cold treatment interaction
         * `background_dmel_genelist_DU.txt` - All expressed genes. Used as the 'background list' in DAVID        

## Phenotype analyses

Analyse CT min data.

* `CTmin.R`

## Gene expression analyses

Take read counts and identifies gene expression changes.

* `CT_edgeR.R` 

The output of `CT_edgeR.R` then needs to be brought together with `CT_SB_edgeR_tidier.py` as follows: 

```
python3 CT_SB_edgeR_tidier.py \
-a CT_analysis/EdgeR_out/TTT_SB_6_contrast.csv \
-b CT_analysis/EdgeR_out/TTT_SB_19_contrast.csv \
-f CT_analysis/EdgeR_out/TTT_female_CT_contrast.csv \
-m CT_analysis/EdgeR_out/TTT_male_CT_contrast.csv \
-i CT_analysis/EdgeR_out/TTT_CT_sex_temp_inter.csv \
-g data/Gene_expression/gff/D.mont_freeze_v1.4.gff.txt \
-v data/Gene_expression/SAGD_SBgenes/SAGD_00026.csv \
-o CT_analysis/EdgeR_out/TTT_CT2_sig

```

Note the joined file `CT_analysis/EdgeR_out/TTT_CT2_sig_CTSB_tidied.csv` is provided for convenience along with a description file `CT_analysis/EdgeR_out/README_TTT_CT2_sig_CTSB_tidied.txt`

Plot gene expression data

* `CT2_plots_and_hists.R`
* `CT2_plots_heatmap.R`

## Functional processes

Plot gene expression of functional clusters
* CT2_plots_heatmap_GOs.R


## Infomation on running scripts

* All scripts should be run from the directory they are in. Output will be stored as the code is run. 
* All python scripts were made using python 3.5. All contain help information which can be displayed by specifying -h as a command line argument.



