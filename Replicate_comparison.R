# Mouse brain samples (MOB + Hippocampus, standard Visium and Ex-ST) Stereoscope Analysis

## Load packages
library(STutility)
library(ggplot2)
library(magrittr)
library(ggpubr)

## Load Visium data - space ranger output files
samples <- list.files(path = "/home/st-analysis_home/zaneta.andrusivova/projects/expST/220629_omni_mouse_brain/spaceranger_output_allsamples/", pattern = "raw_feature_bc_matrix.h5", recursive = TRUE, full.names = TRUE)
imgs <- list.files(path = "/home/st-analysis_home/zaneta.andrusivova/projects/expST/220629_omni_mouse_brain/spaceranger_output_allsamples/", pattern = "tissue_hires_image.png", recursive = TRUE, full.names = TRUE)
spotfiles <- list.files(path = "/home/st-analysis_home/zaneta.andrusivova/projects/expST/220629_omni_mouse_brain/spaceranger_output_allsamples/", pattern = "tissue_positions_list.csv", recursive = TRUE, full.names = TRUE)
json <- list.files(path = "/home/st-analysis_home/zaneta.andrusivova/projects/expST/220629_omni_mouse_brain/spaceranger_output_allsamples/", pattern = "scalefactors_json.json", recursive = TRUE, full.names = TRUE)


## Create an infoTable, add information about sample assay, sample region and sample ID
infoTable <- data.frame(samples, imgs, spotfiles, json, stringsAsFactors = F,
                        assay = c("Visium", "Ex-ST", "Ex-ST", "Ex-ST", "Ex-ST", "Visium" ),
                        anat = c("Hippocampus", "Hippocampus", "Hippocampus", "MOB", "MOB", "MOB"),
                        sample_id = c("V_HC", "Ex-ST_HC1", "Ex-ST_HC2", "Ex-ST_MOB1", "Ex-ST_MOB2", "V_MOB"))


## Create Seurat object
se <- InputFromTable(infoTable, disable.subset = TRUE)


## Remove spots with less than 100 genes
se <- SubsetSTData(se, expression = nFeature_RNA > 100)


## Violin plots for Ex-ST replicates
## Subset Seurat object for Ex-ST MOB or Hippocampus
se_exp <- SubsetSTData(se, assay == "Ex-ST")
se_mob <- SubsetSTData(se_exp, anat == "MOB")
se_HC <- SubsetSTData(se_exp, anat == "Hippocampus")

## Plot UMI/nGenes per replicate (MOB)

pdf("UMI_MOB.pdf")
VlnPlot(se_mob, features = "nCount_RNA", pt.size = 0, group.by = "sample_id")
dev.off()


pdf("nGenes_MOB.pdf")
VlnPlot(se_mob, features = "nFeature_RNA", pt.size = 0, group.by = "sample_id")
dev.off()


## scatter plot - correlation between replicates
mob1 <- SubsetSTData(se_mob, sample_id == "Ex-ST_MOB1")
mob2 <- SubsetSTData(se_mob, sample_id == "Ex-ST_MOB2")



mob1_data <- GetAssayData(mob1, slot = "counts", assay = "RNA")
mob1_umi <- data.frame(rowSums(mob1_data), rowMeans(mob1_data > 0), row.names = rownames(mob1_data))
mob1_umi$log_counts_mob1 <- log1p(mob1_umi$rowSums.mob1_data.)

colnames(mob1_umi)[1] <- "rowSums_mob1"
colnames(mob1_umi)[2] <- "rowMeans_mob1"


mob2_data <- GetAssayData(mob2, slot = "counts", assay = "RNA")
mob2_umi <- data.frame(rowSums(mob2_data), rowMeans(mob2_data > 0), row.names = rownames(mob2_data))
mob2_umi$log_counts_mob2 <- log1p(mob2_umi$rowSums.mob2_data.)

colnames(mob2_umi)[1] <- "rowSums_mob2"
colnames(mob2_umi)[2] <- "rowMeans_mob2"


## merge 2 data frames by rownames
mob_final_df <- merge(mob1_umi, mob2_umi, by= 0, all=T, row.names = F)
rownames(mob_final_df) <- mob_final_df[,1]


ggplot(mob_final_df, aes(x=log_counts_mob1, y=log_counts_mob2)) + geom_point() +
      stat_cor(method = "pearson", label.x =  0.1, label.y = 10.5) +
      geom_smooth(method='lm', color='red') +
      ggtitle(label = "log1p(UMI counts)") +
      theme_bw() + theme_light()
