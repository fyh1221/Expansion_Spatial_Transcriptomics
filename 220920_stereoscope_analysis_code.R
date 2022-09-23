# Mouse brain samples (MOB + Hippocampus, standard Visium and Ex-ST) Stereoscope Analysis

## Load packages
library(STutility)
library(ggplot2)
library(magrittr)

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



## Set paths for stereoscope output directories
stsc.dir <- "/home/st-analysis_home/zaneta.andrusivova/projects/expST/220629_omni_mouse_brain/stsc/results/ST_data_both_assays/"

## Path for where the plots will be exported to
export_path = "/home/st-analysis_home/zaneta.andrusivova/projects/expST/220629_omni_mouse_brain/stsc/"


## Import, format and add sterescope output as an assay to Seurat object
stsc.files <- list.files(stsc.dir, recursive = TRUE, full.names = TRUE, pattern = '^W')
stsc.output <- lapply(stsc.files, read.delim, row.names = 1)
stsc.output <- lapply(stsc.output, t)

ids <- sapply(stsc.output, function(x) {
  unlist(strsplit(colnames(x)[1], "_"))[2]
})

stsc.files <- setNames(stsc.files, nm = ids)
stsc.files <- stsc.files[paste0(1:24)]


## Make one matrix out of the list of different sections
stsc.output2 <- do.call(cbind, stsc.output)

## Create empty matrix
emptyMat <- matrix(0, nrow = nrow(stsc.output2), ncol = ncol(se))
colnames(emptyMat) <- colnames(se)
rownames(emptyMat) <- rownames(stsc.output2)

## Add stereoscope results
emptyMat[, intersect(colnames(se), colnames(stsc.output2))] <- stsc.output2[, intersect(colnames(se), colnames(stsc.output2))]

## Create assay object
stereoscope <- CreateAssayObject(data = emptyMat)

## Add stereoscope as an assay to Seurat object
se[["stereoscope"]] <- stereoscope

## Set stereoscope assay as default
DefaultAssay(se) <- "stereoscope"


## Generate plots for different cell types
pdf( "Excitatory_neurons_hippocampus_CA3.pdf",  width = 15, height = 11, bg = NA)
ST.FeaturePlot(se, features = "Excitatory.neurons..hippocampus.CA3", ncol = 3, pt.size = 2, pt.border = F,
               cols = c("light blue", "light yellow", "orange", "red", "dark red"))
dev.off()





# Density plot, maximum cell proportion value per spot

## Export metadat from Seurat object
metadata <- se@meta.data

## Create data frame with maximum value of cell proportion for each spot
stsc_t <- t(stsc.output2)

df <- data.frame(stsc_t)
df$max <- apply(df, 1, max)

df_filter <- df[rownames(df), df$max]
df_filter$max <- df$max[match(rownames(df_filter), rownames(df))]

## Merge data frame with metadata exported from Seurat object, round max cell type proportion column
final_df <- merge(metadata, df_filter, by = 0)
final_df$roundmax <- round(final_df$max ,digit=2)
final_df$percent_max <- final_df$roundmax*100

## Subset the final data frame for one sample type (MOB or Hippocampus)
final_df_MOB <- subset(final_df, final_df$anat == "MOB")

## Set colors to use for plotting
my_colors <- c("#F8DFDA", "#D8E9F5")

## Density plot
final_df_MOB %>% 
  ggplot(aes(color=assay, x= max, fill = assay)) + 
  geom_density(alpha = 0.2) + 
  xlab("Max % of cell type") +
  ylab("Density") +
  theme_classic() +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values = my_colors)

