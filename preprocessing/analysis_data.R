######################################
# Install and require needed package #
######################################
install.packages("ggplot2")
install.packages("Matrix")
install.packages("rjson")
install.packages("cowplot")
install.packages("RColorBrewer")
install.packages("Seurat")
install.packages("grid")
install.packages("readbitmap")
install.packages("dplyr")
install.packages("data.table")
install.packages("doSNOW")
install.packages("hdf5r")
install.packages('remotes')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("rhdf5")

require(remotes)
remotes::install_github("satijalab/seurat-data")
require(SeuratData)
require(patchwork)
require(ggplot2)
require(pbapply)
require(Seurat)
require(stringr)
require(ggplot2)
require(Matrix)
require(rjson)
require(cowplot)
require(RColorBrewer)
require(grid)
require(readbitmap)
require(Seurat)
require(dplyr)
require(hdf5r)
require(data.table)
require(doSNOW)
require(tidyr)
require(tibble)
require(gridExtra)




#############################
# Load data for four sample #
#############################
corrected.sample1.data <- Read10X(
  data.dir = "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr1/outs/filtered_feature_bc_matrix/"
)

corrected.sample1 <- CreateSeuratObject(counts = corrected.sample1.data,
                                        project = "corrected-S3647Nr1")
corrected.sample1 <- AddMetaData(corrected.sample1, metadata="corrected-S3647Nr1", col.name = "sample")


corrected.sample2.data <- Read10X(
  data.dir = "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr2/outs/filtered_feature_bc_matrix/"
)

corrected.sample2 <- CreateSeuratObject(counts = corrected.sample2.data,
                                        project = "corrected-S3647Nr2")
corrected.sample2 <- AddMetaData(corrected.sample2, metadata="corrected-S3647Nr2", col.name = "sample")


corrected.sample3.data <- Read10X(
  data.dir = "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr3/outs/filtered_feature_bc_matrix/"
)

corrected.sample3 <- CreateSeuratObject(counts = corrected.sample3.data,
                                        project = "corrected-S3647Nr3")
corrected.sample3 <- AddMetaData(corrected.sample3, metadata="corrected-S3647Nr3", col.name = "sample")


corrected.sample4.data <- Read10X(
  data.dir = "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr4/outs/filtered_feature_bc_matrix/"
)

corrected.sample4 <- CreateSeuratObject(counts = corrected.sample4.data,
                                        project = "corrected-S3647Nr4")
corrected.sample4 <- AddMetaData(corrected.sample4, metadata="corrected-S3647Nr4", col.name = "sample")


###################
# Merging samples #
###################
merged.samples <- merge(
  corrected.sample1,
  y = c(corrected.sample2, corrected.sample3, corrected.sample4),
  add.cell.ids = c("corrected-S3647Nr1", "corrected-S3647Nr2",
                   "corrected-S3647Nr3", "corrected-S3647Nr4"),
  project = "merged.sample"
)


####################
# integration data #
####################
# split the dataset into a list
merged.samples.list <- SplitObject(merged.samples, split.by = "sample")

# normalize and identify variable features for each dataset independently
merged.samples.list <- lapply(X = merged.samples.list,
                              FUN = function(x) {
                                x <- NormalizeData(x, normalization.method = "LogNormalize")
                                x <- FindVariableFeatures(x, 
                                                          selection.method = "vst",
                                                          nfeatures = 2000)
                              })


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = merged.samples.list, nfeatures = 2000)


merged.anchors <- FindIntegrationAnchors(object.list = merged.samples.list,
                                         anchor.features = features)

# this command creates an 'integrated' data assay
integrated.data <- IntegrateData(anchorset = merged.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(integrated.data) <- "integrated"


integrated.analysis <- ScaleData(integrated.data, verbose = FALSE)
integrated.analysis <- RunPCA(integrated.analysis, npcs = 30, verbose = FALSE)
integrated.analysis <- RunUMAP(integrated.analysis, reduction = "pca", dims = 1:30)
integrated.analysis <- FindNeighbors(integrated.analysis, reduction = "pca", dims = 1:30)



integrated.analysis.cluster <- FindClusters(integrated.analysis, resolution = 2.5)
# integrated.analysis.cluster <- RunTSNE(integrated.analysis.cluster)


#########################
# Load data about peaks #
#########################
info.peaks <- read.table("/projects/ifpan-annaradli-spatialtrial/data/peaks_annotate_sorted.bed", 
                         header = FALSE,
                         sep = "\t",
                         col.names = c('chromosome',
                                       'start_peak',
                                       'end_peak',
                                       'peak_id',
                                       'score_int(-10*log10pvalue)',
                                       'strand_coverage',
                                       'fold_change_peak_summit',
                                       '-log10pvalue_peak_summit',
                                       '-log10qvalue_peak_summit',
                                       'relative_summit_position_peak_start',
                                       'type_peak',
                                       'chromosome2',
                                       'start_gene',
                                       'end_gene',
                                       'gene_id',
                                       'gene_name',
                                       'strand_gene'))

# For all features
sample.data.all <- integrated.analysis@assays$RNA@counts %>% as.matrix()

# sample.data.all.normalize <- integrated.analysis.cluster@assays$RNA@data %>% as.matrix()

sample.anno.all <- info.peaks[match(str_replace_all(rownames(sample.data.all), "-", "_"), 
                                    str_replace_all(info.peaks$peak_id, "-", "_")),]

sample.info.all <- integrated.analysis@meta.data %>%
  mutate(treatment = ifelse(sample == "corrected-S3647Nr1" | 
                              sample == "corrected-S3647Nr2", "treat", "ctrl"))


# prepare cluster with resulution between 0.1 to 3
for (resolution in seq(0.1, 3, length.out=30)) {
  sample.info.all %>%
    mutate(!!as.name({paste("integrated_snn_res.", 
                            as.character(resolution), 
                            sep = "")}) := FindClusters(integrated.analysis, 
                                                        resolution = resolution)@meta.data[, 5]) -> sample.info.all
  print(resolution)
}


rm(merged.samples,
   merged.samples.list,
   merged.anchors,
   integrated.data)



######################
# Normalization data #
######################
results <- data.frame(median = pbapply(sample.data.all, 1, median))
results$mean <- pbapply(sample.data.all, 1, mean)

sample.anno.all$strand_agree <- as.numeric(paste(sample.anno.all$strand_coverage, "1", sep = "")) == sample.anno.all$strand_gene
filtered.wh <- which(sample.anno.all$strand_agree & results$mean > 0.05)
filtered.anno <- sample.anno.all[filtered.wh,]
filtered.data <- sample.data.all[filtered.wh,]
filtered.info <- sample.info.all

filtered.info$over_2 <- apply(filtered.data, 2, function(x){sum(x > 2)})

threshold <- 500

colfilt.wh <- which(filtered.info$over_2 > threshold)
colfilt.anno <- filtered.anno
colfilt.data <- filtered.data[,colfilt.wh]
colfilt.info <- filtered.info[colfilt.wh,]

normalize <- function(x, range = 1000) {
  order <- order(x, decreasing = T)
  out <- rep(0, length(x))
  out[order[1:1000]] <- 1000:1
  out
}

colfilt.norm.data <- apply(colfilt.data, 1, normalize, threshold) %>% t

stat <- function(x, treat = colfilt.info$treatment, clusters = as.factor(colfilt.info$integrated_snn_res.2)) {
  idx = 1
  p <- vector()
  for(cluster in levels(clusters)) {
    p[idx] <- 1
    if (sum(clusters == cluster & treat == "ctrl") < 2 |
        sum(clusters == cluster & treat == "treat") < 2) {
      
    } else {
      p[idx] <- wilcox.test(x[clusters == cluster & treat == "ctrl"], x[clusters == cluster & treat == "treat"])$p.value
    }
    
    idx = idx + 1
  }
  names(p) <- levels(clusters)
  p
}
colfilt.stat <- pbapply(colfilt.norm.data, 1, stat) %>% t


# Add column names to colfilt.norm.data
colnames(colfilt.norm.data) <- rownames(colfilt.info)

#######################
# Execute correlation # 
#######################
colfilt.norm.data %>% t %>% cor -> colfilt.norm.cor 


cor_interest_gene <- function(gene, n = 10) {
  colfilt.anno %>% 
    select(peak_id, gene_name) %>% 
    mutate(peak_id = str_replace_all(peak_id, "_", "-")) %>% 
    filter(gene_name == gene) %>%
    .[,1] -> tmp.peaks
  
  for (peak in tmp.peaks) {
    colfilt.norm.cor[, peak] %>% 
      as.data.frame() %>%
      rownames_to_column(var = "peak_id") %>%
      rename(peak = ".") %>% 
      .[order(.$peak, decreasing = TRUE),] %>%
      rename(!!as.name(peak) := peak) %>%
      head(n) %>% 
      left_join(., {colfilt.anno %>% 
          select(peak_id, gene_name) %>% 
          mutate(peak_id = str_replace_all(peak_id, "_", "-"))}, by = "peak_id") %>%
      print()
  }
}

cor_interest_gene(gene = "Fzd2", n = 10) 


colfilt.info.peaks <- info.peaks %>%
  filter(peak_id %in% {rownames(colfilt.norm.cor) %>% 
      str_replace_all(., "-", "_")}) 



##################################
# Prepare function: geom_spatial #
##################################
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

palette_cluster <- c("0" = "#b2df8a",
                     "1" = "#e41a1c",
                     "2" = "#377eb8",
                     "3" ="#4daf4a",
                     "4" = "#ff7f00",
                     "5" = "gold", 
                     "6" = "#a65628", 
                     "7" = "#999999", 
                     "8" = "black", 
                     "9" = "grey", 
                     "10" = "white", 
                     "11" = "purple",
                     "12" = "red", 
                     "13" = "blue", 
                     "14" = "pink",
                     "15" = "brown",
                     "16" = "green",
                     "17" = "tomato1",
                     "18" = "yellow3",
                     "19" = "violet",
                     "20" = "yellowgreen",
                     "21" = "lightblue1",
                     "22" = "lightblue4",
                     "23" = "lightgoldenrod3",
                     "24" = "lightpink2",
                     "25" = "magenta",
                     "26" = "limegreen",
                     "27" = "maroon",
                     "28" = "mintcream",
                     "29"=  "oldlace",
                     "30" = "tan",
                     "31" = "chartreuse",
                     "32" = "blue4",
                     "33" = "midnightblue",
                     "34" = "slategray4",
                     "35" = "snow3",
                     "36" = "springgreen",
                     "37" = "plum1")

geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = FALSE,
                          ...) {
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


###################
# Define my paths #
###################
sample_names <- c("corrected-S3647Nr1",  "corrected-S3647Nr2",  "corrected-S3647Nr3",  "corrected-S3647Nr4")

image_paths <- c("/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr1/outs/spatial/tissue_lowres_image.png",
                 "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr2/outs/spatial/tissue_lowres_image.png",
                 "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr3/outs/spatial/tissue_lowres_image.png",
                 "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr4/outs/spatial/tissue_lowres_image.png")


scalefactor_paths <- c("/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr1/outs/spatial/scalefactors_json.json",
                       "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr2/outs/spatial/scalefactors_json.json",
                       "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr3/outs/spatial/scalefactors_json.json",
                       "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr4/outs/spatial/scalefactors_json.json")


tissue_paths <- c("/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr1/outs/spatial/tissue_positions_list.csv",
                  "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr2/outs/spatial/tissue_positions_list.csv",
                  "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr3/outs/spatial/tissue_positions_list.csv",
                  "/projects/ifpan-annaradli-spatialtrial/data/spaceranger_results/corrected_S3647Nr4/outs/spatial/tissue_positions_list.csv")


# cluster_paths <- c("/projects/ifpan-annaradli-spatialtrial/tmp/data/spaceranger_results/corrected_S3647Nr1/outs/analysis/clustering/graphclust/clusters.csv",
#                    "/projects/ifpan-annaradli-spatialtrial/tmp/data/spaceranger_results/corrected_S3647Nr2/outs/analysis/clustering/graphclust/clusters.csv",
#                    "/projects/ifpan-annaradli-spatialtrial/tmp/data/spaceranger_results/corrected_S3647Nr3/outs/analysis/clustering/graphclust/clusters.csv",
#                    "/projects/ifpan-annaradli-spatialtrial/tmp/data/spaceranger_results/corrected_S3647Nr4/outs/analysis/clustering/graphclust/clusters.csv")
# 
# 
# matrix_paths <- c("/projects/ifpan-annaradli-spatialtrial/tmp/data/spaceranger_results/corrected_S3647Nr1/outs/filtered_feature_bc_matrix.h5",
#                   "/projects/ifpan-annaradli-spatialtrial/tmp/data/spaceranger_results/corrected_S3647Nr2/outs/filtered_feature_bc_matrix.h5",
#                   "/projects/ifpan-annaradli-spatialtrial/tmp/data/spaceranger_results/corrected_S3647Nr3/outs/filtered_feature_bc_matrix.h5",
#                   "/projects/ifpan-annaradli-spatialtrial/tmp/data/spaceranger_results/corrected_S3647Nr4/outs/filtered_feature_bc_matrix.h5")



###############################
# Read in down sampled images #
###############################
images_cl <- list()

for (i in 1:length(sample_names)) {
  images_cl[[i]] <- read.bitmap(image_paths[i])
}

height <- list()

for (i in 1:length(sample_names)) {
  height[[i]] <-  data.frame(height = nrow(images_cl[[i]]))
}

height <- bind_rows(height)

width <- list()

for (i in 1:length(sample_names)) {
  width[[i]] <- data.frame(width = ncol(images_cl[[i]]))
}

width <- bind_rows(width)


###############################
# Convert the Images to grobs #
###############################
grobs <- list()
for (i in 1:length(sample_names)) {
  grobs[[i]] <- rasterGrob(images_cl[[i]], width=unit(1,"npc"), height=unit(1,"npc"))
}

images_tibble <- tibble(sample=factor(sample_names), grob=grobs)
images_tibble$height <- height$height
images_tibble$width <- width$width
scales <- list()

for (i in 1:length(sample_names)) {
  scales[[i]] <- rjson::fromJSON(file = scalefactor_paths[i])
}




# ####################
# # Read in clusters #
# ####################
# clusters <- list()
# for (i in 1:length(sample_names)) {
#   clusters[[i]] <- read.csv(cluster_paths[i])
# }


#############################################################
# Combine clusters and tissue information for easy plotting #
#############################################################
bcs <- list()

for (i in 1:length(sample_names)) {
  bcs[[i]] <- read.csv(tissue_paths[i],col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
  bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_lowres_scalef    # scale tissue coordinates for lowres image
  bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_lowres_scalef
  bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
  bcs[[i]]$height <- height$height[i]
  bcs[[i]]$width <- width$width[i]
}

names(bcs) <- sample_names



################################
# Merge All the Necessary Data #
################################
bcs_merge <- bind_rows(bcs, .id = "sample")
bcs_merge <- merge(bcs_merge,umi_sum, by = c("barcode", "sample"))
bcs_merge <- merge(bcs_merge,gene_sum, by = c("barcode", "sample"))


###########################
# Function to create plot #
###########################

###############################################
# Cluster Assignments per Tissue Covered Spot #
###############################################
plot_interest_cluster <- function(data_cluster, size = 1.2, interest_cluster){
  plots <- list()
  
  for (i in 1:length(sample_names)) {
    
    plots[[i]] <- bcs_merge %>% 
      left_join(., {data_cluster %>% 
          .$seurat_clusters %>% 
          as.data.frame() %>% 
          rename(cluster = ".") %>% 
          rownames_to_column(var = "sample.barcode") %>% 
          separate("sample.barcode", c("sample", "barcode"), "_")}, by = c("sample", "barcode")) %>%
      mutate(interest_cluster = ifelse(cluster == interest_cluster, 1, 0)) %>%
      filter(sample ==sample_names[i]) %>%
      filter(tissue == "1") %>% 
      ggplot(aes(x=imagecol,y=imagerow,fill=factor(interest_cluster))) +
      geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape = 21, colour = "black", size = size, stroke = 0.25)+
      coord_cartesian(expand=FALSE) +
      scale_fill_manual(values = c("gray45", "red"))+
      xlim(0,max(bcs_merge %>% 
                   filter(sample ==sample_names[i]) %>% 
                   select(width)))+
      ylim(max(bcs_merge %>% 
                 filter(sample ==sample_names[i]) %>% 
                 select(height)),0)+
      xlab("") +
      ylab("") +
      ggtitle(sample_names[i])+
      labs(fill = "Cluster")+
      guides(fill = guide_legend(override.aes = list(size=3)))+
      theme_set(theme_bw(base_size = 10))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank()) +
      NoLegend()
  }
  
  plot_grid(plotlist = plots)
  
}


plot_clusters <- function(data_cluster, size = 1.2){
  plots <- list()
  
  for (i in 1:length(sample_names)) {
    
    plots[[i]] <- bcs_merge %>% 
      left_join(., {data_cluster %>% 
          .$seurat_clusters %>% 
          as.data.frame() %>% 
          rename(cluster = ".") %>% 
          rownames_to_column(var = "sample.barcode") %>% 
          separate("sample.barcode", c("sample", "barcode"), "_")}, by = c("sample", "barcode")) %>%
      filter(sample ==sample_names[i]) %>%
      filter(tissue == "1") %>% 
      ggplot(aes(x=imagecol,y=imagerow,fill=factor(cluster))) +
      geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape = 21, colour = "black", size = size, stroke = 0.25)+
      coord_cartesian(expand=FALSE) +
      scale_fill_manual(values = palette_cluster)+
      xlim(0,max(bcs_merge %>% 
                   filter(sample ==sample_names[i]) %>% 
                   select(width)))+
      ylim(max(bcs_merge %>% 
                 filter(sample ==sample_names[i]) %>% 
                 select(height)),0)+
      xlab("") +
      ylab("") +
      ggtitle(sample_names[i])+
      labs(fill = "Cluster")+
      guides(fill = guide_legend(override.aes = list(size=3)))+
      theme_set(theme_bw(base_size = 10))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank()) +
      NoLegend()
  }
  
  plot_grid(plotlist = plots)
  
}


plot_feature <- function(data_cluster, peak_id, size){
  plots <- list()
  
  normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
  }
  
  for (i in 1:length(sample_names)) {
    
    plots[[i]] <- colfilt.norm.data %>% 
      .[peak_id,] %>%
      as.data.frame() %>% 
      rename(value = ".") %>% 
      rownames_to_column(var = "sample.barcode") %>% 
      separate("sample.barcode", c("sample", "barcode"), sep = "_") %>% 
      left_join(., data_cluster, by = c("barcode", "sample")) %>% 
      mutate(percentile0.99 = as.numeric(quantile(value, probs = c(0.99))),
             percentile0.05 = as.numeric(quantile(value, probs = c(0.05)))) %>% 
      mutate(value = ifelse(value < percentile0.99, value, percentile0.99),
             value = ifelse(value > percentile0.05, value, percentile0.05)) %>%
      mutate(value = normalize(value)) %>%
      filter(sample == sample_names[i]) %>% 
      ggplot(aes(x=imagecol,y=imagerow,fill= value)) +
      geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape = 21, colour = "black", size = size, stroke = 0.25)+
      coord_cartesian(expand=FALSE)+
      scale_fill_gradientn(colours = myPalette(100), limits=c(0,1))+
      xlim(0,max(bcs_merge %>% 
                   filter(sample ==sample_names[i]) %>% 
                   select(width)))+
      ylim(max(bcs_merge %>% 
                 filter(sample ==sample_names[i]) %>% 
                 select(height)),0)+
      xlab("") +
      ylab("") +
      ggtitle(sample_names[i])+
      theme_set(theme_bw(base_size = 10))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank())
  }
  
  plot_grid(plotlist = plots)
  
}


plot_feature_cluster <- function(data_cluster, peak_id, interest_cluster, size){
  plots <- list()
  
  normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
  }
  
  for (i in 1:length(sample_names)) {
    
    plots[[i]] <- colfilt.norm.data %>% 
      .[peak_id,] %>%
      as.data.frame() %>% 
      rename(value = ".") %>% 
      rownames_to_column(var = "sample.barcode") %>% 
      separate("sample.barcode", c("sample", "barcode"), sep = "_") %>% 
      left_join(., data_cluster, by = c("barcode", "sample")) %>% 
      mutate(percentile0.99 = as.numeric(quantile(value, probs = c(0.99))),
             percentile0.05 = as.numeric(quantile(value, probs = c(0.05)))) %>% 
      mutate(value = ifelse(value < percentile0.99, value, percentile0.99),
             value = ifelse(value > percentile0.05, value, percentile0.05)) %>%
      mutate(value = normalize(value)) %>%
      filter(sample == sample_names[i]) %>% 
      filter(cluster == interest_cluster) %>% 
      ggplot(aes(x=imagecol,y=imagerow,fill= value)) +
      geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape = 21, colour = "black", size = size, stroke = 0.25)+
      coord_cartesian(expand=FALSE)+
      scale_fill_gradientn(colours = myPalette(100), limits=c(0,1))+
      xlim(0,max(bcs_merge %>% 
                   filter(sample ==sample_names[i]) %>% 
                   select(width)))+
      ylim(max(bcs_merge %>% 
                 filter(sample ==sample_names[i]) %>% 
                 select(height)),0)+
      xlab("") +
      ylab("") +
      ggtitle(sample_names[i])+
      theme_set(theme_bw(base_size = 10))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank())
  }
  
  plot_grid(plotlist = plots)
  
}
