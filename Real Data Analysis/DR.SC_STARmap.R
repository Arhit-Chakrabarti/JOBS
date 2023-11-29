################################################################################
# Create the data from the STARMap dataset for LIGHT SAMPLE 1
################################################################################
data_list <- c("light_replicate_1")

# load gene count matrices
gene_data <-NULL
gene_data$data <- list()
gene_data$location <- list() # To compute the spatial locations

for (i in 1:length(data_list)){
  data <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".count.csv", sep = ""), header = TRUE)
  gene_data$data[[i]] <- t(data)
  loc <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".loc.csv", sep = ""), header = TRUE)
  gene_data$location[[i]] <- loc
}

gene_data$name <- read.csv(paste("./STARmap/", data_list[1],"/",data_list[1],".gene_symbol.csv", sep = ""), header = TRUE)[,1]


for(i in 1:length(gene_data$data)){
  rownames(gene_data$data[[i]]) <- gene_data$name
}

genes <- list()
l = 1
# GENES FOR LIGHT REPLICATE, CONSIST OF THE COMMON 33 GENES
genes[[l]] = c("Plcxd2",  "Cux2",    "Egr1",    "eRNA3",   "Slc17a7", "Mgp",     "Fos",
               "Prok2",   "Arx", "Nrn1",    "Nectin3", "Kcna1",   "Pcp4",    "Bgn",
               "Egr2",    "Rorb",    "Pthlh",   "Acss1",  "Tnfaip6", "Penk",    "Nr4a2",
               "Adgrl2",  "Pde1a",   "Deptor",  "Homer1",  "Otof",    "Bcl6",
               "Fam19a1", "Nov",     "Sema3e",  "Cpne5",   "Enpp2",   "Igtp")

counts <- gene_data$data[[l]][which(rownames(gene_data$data[[l]]) %in% genes[[l]]), ]
if(!require(Matrix))install.packages("Matrix"); suppressPackageStartupMessages(library(Matrix))
counts <- as(counts, "sparseMatrix")
colnames(counts) <- 1:ncol(counts)
if(!require(Seurat))install.packages("Seurat"); suppressPackageStartupMessages(library(Seurat))

meta_data <- data.frame(row = as.numeric(gene_data$location[[l]][,2]), 
                        col = as.numeric(gene_data$location[[l]][,1]),
                        annotation = 1:ncol(counts))

row.names(meta_data) <- 1:ncol(counts)
## create Seurat object
sce <- CreateSeuratObject(counts = counts, meta.data = meta_data)
# standard log-normalization
sce <- NormalizeData(sce, verbose = F)
seu <- FindVariableFeatures(sce, nfeatures = nrow(counts), verbose = F)
library(DR.SC)
# Perform Spatial clustering
seus <- DR.SC(seu, K = 2:20, platform = 'ST', verbose=T)

# Plot the DR.SC Clustering
plot.DRSC <- spatialPlotClusters(seus)
DR.SC.clusters <- plot.DRSC$data[,3]
plot.DRSC 

dataa = seu@assays$RNA@data

cluster_labels <- read.csv("./STARmap/light_replicate_1/light_replicate_1.cov.csv")

unique(cluster_labels$ClusterName)
library(tidyverse)
cluster_labels_exc <- cluster_labels %>% filter(ClusterName %in% c("eL2/3", "eL5", "eL6", "eL4") ) %>% dplyr::select(CellID, ClusterID, ClusterName)

est_labels <- 0

for(i in 1:length(cluster_labels_exc$CellID)){
  est_labels[i] <- DR.SC.clusters[which(as.numeric(rownames(gene_data$location[[1]])) %in% cluster_labels_exc$CellID[i])]
}

true_cluster <- as.factor(cluster_labels_exc$ClusterName)

library(sfheaders)
library(sf)
data_list<-c("light_replicate_1")

# load gene count matrices
gene_data <-NULL
gene_data$data <- list()
ncells <- rep(NA, length(data_list))
gene_data$location <- list() # To compute the spatial locations
gene_data$geometry <- list() # Load the geometry

for (i in 1:length(data_list)){
  data <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".count.csv", sep = ""), header = TRUE)
  gene_data$data[[i]] <- t(data)
  ncells[i] <- nrow(data)
}

# load gene names
gene_data$name <- read.csv(paste("./STARmap/", data_list[1],"/",data_list[1],".gene_symbol.csv", sep = ""), header = TRUE)[,1]

for (j in 1:length(data_list)){
  fold <- paste("./STARmap/", data_list[j], "/polygon/", sep = "")
  flst <- list.files(fold)
  all_locations <- NULL
  all_geometry <- NULL
  for (i in 1:length(flst)){
    fname <- paste(fold, "polygon", i, ".txt", sep = "")
    xym <- read.table(fname, header = F)
    yxm<-xym[,c("V2","V1")]
    p<-sfc_polygon(yxm)
    loc <- st_centroid(p)
    all_geometry <- rbind(all_geometry, p)
    all_locations <- rbind(all_locations, st_coordinates(loc))
    
  }
  gene_data$location[[j]] <- all_locations
  gene_data$geometry[[j]] <- all_geometry 
}

indexes = as.numeric(1:nrow(gene_data$location[[1]]))
locs_new = data.frame(V1 = gene_data$location[[1]][indexes, 2],
                      V2 = gene_data$location[[1]][indexes, 1])


Polygon_Data = st_sf(data.frame(all_geometry))
Data_new <- NULL

for(i in 1:length(indexes)){
  Data_new <-rbind(Data_new, Polygon_Data$all_geometry[indexes[i]])
}

Data_new2 = data.frame(Cluster = factor(DR.SC.clusters), Data_new)
Data_new2 = st_sf(data.frame(Data_new2))
library(tidyverse)

myvalues = c("1" = "#F8766D",
             "2" = "#00BA38",
             "3" = "#619CFF", 
             "4" = "blueviolet",
             "5" = "cyan4",
             "6" = "#E6AB02",
             "7" = "#E36EF6",
             "8" = "bisque4",
             "9" = "coral4",
             "10" = "darkslateblue")
library(cowplot)

plot_Light1 = ggplot() +
  geom_sf(data = Data_new2, aes(fill = Cluster), alpha = 0.8)  +
  theme_minimal() +
  labs(title = "", subtitle = paste0("ARI = ", round(aricode::ARI(true_cluster, est_labels), 4))) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        plot.title = element_text(face = "bold", size=14),
        plot.subtitle = element_text(face = "bold", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(face = "bold", size=18),
        legend.text = element_text(size = 16)
  ) + 
  guides(colour = guide_legend(title="clusters",
                               title.hjust = 0.5)) +
  scale_fill_manual(values = myvalues)

plot_Light1


################################################################################
# Create the data from the STARMap dataset for LIGHT SAMPLE 2
################################################################################
data_list <- c("light_replicate_2")

# load gene count matrices
gene_data <-NULL
gene_data$data <- list()
gene_data$location <- list() # To compute the spatial locations

for (i in 1:length(data_list)){
  data <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".count.csv", sep = ""), header = TRUE)
  gene_data$data[[i]] <- t(data)
  loc <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".loc.csv", sep = ""), header = TRUE)
  gene_data$location[[i]] <- loc
}

gene_data$name <- read.csv(paste("./STARmap/", data_list[1],"/",data_list[1],".gene_symbol.csv", sep = ""), header = TRUE)[,1]


for(i in 1:length(gene_data$data)){
  rownames(gene_data$data[[i]]) <- gene_data$name
}

genes <- list()
l = 1
# GENES FOR LIGHT REPLICATE, CONSIST OF THE COMMON 33 GENES
genes[[l]] = c("Plcxd2",  "Cux2",    "Egr1",    "eRNA3",   "Slc17a7", "Mgp",     "Fos",
               "Prok2",   "Arx", "Nrn1",    "Nectin3", "Kcna1",   "Pcp4",    "Bgn",
               "Egr2",    "Rorb",    "Pthlh",   "Acss1",  "Tnfaip6", "Penk",    "Nr4a2",
               "Adgrl2",  "Pde1a",   "Deptor",  "Homer1",  "Otof",    "Bcl6",
               "Fam19a1", "Nov",     "Sema3e",  "Cpne5",   "Enpp2",   "Igtp")

counts <- gene_data$data[[l]][which(rownames(gene_data$data[[l]]) %in% genes[[l]]), ]
if(!require(Matrix))install.packages("Matrix"); suppressPackageStartupMessages(library(Matrix))
counts <- as(counts, "sparseMatrix")
colnames(counts) <- 1:ncol(counts)
if(!require(Seurat))install.packages("Seurat"); suppressPackageStartupMessages(library(Seurat))

meta_data <- data.frame(row = as.numeric(gene_data$location[[l]][,2]), 
                        col = as.numeric(gene_data$location[[l]][,1]),
                        annotation = 1:ncol(counts))

row.names(meta_data) <- 1:ncol(counts)
## create Seurat object
sce <- CreateSeuratObject(counts = counts, meta.data = meta_data)
# standard log-normalization
sce <- NormalizeData(sce, verbose = F)
seu <- FindVariableFeatures(sce, nfeatures = nrow(counts), verbose = F)
library(DR.SC)
# Perform Spatial clustering
seus <- DR.SC(seu, K = 2:20, platform = 'ST', verbose=T)

# Plot the DR.SC Clustering
plot.DRSC <- spatialPlotClusters(seus)
DR.SC.clusters <- plot.DRSC$data[,3]
plot.DRSC 

dataa = seu@assays$RNA@data

cluster_labels <- read.csv("./STARmap/light_replicate_2/light_replicate_2.cov.csv")

unique(cluster_labels$ClusterName)
library(tidyverse)
cluster_labels_exc <- cluster_labels %>% filter(ClusterName %in% c("eL2/3", "eL5", "eL6", "eL4") ) %>% dplyr::select(CellID, ClusterID, ClusterName)

est_labels <- 0

for(i in 1:length(cluster_labels_exc$CellID)){
  est_labels[i] <- DR.SC.clusters[which(as.numeric(rownames(gene_data$location[[1]])) %in% cluster_labels_exc$CellID[i])]
}

true_cluster <- as.factor(cluster_labels_exc$ClusterName)

library(sfheaders)
library(sf)
data_list<-c("light_replicate_2")

# load gene count matrices
gene_data <-NULL
gene_data$data <- list()
ncells <- rep(NA, length(data_list))
gene_data$location <- list() # To compute the spatial locations
gene_data$geometry <- list() # Load the geometry

for (i in 1:length(data_list)){
  data <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".count.csv", sep = ""), header = TRUE)
  gene_data$data[[i]] <- t(data)
  ncells[i] <- nrow(data)
}

# load gene names
gene_data$name <- read.csv(paste("./STARmap/", data_list[1],"/",data_list[1],".gene_symbol.csv", sep = ""), header = TRUE)[,1]

for (j in 1:length(data_list)){
  fold <- paste("./STARmap/", data_list[j], "/polygon/", sep = "")
  flst <- list.files(fold)
  all_locations <- NULL
  all_geometry <- NULL
  for (i in 1:length(flst)){
    fname <- paste(fold, "polygon", i, ".txt", sep = "")
    xym <- read.table(fname, header = F)
    yxm<-xym[,c("V2","V1")]
    p<-sfc_polygon(yxm)
    loc <- st_centroid(p)
    all_geometry <- rbind(all_geometry, p)
    all_locations <- rbind(all_locations, st_coordinates(loc))
    
  }
  gene_data$location[[j]] <- all_locations
  gene_data$geometry[[j]] <- all_geometry 
}

indexes = as.numeric(1:nrow(gene_data$location[[1]]))
locs_new = data.frame(V1 = gene_data$location[[1]][indexes, 2],
                      V2 = gene_data$location[[1]][indexes, 1])


Polygon_Data = st_sf(data.frame(all_geometry))
Data_new <- NULL

for(i in 1:length(indexes)){
  Data_new <-rbind(Data_new, Polygon_Data$all_geometry[indexes[i]])
}


Data_new2 = data.frame(Cluster = factor(DR.SC.clusters), Data_new)
Data_new2 = st_sf(data.frame(Data_new2))
library(tidyverse)

myvalues = c("1" = "#F8766D",
             "2" = "#00BA38",
             "3" = "#619CFF", 
             "4" = "blueviolet",
             "5" = "cyan4",
             "6" = "#E6AB02",
             "7" = "#E36EF6",
             "8" = "bisque4",
             "9" = "coral4",
             "10" = "darkslateblue")
library(cowplot)

plot_Light2 = ggplot() +
  geom_sf(data = Data_new2, aes(fill = Cluster), alpha = 0.8)  +
  theme_minimal() +
  labs(title = "", subtitle = paste0("ARI = ", round(aricode::ARI(true_cluster, est_labels), 4))) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        plot.title = element_text(face = "bold", size=14),
        plot.subtitle = element_text(face = "bold", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(face = "bold", size=18),
        legend.text = element_text(size = 16)
  ) + 
  guides(colour = guide_legend(title="clusters",
                               title.hjust = 0.5)) +
  scale_fill_manual(values = myvalues)

plot_Light2

