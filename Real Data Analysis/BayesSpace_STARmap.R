################################################################################
# Create the data from the STARMap dataset for LIGHT SAMPLE 1
################################################################################
data_list <- c("light_replicate_1")

# load gene count matrices
gene_data <-NULL
gene_data$data <- list()
ncells <- rep(NA, length(data_list))
gene_data$location <- list() # To compute the spatial locations

for (i in 1:length(data_list)){
  data <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".count.csv", sep = ""), header = TRUE)
  gene_data$data[[i]] <- t(data)
  ncells[i] <- nrow(data)
  loc <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".loc.csv", sep = ""), header = TRUE)
  gene_data$location[[i]] <- loc
}

gene_data$name <- read.csv(paste("./STARmap/", data_list[1],"/",data_list[1],".gene_symbol.csv", sep = ""), header = TRUE)[,1]


library(BayesSpace)
colnames(loc) <- c("row", "col")
rowData <- gene_data$name
colData <- loc
counts <- gene_data$data[[i]]

sce <- SingleCellExperiment(assays=list(counts=as(counts, "dgCMatrix")),
                            rowData=rowData,
                            colData=colData)

sce2 <- sce[, colSums(counts(sce)) > 0]

sce2 <- spatialPreprocess(sce2, platform="ST", 
                          n.PCs=7, n.HVGs=50, log.normalize = TRUE,
                          assay.type = "counts")

STARmap <- spatialCluster(sce2, q=4, platform="ST", d=7,
                          init.method="mclust", model="t", gamma=2,
                          nrep=2000, burn.in=1000,
                          save.chain=TRUE)
library(tidyverse)

cluster_labels <- read.csv("./STARmap/light_replicate_1/light_replicate_1.cov.csv")

cluster_labels_exc <- cluster_labels %>% filter(ClusterName %in% c("eL2/3", "eL5", "eL6", "eL4") ) %>% dplyr::select(CellID, ClusterID, ClusterName)

indexes <- which(colSums(counts(sce)) > 0)
cluster_labels_exc <- cluster_labels_exc[cluster_labels_exc$CellID %in% intersect(indexes, cluster_labels_exc$CellID), ]
est_labels <- 0

for(i in 1:length(cluster_labels_exc$CellID)){
  est_labels[i] <- STARmap$spatial.cluster[which(indexes %in% cluster_labels_exc$CellID[i])]
}

true_cluster <- as.factor(cluster_labels_exc$ClusterName)

aricode::ARI(true_cluster, est_labels)

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

locs_new = data.frame(V1 = gene_data$location[[1]][indexes, 2],
                      V2 = gene_data$location[[1]][indexes, 1])

Polygon_Data = st_sf(data.frame(all_geometry))
Data_new <- NULL

for(i in 1:length(indexes)){
  Data_new <-rbind(Data_new, Polygon_Data$all_geometry[indexes[i]])
}

Data_new2 = data.frame(Cluster = factor(STARmap$spatial.cluster), Data_new)
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
ncells <- rep(NA, length(data_list))
gene_data$location <- list() # To compute the spatial locations

for (i in 1:length(data_list)){
  data <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".count.csv", sep = ""), header = TRUE)
  gene_data$data[[i]] <- t(data)
  ncells[i] <- nrow(data)
  loc <- read.csv(paste("./STARmap/", data_list[i],"/",data_list[i],".loc.csv", sep = ""), header = TRUE)
  gene_data$location[[i]] <- loc
}

gene_data$name <- read.csv(paste("./STARmap/", data_list[1],"/",data_list[1],".gene_symbol.csv", sep = ""), header = TRUE)[,1]


library(BayesSpace)
colnames(loc) <- c("row", "col")
rowData <- gene_data$name
colData <- loc
counts <- gene_data$data[[i]]

sce <- SingleCellExperiment(assays=list(counts=as(counts, "dgCMatrix")),
                            rowData=rowData,
                            colData=colData)

sce2 <- sce[, colSums(counts(sce)) > 0]

sce2 <- spatialPreprocess(sce2, platform="ST", 
                          n.PCs=7, n.HVGs=50, log.normalize = TRUE,
                          assay.type = "counts")

STARmap <- spatialCluster(sce2, q=4, platform="ST", d=7,
                          init.method="mclust", model="t", gamma=2,
                          nrep=2000, burn.in=1000,
                          save.chain=TRUE)
library(tidyverse)

cluster_labels <- read.csv("./STARmap/light_replicate_2/light_replicate_2.cov.csv")

cluster_labels_exc <- cluster_labels %>% filter(ClusterName %in% c("eL2/3", "eL5", "eL6", "eL4") ) %>% dplyr::select(CellID, ClusterID, ClusterName)

indexes <- which(colSums(counts(sce)) > 0)
cluster_labels_exc <- cluster_labels_exc[cluster_labels_exc$CellID %in% intersect(indexes, cluster_labels_exc$CellID), ]
est_labels <- 0

for(i in 1:length(cluster_labels_exc$CellID)){
  est_labels[i] <- STARmap$spatial.cluster[which(indexes %in% cluster_labels_exc$CellID[i])]
}

true_cluster <- as.factor(cluster_labels_exc$ClusterName)

aricode::ARI(true_cluster, est_labels)

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

locs_new = data.frame(V1 = gene_data$location[[1]][indexes, 2],
                      V2 = gene_data$location[[1]][indexes, 1])

Polygon_Data = st_sf(data.frame(all_geometry))
Data_new <- NULL

for(i in 1:length(indexes)){
  Data_new <-rbind(Data_new, Polygon_Data$all_geometry[indexes[i]])
}

Data_new2 = data.frame(Cluster = factor(STARmap$spatial.cluster), Data_new)
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
