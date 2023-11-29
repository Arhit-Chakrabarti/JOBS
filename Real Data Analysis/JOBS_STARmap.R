data_list <- c("light_replicate_1",
               "light_replicate_2")

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


gene_data$reads_per_cell <-list()
median_expression <- 0

for(i in 1:length(gene_data$data)){
  ind <- which(apply(gene_data$data[[i]], 2, sum) > 100 & apply(gene_data$data[[i]], 2, sum) < 2000)
  gene_data$data[[i]] <- gene_data$data[[i]][ , ind]
  gene_data$location[[i]] <- gene_data$location[[i]][ind, ]
  gene_data$reads_per_cell[[i]] <- apply(gene_data$data[[i]], 2, sum)
  median_expression[i] <- median(gene_data$reads_per_cell[[i]])
  rownames(gene_data$data[[i]]) <- gene_data$name
}


gene_data$N <-  replicate(length(gene_data$data), "list", simplify = FALSE)

for(l in 1:length(gene_data$data)){
  gene_data$N[[l]] <- matrix(NA, nrow = dim(gene_data$data[[l]])[1], ncol = dim(gene_data$data[[l]])[2])
  for (i in 1:ncol(gene_data$N[[l]])){
    # the formulae used in STARmap protocol
    gene_data$N[[l]][,i] <- log(1 + median_expression[l]*((gene_data$data[[l]][,i] + 0.01)/sum(gene_data$data[[l]][,i])))
  }
  rownames(gene_data$N[[l]]) <- gene_data$name
}

genes <- list()

for(l in 1:length(gene_data$data)){
  counts <- gene_data$data[[l]]
  if(!require(Matrix))install.packages("Matrix"); suppressPackageStartupMessages(library(Matrix))
  counts <- as(counts, "sparseMatrix")
  colnames(counts) <- 1:ncol(counts)
  if(!require(Seurat))install.packages("Seurat"); suppressPackageStartupMessages(library(Seurat))
  
  meta_data <- data.frame(row = as.numeric(gene_data$location[[l]][,1]), 
                          col = as.numeric(gene_data$location[[l]][,2]),
                          annotation = 1:ncol(counts))
  
  row.names(meta_data) <- 1:ncol(counts)
  ## create Seurat object
  sce <- CreateSeuratObject(counts = counts, meta.data = meta_data)
  # standard log-normalization
  sce <- NormalizeData(sce, verbose = F)
  seu <- FindVariableFeatures(sce, nfeatures = 160, verbose = F)
  if(!require(DR.SC))install.packages("DR.SC"); suppressPackageStartupMessages(library(DR.SC))
  seus <- FindSVGs(seu, nfeatures = 160, verbose = F)
  genes[[l]] = seus@assays$RNA@var.features[1:50]
}

common.genes = intersect(genes[[1]], genes[[2]])

dataa <- list()
locs <- list()
nearest_neighbors <- list()

for(l in 1:length(gene_data$data)){
  dataa[[l]] <- gene_data$N[[l]][common.genes, ]
  if(!require(fields))install.packages("fields"); suppressPackageStartupMessages(library(fields))
  if(!require(NPVecchia))install.packages("NPVecchia"); suppressPackageStartupMessages(library(NPVecchia))
  
  order <- orderMaxMinFaster(gene_data$location[[l]])
  
  # Reorder data and location by maximin ordering
  dataa[[l]] <- dataa[[l]][, order]
  locs[[l]] <- gene_data$location[[l]][order, ]
  
  n_locs = ncol(dataa[[l]]) # Number of spatial locations
  num_reps = nrow(dataa[[l]]) # Number of replicates
  
  temp <- rdist(locs[[l]])
  
  nearest_neighbors[[l]] <- find_nn_dist(temp, 50)
}

counts <- list()
for(l in 1:length(gene_data$data)){
  counts[[l]] <- gene_data$data[[l]][common.genes,]  
}


# Source all the required libraries
suppressPackageStartupMessages(library(MixMatrix))
suppressPackageStartupMessages(library(MCMCpack))
suppressPackageStartupMessages(library(NPVecchia))
suppressPackageStartupMessages(library(fields)) # Contains rdist function
suppressPackageStartupMessages(library(adaptMCMC))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mvtnorm))


# Quiet Function 
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

suppressPackageStartupMessages(source("functions_replicates.R"))

n_locs = c(ncol(dataa[[1]]),
           ncol(dataa[[2]])) # Number of spatial locations

num_reps = nrow(dataa[[1]]) # Number of replicates

m = 10 # Number of nearest neighbors 
N = nrow(dataa[[1]]) # Define the number of replicates

nu = 3 * N # Degrees of freedom for the IW prior

Psi = diag(1, num_reps) # Set scaling parameter for IW prior
# Lambda = as.matrix(forceSymmetric(riwish(v = nu, S = Psi))) # Starting value for the IW sampling (taken to be equal to true value)
Lambda = as.matrix(forceSymmetric(riwish(v = nu, S = Psi))) # Starting value for the IW sampling (taken to be equal to true value)

n_iterations = 1500

# Define lists to store the samples
U_samples <- replicate(length(dataa), list(), simplify = FALSE)
d_samples <- replicate(length(dataa), list(), simplify = FALSE)

post <- list()

theta_samples <- matrix(NA, nrow = n_iterations, ncol = 3)

Lambda_samples <- list()

#A scale matrix similar to this worked well for our application
scale_mat <- matrix(c(0.05, -0.04, 0, -0.04, 0.05, 0, 0, 0, 0.01), nc = 3)

#number of samples
nruns <- 2
log_likelihood <- list()
init_theta = c(1,-1,0)

start_time <- Sys.time()
for(i in 1:n_iterations){
  # Printing the iterations
  if(i == 1){
    cat(paste0("Iteration: ", i, "\n"))
  }
  if(i %% floor((10/100)*(n_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", i, "\n"))
  }
  # Find the optimal thetas as a function of Lambda from the previous iteration
  thetas_new <- quiet(adaptMCMC::MCMC(minus_loglikeli_my_new, datumT = dataa, 
                                      NNarrayT = nearest_neighbors,
                                      Lambda = Lambda,
                                      m = m,
                                      init = init_theta, negativ = FALSE,
                                      scale = scale_mat, adapt = TRUE, 
                                      acc.rate = 0.234, n = nruns, 
                                      showProgressBar = FALSE)$samples[nruns, ])
  
  theta_samples[i, ] = thetas_new
  for(r in 1:length(dataa)){
    # Convert the hyper-parameters thetas to priors which is a function of Lambda from previous iteration
    priors = thetas_to_priors_my(thetas_new, n = nrow(nearest_neighbors[[r]]), m = m)
    # get initial posterior sample
    post[[r]] = get_posts_my(datum = dataa[[r]], priors = priors, Lambda = Lambda, NNarray = nearest_neighbors[[r]])
    
    # Sample posterior sparse matrix U as a function of Lambda from previous iteration
    U_samples[[r]][[i]] = samp_posts_my_new(post[[r]], NNarray = nearest_neighbors[[r]], bayesian = TRUE)
    # Store the diagonal elements of U
    d_samples[[r]][[i]] = diag(U_samples[[r]][[i]])
  }
  
  
  SS = matrix(0, nrow = N, ncol = N) # Define the sample sum of square matrix
  # ALTERNATE
  for(r in 1:length(dataa)){
    for(ind in 1:n_locs[r]){
      if(ind == 1){
        mean.normal = rep(0, num_reps)
        d = 1/((d_samples[[r]][[i]][ind])^2)
        SS = SS + tcrossprod(dataa[[r]][, ind])/d
      }else{
        X = as.matrix(post[[r]][[5]][[ind]])
        d = 1/((d_samples[[r]][[i]][ind])^2)
        u =  sqrt(d) * U_samples[[r]][[i]][as.numeric(post[[r]][[6]][[ind]]), ind]  
        mean.normal = as.numeric(X %*% u)
        SS = SS + (tcrossprod(dataa[[r]][, ind] - mean.normal))/d
      }
    }
  }
  
  # Draw sample from IW i.e. from full conditional [Lambda|-]
  Lambda = riwish(v = nu + sum(n_locs), S = as.matrix(forceSymmetric(Psi + SS)))
  # Store Lambda as sample
  Lambda_samples[[i]] = as.matrix(forceSymmetric(Lambda))
  log_like <- replicate(length(dataa), list(), simplify = FALSE)
  
  for(r in 1:length(dataa)){
    for(ind in 1:n_locs[[r]]){
      if(ind == 1){
        mean.normal = rep(0, num_reps)
        d = 1/((d_samples[[r]][[i]][ind])^2)
      }else{
        X = as.matrix(post[[r]][[5]][[ind]])
        d = 1/((d_samples[[r]][[i]][ind])^2)
        u =  sqrt(d) * U_samples[[r]][[i]][as.numeric(post[[r]][[6]][[ind]]), ind]  
        mean.normal = as.numeric(X %*% u)
        
      }
      
      log_like[[r]][[ind]] = dmvnorm(x = dataa[[r]][, ind], mean = mean.normal, sigma = d * Lambda_samples[[i]], log = TRUE)
    }
    log_likelihood[[i]] = sum(unlist(log_like))
  }
}# End of Gibbs Sampling
end_time <- Sys.time()
time.taken = end_time - start_time

U1_post = Matrix(0, nrow = n_locs[1], ncol = n_locs[1], sparse = TRUE)
U2_post = Matrix(0, nrow = n_locs[2], ncol = n_locs[2], sparse = TRUE)

Lambda_post = matrix(0, nrow = num_reps, ncol = num_reps)
theta_post = 0

burn = 500 # Number of burn-in

thin = 1

samples <- seq(from = (burn + 1), to = n_iterations, by = thin)

for(i in samples){
  U1_post = U1_post + U_samples[[1]][[i]]  # Sum of all matrices from Posterior of U1
  U2_post = U2_post + U_samples[[2]][[i]]  # Sum of all matrices from Posterior of U2
  Lambda_post = Lambda_post + Lambda_samples[[i]] # Sum of all matrices from Posterior of Lambda
}

U1_post = U1_post/length(samples) # Posterior Mean as an estimate of U1
U2_post = U2_post/length(samples) # Posterior Mean as an estimate of U2

Sigma1_post = solve(t(U1_post),  solve(U1_post)) # Posterior estimate of Sigma1
Sigma2_post = solve(t(U2_post),  solve(U2_post)) # Posterior estimate of Sigma2

Lambda_post = Lambda_post/length(samples) # Posterior Mean as an estimate of Lambda

log_like <- unlist(log_likelihood[samples])

library(tidyverse)
ll_plot <- data.frame(x = 1:length(log_like), y = log_like) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Traceplot of log-likelihood", x = "Iteration post burn-in", y = "") + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )

library(forecast)
ACF_plot <-  ggAcf(x = log_like, lag.max = 40) + ggtitle("ACF of log-likelihood") + labs(y = "") +
  # ylim(c(-0.35, 0.65)) +
  theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )


################################################################################
# SPECTRAL CLUSTERING AND PLOT FOR LIGHT SAMPLE 1
################################################################################
library(Matrix)
# Change the Spatial Correslation matrices for each of the four samples
S <- cov2cor(as.matrix(Sigma1_post))
# Convert correlation matrix to lie between 0 and 1
S <- (S + 1)/2
# Get the Normalized Graph Laplacian matrix
D = diag(rowSums(S))
L = diag(1, nrow = nrow(S)) - solve(D,S)
# Perform eigen analysis on the Normalized Graph Laplacian matrix
eigL = eigen(L)
n_locs <- ncol(S)

ten_eigen_values <- sort(eigL$values[(n_locs-10):(n_locs-1)])
library(kneedle)

knees <- kneedle(x = 1:10, y = ten_eigen_values, decreasing = FALSE, concave = TRUE)[1]

# Extract the eigen vectors corresponding to the (k-1) smallest eigen values 
V = eigL$vectors[, (n_locs-knees):(n_locs-1) ]
# To find the number of cluster
cent = 2:20
WSS <- 0

# Get the Within Sum of Squares for different number of clusters

for(ind in 1:length(cent)){
  K.means.run <- list()
  SS <- 0
  for(i in 1:100){
    K.means.run[[i]] = kmeans(V, centers = cent[ind], iter.max = 30)
    SS[i] = K.means.run[[i]]$betweenss/ K.means.run[[i]]$totss
  }
  WSS[ind] <- K.means.run[[which.max(SS)]]$tot.withinss
}

center <- kneedle(x = cent, y = WSS, decreasing = TRUE, concave = FALSE)[1]
SS <- 0

K.means.run <- list()
for(i in 1:100){
  K.means.run[[i]] = kmeans(V, centers = center, iter.max = 30)
  SS[i] = K.means.run[[i]]$betweenss/ K.means.run[[i]]$totss
}

K.means.post <- K.means.run[[which.max(SS)]]$cluster

indexes = as.numeric(rownames(locs[[1]]))

cluster_labels <- read.csv("./STARmap/light_replicate_1/light_replicate_1.cov.csv")

cluster_labels_exc <- cluster_labels %>% filter(ClusterName %in% c("eL2/3", "eL5", "eL6", "eL4") ) %>% dplyr::select(CellID, ClusterID, ClusterName)

cluster_labels_exc <- cluster_labels_exc[cluster_labels_exc$CellID %in% intersect(indexes, cluster_labels_exc$CellID), ]
est_labels <- 0

for(i in 1:length(cluster_labels_exc$CellID)){
  est_labels[i] <- K.means.post[which(indexes %in% cluster_labels_exc$CellID[i])]
}

true_cluster <- as.factor(cluster_labels_exc$ClusterName)

# load gene count matrices
gene_data <-NULL


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

Data_new_final = data.frame(Cluster = factor(K.means.post), Data_new)
Data_new_final = st_sf(data.frame(Data_new_final))

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
  geom_sf(data = Data_new_final, aes(fill = Cluster), alpha = 0.8)  +
  theme_minimal() +
  labs(title = "", subtitle = paste0("ARI = ", round(aricode::ARI(true_cluster, est_labels), 4))) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        plot.title = element_text(face = "bold", size=14),
        plot.subtitle = element_text(face = "bold", size=14),
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
# GENE-NETWORK ESTIMATION
################################################################################
Z1 = dataa[[1]] %*% U1_post
Z2 = dataa[[2]] %*% U2_post

Z = as.matrix(cbind(Z1, Z2))

library(huge)
out = huge(t(Z), method = "glasso", nlambda = 100, cov.output = TRUE)
out.select = huge.select(out, criterion = "stars")

plot(out.select)

lambda.best <- out.select$opt.lambda
lambda.best

Lambda_corr_post = -cov2cor(solve(Lambda_post))
diag(Lambda_corr_post) = 1

Precision_mat = out.select$opt.icov

head(Precision_mat, c(6,6))

Precision_mat[!lower.tri(Precision_mat, diag = FALSE)] <- 0
cell.connection = apply(Precision_mat, 1, function(x){which(x != 0, arr.ind = TRUE)})

connections <- list()
for(i in 1:length(cell.connection)){
  if(! identical(cell.connection[[i]], integer(0))){
    connections[[i]] <- cbind(cell.connection[[i]], i)
  }
  
}
connections <- Filter(Negate(is.null), connections)
network = NULL
for(i in 1:length(connections)){
  network = rbind(network, connections[[i]])
}

colnames(network) <- c("gene number1", "gene number2")

# Get the genes that form a network and those who are isolated 
gene.numbers = as.vector(t(network))
gene.network = rownames(dataa[[1]])[as.vector(t(network))]
gene.no.network = rownames(dataa[[1]])[setdiff(x = 1:nrow(dataa[[1]]), y = unique(gene.numbers))]

library(igraph)
# Create the Graph
g <- graph(edges = gene.network, isolates = gene.no.network, directed = FALSE)
set.seed(123465)
par(mfrow = c(1,1))
# Plot the Graph
library(qgraph)
e <- get.edgelist(g,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                       area= 10 *(vcount(g)^2),
                                       repulse.rad= (vcount(g)^3),
                                       max.delta = 5 * (vcount(g)))
plot(g,layout=l,vertex.size = 16, vertex.label.cex = 1.1, vertex.shape = "none", edge.color="gray8")

################################################################################
# SPATIAL GENE EXPRESSION BOX-PLOT (CO-EXPRESSED GENES)
################################################################################
gene.interested = unique(gene.network)

plots_boxplot <- list()

for(clus in 1:length(unique(K.means.post))){
  long_data = gather(as.data.frame(t(dataa[[1]][which(rownames(dataa[[1]]) %in% gene.interested), which(K.means.post == clus)])), "gene", "expr")
  
  plots_boxplot[[clus]] <- long_data %>% ggplot(aes(x = gene, y = expr, fill = gene)) + 
    geom_boxplot(width = 0.3) + theme_minimal() + 
    labs(x = "", y = "Expression", title = paste0("Cluster ", clus)) + 
    theme(text = element_text(size = 14),
          axis.text.x = element_text(vjust = .5, face = "bold"),
          # axis.text.x = element_blank(),
          axis.text.y = element_text(vjust = .5, face = "bold"),
          plot.title = element_text(face = "bold", size=14),
          plot.subtitle = element_text(face = "bold", size=12),
          panel.grid.major = element_blank(),
          legend.title = element_text(face = "bold", size=18),
          legend.text = element_text(size = 16),
          legend.position = "none"
    ) 
} 

library(gridExtra)
n <- length(plots_boxplot)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plots_boxplot, ncol=nCol))

################################################################################
# SPATIAL GENE EXPRESSION BOX-PLOT (TOP 10 SPATIALLY VARYING GENES)
################################################################################
gene.interested = c("Bgn",     "Cux2",    "Egr1",    "eRNA3",   "Fos",     "Mgp",
                   "Nov",     "Pcp4",    "Plcxd2",   "Tnfaip6")

plots_boxplot <- list()

for(clus in 1:length(unique(K.means.post))){
  long_data = gather(as.data.frame(t(dataa[[1]][which(rownames(dataa[[1]]) %in% gene.interested), which(K.means.post == clus)])), "gene", "expr")
  
  plots_boxplot[[clus]] <- long_data %>% ggplot(aes(x = gene, y = expr, fill = gene)) + 
    geom_boxplot(width = 0.3) + theme_minimal() + 
    labs(x = "", y = "Expression", title = paste0("Cluster ", clus)) + 
    theme(text = element_text(size = 14),
          axis.text.x = element_text(vjust = .5, face = "bold"),
          # axis.text.x = element_blank(),
          axis.text.y = element_text(vjust = .5, face = "bold"),
          plot.title = element_text(face = "bold", size=14),
          plot.subtitle = element_text(face = "bold", size=12),
          panel.grid.major = element_blank(),
          legend.title = element_text(face = "bold", size=18),
          legend.text = element_text(size = 16),
          legend.position = "none"
    ) 
} 

library(gridExtra)
n <- length(plots_boxplot)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plots_boxplot, ncol=nCol))



################################################################################
# SPECTRAL CLUSTERING AND PLOT FOR LIGHT SAMPLE 2
################################################################################
# Change the Spatial Correslation matrices for each of the four samples
S <- cov2cor(as.matrix(Sigma2_post))
# Convert correlation matrix to lie between 0 and 1
S <- (S + 1)/2
# Get the Normalized Graph Laplacian matrix
D = diag(rowSums(S))
L = diag(1, nrow = nrow(S)) - solve(D,S)
# Perform eigen analysis on the Normalized Graph Laplacian matrix
eigL = eigen(L)
n_locs <- ncol(S)

ten_eigen_values <- sort(eigL$values[(n_locs-10):(n_locs-1)])
library(kneedle)

knees <- kneedle(x = 1:10, y = ten_eigen_values, decreasing = FALSE, concave = TRUE)[1]

# Extract the eigen vectors corresponding to the (k-1) smallest eigen values 
V = eigL$vectors[, (n_locs-knees):(n_locs-1) ]
# To find the number of cluster
cent = 2:20
WSS <- 0

# Get the Within Sum of Squares for different number of clusters

for(ind in 1:length(cent)){
  K.means.run <- list()
  SS <- 0
  for(i in 1:100){
    K.means.run[[i]] = kmeans(V, centers = cent[ind], iter.max = 30)
    SS[i] = K.means.run[[i]]$betweenss/ K.means.run[[i]]$totss
  }
  WSS[ind] <- K.means.run[[which.max(SS)]]$tot.withinss
}

center <- kneedle(x = cent, y = WSS, decreasing = TRUE, concave = FALSE)[1]
SS <- 0

K.means.run <- list()
for(i in 1:100){
  K.means.run[[i]] = kmeans(V, centers = center, iter.max = 30)
  SS[i] = K.means.run[[i]]$betweenss/ K.means.run[[i]]$totss
}

K.means.post <- K.means.run[[which.max(SS)]]$cluster

indexes = as.numeric(rownames(locs[[2]]))

cluster_labels <- read.csv("./STARmap/light_replicate_2/light_replicate_2.cov.csv")

cluster_labels_exc <- cluster_labels %>% filter(ClusterName %in% c("eL2/3", "eL5", "eL6", "eL4") ) %>% dplyr::select(CellID, ClusterID, ClusterName)

cluster_labels_exc <- cluster_labels_exc[cluster_labels_exc$CellID %in% intersect(indexes, cluster_labels_exc$CellID), ]
est_labels <- 0

for(i in 1:length(cluster_labels_exc$CellID)){
  est_labels[i] <- K.means.post[which(indexes %in% cluster_labels_exc$CellID[i])]
}

true_cluster <- as.factor(cluster_labels_exc$ClusterName)

# load gene count matrices
gene_data <-NULL


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

Data_new_final = data.frame(Cluster = factor(K.means.post), Data_new)
Data_new_final = st_sf(data.frame(Data_new_final))

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
  geom_sf(data = Data_new_final, aes(fill = Cluster), alpha = 0.8)  +
  theme_minimal() +
  labs(title = "", subtitle = paste0("ARI = ", round(aricode::ARI(true_cluster, est_labels), 4))) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        plot.title = element_text(face = "bold", size=14),
        plot.subtitle = element_text(face = "bold", size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(face = "bold", size=18),
        legend.text = element_text(size = 16)
  ) + 
  guides(colour = guide_legend(title="clusters",
                               title.hjust = 0.5)) +
  scale_fill_manual(values = myvalues)

plot_Light2



################################################################################
# GENE-CORRELATION HEATMAP
################################################################################
library(corrplot)
Cor_mat = cov2cor(Lambda_post)

rownames(Cor_mat) <- rownames(dataa[[1]])
colnames(Cor_mat) <- rownames(dataa[[1]])
corrplot(Cor_mat, method = 'color', order = 'hclust', type = "lower", diag = FALSE, tl.cex = 1.0, is.corr = FALSE)

