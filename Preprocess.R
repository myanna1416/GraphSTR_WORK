#preprocess 
library(transport)
library(ROptimalTransport)
library(torch)
library(Seurat)# Load data and create a Seurat object, analogous to using scanpy
library(Matrix) 
library(FNN) 
filter_with_overlap_gene <- function(adata, adata_sc) {
  # Assuming adata and adata_sc are Seurat objects
  
  # Check for the presence of 'highly_variable' in both datasets
  if (!"highly_variable" %in% colnames(adata@meta.data)) {
    stop("'highly_variable' are not existed in adata!")
  } else {
    adata <- adata[, adata@meta.data$highly_variable, drop = FALSE]
  }
  
  if (!"highly_variable" %in% colnames(adata_sc@meta.data)) {
    stop("'highly_variable' are not existed in adata_sc!")
  } else {
    adata_sc <- adata_sc[, adata_sc@meta.data$highly_variable, drop = FALSE]
  }
  
  # Find overlap in genes, assuming gene names are row names
  genes <- intersect(rownames(adata), rownames(adata_sc))
  cat('Number of overlap genes:', length(genes), "\n")
  
  # Store the overlapping genes in metadata
  metadata(adata)$overlap_genes <- genes
  metadata(adata_sc)$overlap_genes <- genes
  
  # Subset both datasets to overlapping genes
  adata <- adata[, genes, drop = FALSE]
  adata_sc <- adata_sc[, genes, drop = FALSE]
  
  return(list(adata, adata_sc))
}
permutation <- function(feature) {
  set.seed(41)
  
  # Generate indices for the feature matrix
  ids <- seq_len(nrow(feature))
  
  # Randomly permute the indices
  ids <- sample(ids)
  
  # Permute the feature matrix rows according to the shuffled indices
  feature_permutated <- feature[ids, , drop = FALSE]
  
  return(feature_permutated)
}
construct_interaction <- function(adata, n_neighbors = 3) {
  # Assuming 'spatial' coordinates are stored similarly to Python's AnnData `obsm`
  position <- adata$obsm$spatial
  
  # Calculate distance matrix using Euclidean distance
  distance_matrix <- as.matrix(dist(position))
  n_spot <- nrow(distance_matrix)
  
  adata$obsm$distance_matrix <- distance_matrix
  
  # Initialize the interaction matrix with zeros
  interaction <- matrix(0, nrow = n_spot, ncol = n_spot)
  
  # Find k-nearest neighbors
  for (i in seq_len(n_spot)) {
    vec <- distance_matrix[i, ]
    distance <- order(vec)  # Get indices that would sort the vector
    for (t in seq_len(n_neighbors)) {
      y <- distance[t + 1]  # Skip the first element as it is the point itself
      interaction[i, y] <- 1
    }
  }
  
  adata$obsm$graph_neigh <- interaction
  
  # Transform adjacency matrix to be symmetrical
  adj <- interaction + t(interaction)  # Transpose and add
  adj <- ifelse(adj > 1, 1, adj)  # Cap values at 1
  
  adata$obsm$adj <- adj
  
  return(adata)
}

construct_interaction_KNN <- function(adata, n_neighbors = 3) {
  # Assuming 'spatial' is a matrix or data frame in adata$obsm
  position <- adata$obsm$spatial
  n_spot <- nrow(position)
  
  # Fit nearest neighbors model and find neighbors
  nbrs <- get.knnx(data = position, query = position, k = n_neighbors + 1)
  
  # Prepare indices for interaction matrix construction
  indices <- nbrs$nn.index
  
  # R uses 1-based indexing, adjusting for that here
  x <- rep(seq_len(n_spot), each = n_neighbors)
  y <- as.vector(t(indices[, -1]))  # Remove the first column and flatten the rest
  
  # Initialize the interaction matrix with zeros
  interaction <- matrix(0, nrow = n_spot, ncol = n_spot)
  
  # Populate the interaction matrix
  for (i in seq_along(x)) {
    interaction[x[i], y[i]] <- 1
  }
  
  adata$obsm$graph_neigh <- interaction
  
  # Transform adjacency matrix to be symmetrical
  adj <- interaction + t(interaction)  # Add transpose of itself
  adj <- ifelse(adj > 1, 1, adj)  # Ensure all values are capped at 1
  
  adata$obsm$adj <- adj
  
  cat('Graph constructed!\n')
  return(adata)
}
library(Seurat)

preprocess <- function(adata) {
  # Find highly variable genes similar to Seurat's approach which aligns closely with "seurat_v3"
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 3000)
  
  # Normalize data to a total expression of 10,000 per cell
  adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Log-transform the data
  adata <- ScaleData(adata, features = rownames(adata), scale.max = 10, verbose = FALSE)
  
  # Print statement to indicate completion of preprocessing
  cat("Preprocessing completed!\n")
  
  return(adata)
}
get_feature <- function(adata, deconvolution = FALSE) {
  
  # Conditionally select data based on 'highly_variable' markers or not
  if (deconvolution) {
    adata_Vars <- adata
  } else {
    # Assume adata is a list-like structure or a Seurat object with var metadata
    highly_variable_indices <- which(adata$var$highly_variable)
    adata_Vars <- adata[, highly_variable_indices]
  }
  
  # Check if the data matrix is sparse and convert to dense if necessary
  if (inherits(adata_Vars$X, "sparseMatrix")) {
    feat <- as.matrix(adata_Vars$X)  # Convert sparse matrix to dense matrix
  } else {
    feat <- adata_Vars$X
  }
  
  # Data augmentation by permuting features
  feat_a <- permutation(feat)  # Ensure permutation function is defined in R
  
  # Store features and augmented features back into adata object
  adata$obsm$feat <- feat
  adata$obsm$feat_a <- feat_a
  
  return(adata)
}
add_contrastive_label <- function(adata) {
  # Number of observations (spots)
  n_spot <- nrow(adata$X)  # Assuming adata is structured such that adata$X contains the data matrix
  
  # Create matrices of ones and zeros
  one_matrix <- matrix(1, nrow = n_spot, ncol = 1)
  zero_matrix <- matrix(0, nrow = n_spot, ncol = 1)
  
  # Concatenate matrices to form the label matrix
  label_CSL <- cbind(one_matrix, zero_matrix)
  
  # Store the contrastive label in the adata object under 'obsm'
  if (!"obsm" %in% names(adata)) {
    adata$obsm <- list()
  }
  adata$obsm$label_CSL <- label_CSL
  
  return(adata)
}

normalize_adj <- function(adj) {
  # Convert the adjacency matrix to a sparse matrix in coordinate format
  adj <- as(adj, "CsparseMatrix")  # Ensuring it is a sparse matrix
  
  # Calculate the sum of each row (degree of nodes)
  rowsum <- rowSums(adj)
  
  # Compute D^(-1/2) for each element
  d_inv_sqrt <- 1 / sqrt(rowsum)
  d_inv_sqrt[is.infinite(d_inv_sqrt)] <- 0  # Replace Inf with 0
  
  # Create a diagonal matrix of D^(-1/2)
  d_mat_inv_sqrt <- Diagonal(x = d_inv_sqrt)
  
  # Perform symmetric normalization: A' = D^(-1/2) * A * D^(-1/2)
  adj_normalized <- d_mat_inv_sqrt %*% adj %*% d_mat_inv_sqrt
  
  # Convert the normalized adjacency matrix back to a regular matrix
  return(as.matrix(adj_normalized))
}
library(Matrix)

preprocess_adj <- function(adj) {
  # Normalize the adjacency matrix using the previously defined normalize_adj function
  adj_normalized <- normalize_adj(adj)
  
  # Add an identity matrix to include self-connections
  # Using diagonal() to create an identity matrix of appropriate size
  adj_normalized <- adj_normalized + Diagonal(x = rep(1, nrow(adj_normalized)))
  
  # Return the normalized adjacency matrix with self-loops
  return(adj_normalized)
}

sparse_mx_to_torch_sparse_tensor <- function(sparse_mx) {
  # Convert the input to a sparse matrix in COO format, if not already
  if (!inherits(sparse_mx, "dgCMatrix")) {
    sparse_mx <- as(sparse_mx, "dgCMatrix")
  }
  
  # Get the indices and values from the sparse matrix
  i <- which(sparse_mx@x > 0, arr.ind = TRUE)  # Find indices with non-zero entries
  values <- sparse_mx@x  # Extract non-zero values
  
  # Convert indices to 0-based for Torch
  indices <- torch_tensor(as.integer(i) - 1, dtype = torch_int64())
  
  # Create the values tensor
  values <- torch_tensor(values, dtype = torch_float32())
  
  # Define the shape of the tensor
  shape <- torch_size(as.integer(dim(sparse_mx)))
  
  # Create the sparse tensor
  sparse_tensor <- torch_sparse_tensor(indices, values, shape)
  
  return(sparse_tensor)
}

preprocess_adj_sparse <- function(adj) {
  # Ensure the matrix is in COO format suitable for operations
  if (!inherits(adj, "dgCMatrix")) {
    adj <- as(adj, "dgCMatrix")
  }
  
  # Add self-loops by adding an identity matrix
  adj_ <- adj + Diagonal(x = rep(1, nrow(adj)))
  
  # Compute the sum of each row (degree of nodes)
  rowsum <- as.array(rowSums(adj_))
  
  # Compute D^(-1/2) for the degree matrix
  degree_mat_inv_sqrt <- Diagonal(x = 1 / sqrt(rowsum))
  degree_mat_inv_sqrt@x[is.infinite(degree_mat_inv_sqrt@x)] <- 0
  
  # Normalize the adjacency matrix symmetrically
  adj_normalized <- degree_mat_inv_sqrt %*% adj_ %*% t(degree_mat_inv_sqrt)
  
  # Convert the normalized adjacency matrix to a sparse tensor
  return(sparse_mx_to_torch_sparse_tensor(adj_normalized))
}
fix_seed <- function(seed) {
  # Set the seed for R's base random number generator, which affects sample() and runif() among others
  set.seed(seed)
  
  # For reproducibility in 'torch' package
  torch::torch_manual_seed(seed)
  
  # Assuming use of GPU via torch; these settings would be set if R torch supports them
  # For now, R's 'torch' package does not directly support setting CUDA seeds and determinism as explicitly as Python
  
  # Setting environment variables related to reproducibility
  Sys.setenv(PYTHONHASHSEED = as.character(seed))
  
  # No direct equivalent in R for CUDA-specific settings, but if necessary, could use system-specific ways to set configurations
  # Sys.setenv(CUBLAS_WORKSPACE_CONFIG = ':4096:8') # Uncomment if such environment variable settings are needed in R
  
  # Note: R does not use cuDNN in the same manner, so these options are illustrative rather than functional
  if (requireNamespace("torch", quietly = TRUE)) {
    torch::torch_set_deterministic(TRUE)
  }
}

