#preprocess 


#' Filter datasets by overlapping highly variable genes
#'
#' This function takes two Seurat object datasets and filters both by the overlapping
#' "highly variable" genes present in both. It stops with an error if the "highly variable"
#' attribute does not exist in either dataset.
#'
#' @param adata A Seurat object containing single-cell data, expected to have a "highly_variable" attribute in `meta.data`.
#' @param adata_sc A second Seurat object containing single-cell data to be compared with `adata`, also expected to have a "highly_variable" attribute in `meta.data`.
#' @return A list of two Seurat objects, `adata` and `adata_sc`, filtered to contain only genes that are marked as highly variable in both datasets.
#' @examples
#' # Assuming 's1' and 's2' are Seurat objects with 'highly_variable' metadata:
#' results <- filter_with_overlap_gene(s1, s2)
#' print(results[[1]])  # Print filtered 's1'
#' print(results[[2]])  # Print filtered 's2'
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


#' Permute the rows of a matrix
#'
#' This function shuffles the rows of a given matrix randomly. It is useful for creating
#' permutation samples of dataset rows. The function sets a seed for reproducibility.
#'
#' @param feature A matrix whose rows are to be permuted.
#' @return A matrix with the same dimensions as `feature`, but with rows randomly permuted.
#' @examples
#' data_matrix <- matrix(1:12, ncol=3, byrow=TRUE)
#' print("Original matrix:")
#' print(data_matrix)
#' print("Permuted matrix:")
#' permuted_matrix <- permutation(data_matrix)
#' print(permuted_matrix)
permute_features <- function(feature) {
  set.seed(41) #later check on this seed number 
  
  # Generate indices for the feature matrix
  ids <- seq_len(nrow(feature))
  
  # Randomly permute the indices
  ids <- sample(ids)
  
  # Permute the feature matrix rows according to the shuffled indices
  feature_permutated <- feature[ids, , drop = FALSE]
  
  return(feature_permutated)
}


#' Construct an interaction matrix based on spatial proximity
#'
#' This function constructs an interaction matrix for given spatial data by calculating 
#' the Euclidean distance between points and identifying the k-nearest neighbors for each point.
#' The resulting interaction matrix is symmetric, with ones indicating neighbor relationships
#' and zeros otherwise. This matrix can be used for graph-based analyses.
#'
#' @param adata An object containing spatial coordinates in `adata@misc$spatial`.
#' @param n_neighbors The number of nearest neighbors to identify for each point; default is 3.
#' @return The input `adata` object, modified to include the new `distance_matrix`, `graph_neigh`, 
#'         and `adj` (symmetrical adjacency matrix) within `adata@misc`.
#' @examples
#' adata <- list(misc = list(spatial = matrix(runif(20), ncol = 2)))
#' adata <- construct_interaction(adata, n_neighbors = 3)
#' print(adata@misc$adj)
construct_interaction <- function(adata, n_neighbors = 3) {
  position <- GetTissueCoordinates(adata,scale = NULL)
  position <- position[,c(2,1)]
  adata@misc$spatial <- position
  
  # Calculate distance matrix using Euclidean distance
  distance_matrix <- as.matrix(dist(position))
  n_spot <- nrow(distance_matrix)
  
  adata@misc$distance_matrix <- distance_matrix
  
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
  
  adata@misc$graph_neigh <- interaction
  
  # Transform adjacency matrix to be symmetrical
  adj <- interaction + t(interaction)  # Transpose and add
  adj <- ifelse(adj > 1, 1, adj)  # Cap values at 1
  
  adata@misc$adj <- adj
  
  return(adata)
}


#' Preprocess single-cell RNA-seq data
#'
#' This function preprocesses single-cell RNA-seq data by performing several steps:
#' 1. Identifying highly variable genes using the variance stabilizing transformation (VST) method.
#' 2. Normalizing data to a total expression of 10,000 per cell using log-normalization.
#' 3. Log-transforming the data and scaling it.
#' These steps prepare the data for further analysis such as clustering or principal component analysis.
#'
#' @param adata A Seurat object containing single-cell RNA-seq data.
#' @return A Seurat object with preprocessing steps applied, including updated fields for variable
#'         features and normalized, scaled data.
#' @examples
#' # Load example data
#' data("pbmc_small")
#' # Preprocess data
#' pbmc_small <- preprocess(pbmc_small)
#' # View processed data
#' print(pbmc_small)
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




#' Extract and process features from single-cell RNA-seq data
#'
#' This function extracts features from an `adata` object, focusing on highly variable genes unless
#' deconvolution is required. It optionally converts sparse matrix representations to dense, and
#' performs data augmentation by permuting the features.
#'
#' @param adata A data object, typically a Seurat object, containing single-cell RNA-seq data with
#'              predefined 'highly_variable' markers in the `var` metadata.
#' @param deconvolution A logical indicating whether to use all data or only highly variable genes
#'                      for feature extraction. Default is `FALSE`, meaning only highly variable genes are used.
#' @return The `adata` object with additional slots for features (`feat`) and permuted features (`feat_a`)
#'         in the `misc` component.
#' @examples
#' # Assuming 'sce' is a Seurat object with 'highly_variable' metadata:
#' sce <- get_feature(sce)
#' print(sce@misc$feat)  # Features
#' print(sce@misc$feat_a)  # Augmented features
get_feature <- function(adata, deconvolution = FALSE) {
  
  # Conditionally select data based on 'highly_variable' markers or not
  if (deconvolution) {
    adata_Vars <- adata
  } else {
    # Assume adata is a list-like structure or a Seurat object with var metadata
    highly_variable_features <- VariableFeatures(adata)
    #adata_Vars <- adata[highly_variable_features,]
    adata_Vars <- subset(adata, features = highly_variable_features)
  }
  
  assay_data <- GetAssayData(adata_Vars, slot = "counts")
  
  # Check if the data matrix is sparse and convert to dense if necessary
  if (inherits(assay_data, "dgCMatrix")) {  # 'dgCMatrix' is the class for sparse matrices in R
    feat <- as.matrix(assay_data)  # Convert sparse matrix to dense matrix
  } else {
    feat <- assay_data  # Directly use the matrix if it is already dense
  }
  
  # Data augmentation by permuting features
  feat_a <- permute_features(feat)  
  
  # Store features and augmented features back into adata object
  adata@misc$feat <- feat
  adata@misc$feat_a <- feat_a
  #adata[['feat']] <- CreateAssayObject(feat)  # Create new assay for original features
  #adata[['feat_a']] <- CreateAssayObject(feat_a) # Create new assay for augmented features
  
  return(adata)
}

#' Add contrastive labels to a data object
#'
#' This function creates and appends a contrastive label matrix to the 'misc' component of 
#' the provided data object. Each label matrix consists of two columns, one filled with ones
#' and the other with zeros, representing binary class labels for contrastive learning tasks.
#'
#' @param adata A data object, typically structured similarly to a Seurat or anndata object, 
#'              expected to contain a data matrix in `adata$X`.
#' @return The `adata` object with a new `label_CSL` matrix appended to the `misc` component.
#'         This matrix has one column of ones and one column of zeros.
#' @examples
#' # Assume 'adata' is a data object with a matrix `adata$X` containing observational data
#' adata <- list(X = matrix(rnorm(20), ncol = 2, nrow = 10))
#' adata <- add_contrastive_label(adata)
#' print(adata@misc$label_CSL)  # Print the contrastive labels
add_contrastive_label <- function(adata) {
  
  #replaceing Adata$X with GetAssayData(adata)
  # Number of observations (spots)
  n_spot <- nrow(GetAssayData(adata))  # Assuming adata is structured such that adata$X contains the data matrix
  
  # Create matrices of ones and zeros
  one_matrix <- matrix(1, nrow = n_spot, ncol = 1)
  zero_matrix <- matrix(0, nrow = n_spot, ncol = 1)
  
  # Concatenate matrices to form the label matrix
  label_CSL <- cbind(one_matrix, zero_matrix)
  
  # Store the contrastive label in the adata object under 'misc'
  if (!"misc" %in% names(adata)) {
    adata@misc <- list()
  }
  adata@misc$label_CSL <- label_CSL
  
  return(adata)
}

#' Normalize an adjacency matrix
#'
#' This function normalizes an adjacency matrix using the symmetric normalization method.
#' The function converts the matrix to a sparse format, calculates the degree of each node,
#' computes the inverse square root of the degree matrix, and applies the normalization
#' symmetrically to the adjacency matrix.
#'
#' @param adj A numeric matrix representing the adjacency matrix of a graph.
#' @return A matrix representing the normalized adjacency matrix, which is symmetrically
#'         normalized using the formula D^(-1/2) * A * D^(-1/2), where D is the diagonal degree matrix.
#' @examples
#' # Create an example adjacency matrix
#' adj_matrix <- matrix(c(0, 1, 0, 0, 
#'                        1, 0, 1, 0,
#'                        0, 1, 0, 1,
#'                        0, 0, 1, 0), nrow = 4, byrow = TRUE)
#' # Normalize the adjacency matrix
#' normalized_adj <- normalize_adj(adj_matrix)
#' print(normalized_adj)
normalize_adj <- function(adj) {
  # Convert the adjacency matrix to a sparse matrix in coordinate format
  adj <- Matrix::as(adj, "CsparseMatrix")  # Ensuring it is a sparse matrix
  
  # Calculate the sum of each row (degree of nodes)
  rowsum <- rowSums(adj)
  
  # Compute D^(-1/2) for each element
  d_inv_sqrt <- 1 / sqrt(rowsum)
  d_inv_sqrt[is.infinite(d_inv_sqrt)] <- 0  # Replace Inf with 0
  
  # Create a diagonal matrix of D^(-1/2)
  d_mat_inv_sqrt <- Matrix::Diagonal(x = d_inv_sqrt, sparse = TRUE)
  
  # Perform symmetric normalization: A' = D^(-1/2) * A * D^(-1/2)
  adj_normalized <- d_mat_inv_sqrt %*% adj %*% d_mat_inv_sqrt
  
  # Convert the normalized adjacency matrix back to a regular matrix
  return(as.matrix(adj_normalized))
}

#' Preprocess an adjacency matrix by adding self-connections
#'
#' This function takes a given adjacency matrix, normalizes it using the `normalize_adj` function,
#' and adds self-connections by including an identity matrix. This preprocessing is common in
#' graph neural networks to enhance node features by including self-loops.
#'
#' @param adj A numeric matrix representing the adjacency matrix of a graph.
#' @return A matrix representing the preprocessed adjacency matrix, which includes self-loops
#'         and is normalized using the formula D^(-1/2) * A * D^(-1/2) with added identity matrix.
#' @examples
#' # Create an example adjacency matrix
#' adj_matrix <- matrix(c(0, 1, 0, 0,
#'                        1, 0, 1, 0,
#'                        0, 1, 0, 1,
#'                        0, 0, 1, 0), nrow = 4, byrow = TRUE)
#' # Preprocess the adjacency matrix
#' preprocessed_adj <- preprocess_adj(adj_matrix)
#' print(preprocessed_adj)
preprocess_adj <- function(adj) {
  # Normalize the adjacency matrix using the previously defined normalize_adj function
  adj_normalized <- normalize_adj(adj)
  
  # Add an identity matrix to include self-connections
  # Using diagonal() to create an identity matrix of appropriate size
  adj_normalized <- adj_normalized + Diagonal(x = rep(1, nrow(adj_normalized)))
  
  # Return the normalized adjacency matrix with self-loops
  return(adj_normalized)
}

#' Convert a sparse matrix to a Torch sparse tensor
#'
#' This function converts a given sparse matrix into a COO format sparse tensor compatible
#' with the Torch library in R. The function ensures that the matrix is in a suitable sparse
#' format (dgCMatrix), extracts the non-zero indices and values, and constructs a sparse tensor.
#' Indices are adjusted to zero-based indexing required by Torch.
#'
#' @param sparse_mx A sparse matrix, ideally of class dgCMatrix. If not, it will be
#'                  converted to this format.
#' @return A Torch sparse tensor constructed from the non-zero elements of the input sparse matrix.
#' @examples
#' # Create a sparse matrix
#' library(Matrix)
#' sparse_matrix <- sparseMatrix(i = c(1, 3, 4), j = c(1, 2, 2), x = c(4, 5, 2))
#'
#' # Convert to Torch sparse tensor
#' sparse_tensor <- sparse_mx_to_torch_sparse_tensor(sparse_matrix)
#' print(sparse_tensor)
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

#' Preprocess a sparse adjacency matrix for graph neural networks
#'
#' This function preprocesses a sparse adjacency matrix by first ensuring it is in the 
#' compressed sparse column (CSC) format, then adding self-loops by including an identity matrix,
#' and normalizing it using the symmetric normalization method. The output is suitable for use in 
#' graph-based machine learning models, specifically formatted as a sparse tensor.
#'
#' @param adj A numeric matrix or a sparse matrix in any format, representing the adjacency matrix of a graph.
#' @return A sparse tensor, formatted for use with Torch, representing the normalized adjacency matrix 
#'         with self-loops added. This tensor is suitable for graph neural network models.
#' @examples
#' # Create an example adjacency matrix
#' library(Matrix)
#' adj_matrix <- Matrix(c(0, 1, 0, 0, 
#'                        1, 0, 1, 0,
#'                        0, 1, 0, 1,
#'                        0, 0, 1, 0), nrow = 4, byrow = TRUE, sparse = TRUE)
#'
#' # Preprocess the adjacency matrix
#' processed_adj_tensor <- preprocess_adj_sparse(adj_matrix)
#' print(processed_adj_tensor)
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

#' Set the seed for reproducibility across R and Torch environments
#'
#' This function sets the seed for R's base random number generator and the Torch package,
#' aiming to ensure reproducibility in operations that involve random number generation.
#' It also sets an environment variable for Python hash seeding, although its direct
#' impact in R is limited unless specifically interfaced with Python code.
#' Note: This function does not configure GPU-specific settings as Python might allow,
#' but prepares the environment for deterministic operations where possible.
#'
#' @param seed An integer value to be used as the seed for random number generators.
#' @examples
#' fix_seed(123)
#' # Run some random operations
#' print(runif(5))  # Random numbers from uniform distribution
#' print(torch::torch_randn(c(2, 3)))  # Random tensor
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

