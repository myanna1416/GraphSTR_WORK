#utils

mclust_R <- function(adata, num_cluster, modelNames = 'EEE', used_misc = 'emb_pca', random_seed = 2020) {
  # Clustering using the mclust algorithm
  # The parameters are the same as those in the R package mclust
  
  set.seed(random_seed)
  # Convert numpy array to R matrix
  adata_matrix <- as.matrix(adata@misc[[used_misc]])
  
  # Call Mclust function
  res <- Mclust(data = adata_matrix, G = num_cluster, modelNames = modelNames)
  mclust_res <- res$classification
  
  # Update adata object
  adata$obs$mclust <- as.integer(mclust_res)
  adata$obs$mclust <- as.factor(adata$obs$mclust)
  
  return(adata)
}

clustering <- function(adata, n_clusters = 7, radius = 50, key = 'emb', method = 'mclust', 
                       start = 0.1, end = 3.0, increment = 0.01, refinement = FALSE) {
  # Perform PCA
  pca <- prcomp(adata@misc[[key]], center = TRUE, scale. = TRUE, retx = TRUE)
  embedding <- pca$x
  adata@misc[['emb_pca']] <- embedding
  
  if (method == 'mclust') {
    # Assuming mclust_R() is a placeholder for an actual mclust implementation in R
    adata <- mclust_R(adata, used_misc = 'emb_pca', num_cluster = n_clusters)
    adata@misc[['domain']] <- as.factor(adata@misc[['mclust']])
  } else if (method == 'leiden') {
    # Convert data to igraph object
    dist_matrix <- dist(embedding)
    graph <- graph_from_adjacency_matrix(as_adjacency_matrix(dist_matrix, sparse = TRUE), mode = "undirected", weighted = TRUE)
    # Find resolution parameter
    res <- search_res(adata, n_clusters, use_rep = 'emb_pca', method = method, start = start, end = end, increment = increment)
    # Run Leiden algorithm
    community <- leiden(graph, resolution_parameter = res)
    adata$obs[['domain']] <- factor(community$membership)
  } else if (method == 'louvain') {
    # Convert data to igraph object
    dist_matrix <- dist(embedding)
    graph <- graph_from_adjacency_matrix(as_adjacency_matrix(dist_matrix, sparse = TRUE), mode = "undirected", weighted = TRUE)
    # Find resolution parameter
    res <- search_res(adata, n_clusters, use_rep = 'emb_pca', method = method, start = start, end = end, increment = increment)
    # Run Louvain algorithm
    community <- cluster_louvain(graph, resolution = res)
    adata$obs[['domain']] <- factor(community$membership)
  }
  
  if (refinement) {
    # Implement refinement logic if required
    new_type <- refine_label(adata, radius, key = 'domain')
    adata$obs[['domain']] <- new_type
  }
  
  return(adata)
}

# Helper functions (need to be defined or adapted)
search_res <- function(adata, n_clusters, use_rep, method, start, end, increment) {
  # Dummy placeholder for finding an optimal resolution parameter
  return(1)  # Adjust based on actual analysis
}
#removes the noise from data set
refine_label <- function(adata, radius = 50, key = 'label') {
  n_neigh <- radius
  new_type <- character(0) # c()empty array in R 
  old_type <- as.character(adata$obs[[key]])
  
  # Calculate distance
  position <- adata@misc[['spatial']]
  distance <- ot::dist(position, position, metric = 'euclidean')
  n_cell <- nrow(distance)
  
  for (i in 1:n_cell) {
    vec <- distance[i, ]
    index <- order(vec)
    neigh_type <- character(0)
    for (j in 2:(n_neigh + 1)) {
      neigh_type <- c(neigh_type, old_type[index[j]])
    }
    max_type <- names(sort(table(neigh_type), decreasing = TRUE)[1])
    new_type <- c(new_type, max_type)
  }
  
  new_type <- as.character(new_type)
  # adata$obs[['label_refined']] <- new_type
  
  return(new_type)
}

extract_top_value <- function(map_matrix, retain_percent = 0.1) {
  # Filter out cells with low mapping probability
  
  # Calculate the number of cells to retain
  top_k <- retain_percent * ncol(map_matrix)
  
  # Retain top values for each spot
  output <- map_matrix * (order(order(map_matrix)) >= ncol(map_matrix) - top_k)
  
  return(output)
}

construct_cell_type_matrix <- function(adata_sc) {
  label <- 'cell_type'
  n_type <- length(unique(adata_sc$obs[[label]]))
  zeros <- matrix(0, nrow = nrow(adata_sc), ncol = n_type)
  cell_type <- unique(adata_sc$obs[[label]])
  cell_type <- as.character(cell_type)
  cell_type <- sort(cell_type)
  mat <- as.data.frame(zeros, row.names = adata_sc$obs_names, col.names = cell_type)
  for (cell in adata_sc$obs_names) {
    ctype <- adata_sc$obs[cell, label]
    mat[cell, as.character(ctype)] <- 1
  }
  # res <- colSums(mat)
  return(mat)
}


project_cell_to_spot <- function(adata, adata_sc, retain_percent = 0.1) {
  # Project cell types onto ST data using mapped matrix in adata.misc
  
  # Read map matrix
  map_matrix <- adata@misc[['map_matrix']]  # spot x cell
  
  # Extract top-k values for each spot
  map_matrix <- extract_top_value(map_matrix)  # filtering by spot
  
  # Construct cell type matrix
  matrix_cell_type <- construct_cell_type_matrix(adata_sc)
  matrix_cell_type <- as.matrix(matrix_cell_type)
  
  # Projection by spot-level
  matrix_projection <- map_matrix %*% matrix_cell_type
  
  # Rename cell types
  cell_type <- unique(adata_sc$obs[['cell_type']])
  cell_type <- as.character(cell_type)
  cell_type <- sort(cell_type)
  # cell_type <- gsub(' ', '_', cell_type)
  df_projection <- as.data.frame(matrix_projection, row.names = adata$obs_names, col.names = cell_type)  # spot x cell type
  
  # Normalize by row (spot)
  df_projection <- t(apply(df_projection, 1, function(x) x / sum(x)))
  
  # Add projection results to adata
  adata$obs[, cell_type] <- df_projection
  
  return(NULL)  # No return value (None in Python)
}

search_res <- function(adata, n_clusters, method = 'leiden', use_rep = 'emb', start = 0.1, end = 3.0, increment = 0.01) {
  # Searching corresponding resolution according to given cluster number
  
  cat('Searching resolution...\n')
  label <- 0
  sc.pp.neighbors(adata, n_neighbors = 50, use_rep = use_rep)
  for (res in sort(seq(start, end, by = increment), decreasing = TRUE)) {
    if (method == 'leiden') {
      sc.tl.leiden(adata, random_state = 0, resolution = res)
      count_unique <- length(unique(as.character(adata$obs[['leiden']])))
      cat('resolution =', res, ', cluster number =', count_unique, '\n')
    } else if (method == 'louvain') {
      sc.tl.louvain(adata, random_state = 0, resolution = res)
      count_unique <- length(unique(as.character(adata$obs[['louvain']])))
      cat('resolution =', res, ', cluster number =', count_unique, '\n')
    }
    if (count_unique == n_clusters) {
      label <- 1
      break
    }
  }
  
  stopifnot(label == 1, 'Resolution is not found. Please try a bigger range or smaller step!')
  
  return(res)
}
