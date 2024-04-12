#For data manipulation (similar to pandas)
library(dplyr)
#For PyTorch operations, similar to PyTorch in Python
library(torch)
#R Equivalent of Python's tqdm
library(pbapply)
pblapply(1:100, function(x) Sys.sleep(0.1))
# Creating an array in R
my_array <- array(data = 1:24, dim = c(2, 3, 4))
#R Equivalent of Python's time 
Sys.time()
runif(1, 0, 1)
# For sparse matrix operations
library(Matrix) 
#install.packages("R6")  # Install R6 package 
library(R6)  # Load the R6 library
#install.packages("Dict")
library(Dict)#dict library in python
#install.packages("torch")
library(torch) #PyTorch in R

GraphST <- R6::R6Class(
  "GraphST",
  
  public = list(
    adata = NULL,
    adata_sc = NULL,  # Adding adata_sc for scRNA-seq data
    device = "cpu",
    learning_rate = 0.001,
    learning_rate_sc = 0.01,
    weight_decay = 0.00,
    epochs = 600,
    dim_input = 3000,  
    dim_output = 64,  
    random_seed = 41,
    alpha = 10,
    beta = 1,
    theta = 0.1,
    lamda1 = 10,
    lamda2 = 1,
    deconvolution = FALSE,
    datatype = "10X",
    
    initialize = function(adata, adata_sc = NULL, device = "cpu", learning_rate = 0.001, learning_rate_sc = 0.01,
                          weight_decay = 0.00, epochs = 600, random_seed = 41, alpha = 10, 
                          beta = 1, theta = 0.1, lamda1 = 10, lamda2 = 1, 
                          deconvolution = FALSE, datatype = "10X") {
      self$adata <- adata  
      self$adata_sc <- adata_sc  # Initialize adata_sc
      self$device <- device
      self$learning_rate <- learning_rate
      self$learning_rate_sc <- learning_rate_sc
      self$weight_decay <- weight_decay
      self$epochs <- epochs
      self$random_seed <- random_seed
      self$alpha <- alpha
      self$beta <- beta
      self$theta <- theta
      self$lamda1 <- lamda1
      self$lamda2 <- lamda2
      self$deconvolution <- deconvolution
      self$datatype <- datatype
      
      set.seed(self$random_seed)  # Fix the seed for reproducibility
    }
  )
)

my_adata <- list()  # Assume this is your prepared spatial data
my_adata_sc <- list()  # Assume this is your scRNA-seq data

graph_st_instance <- GraphST$new(adata = my_adata, adata_sc = my_adata_sc)

# Check for 'highly_variable' in adata$var
if (!"highly_variable" %in% keys(self$adata$var)) {
  self$adata <- preprocess(self$adata)
}

# Check for 'adj' in adata$obsm and construct interaction matrices
if (!"adj" %in% keys(self$adata$obsm)) {
  if (self$datatype %in% c("Stereo", "Slide")) {
    self$adata <- construct_interaction_KNN(self$adata)
  } else {
    self$adata <- construct_interaction(self$adata)
  }
}

# Check for 'label_CSL' in adata$obsm and add contrastive labels
if (!"label_CSL" %in% keys(self$adata$obsm)) {
  self$adata <- add_contrastive_label(self$adata)
}

# Check for 'feat' in adata$obsm and extract features
if (!"feat" %in% keys(self$adata$obsm)) {
  self$adata <- get_feature(self$adata)
}

# Convert 'feat' to a torch tensor and move to the specified device
self$features <-
  torch_tensor(as.array(self$adata$obsm$feat), 
               dtype = torch_float32)$to(device = self$device)
# Convert 'feat_a' to a torch tensor and move to the specified device
self$features_a <-
  torch_tensor(as.array(self$adata$obsm$feat_a), 
               dtype = torch_float32)$to(device = self$device)
# Convert 'label_CSL' to a torch tensor and move to the specified device
self$label_CSL <-
  torch_tensor(as.array(self$adata$obsm$label_CSL), 
               dtype = torch_float32)$to(device = self$device)
# For 'adj', we simply reference it as it does not necessarily need conversion for this context
self$adj <- self$adata$obsm$adj
n <- nrow(self$adj)
self$graph_neigh <-
  torch_tensor(as.array(self$adata$obsm$graph_neigh) + diag(rep(1, n)), 
               dtype = torch_float32)$to(device = self$device)
self$dim_input <- dim(self$features)[2]
self$dim_output <- dim_output

self$adj <- preprocess_adj(self$adj)
self$adj <-
  torch_tensor(self$adj, dtype = torch_float32)$to(device = self$device)
#-- you could use this alt. --
#--self$adj <- torch_tensor(as.array(self$adj), dtype = torch_float32)$to(device = self$device)
if (self$deconvolution) {
  self$adata_sc <-
    self$adata_sc$clone(deep = TRUE)  # Assuming adata_sc can be a R6 object
}
# Assuming self$feat_sc and self$feat_sp are numeric matrices or data frames in R

# Replace NA (equivalent to NaN in Python) with 0 for self$feat_sc
self$feat_sc[is.na(self$feat_sc)] <- 0

# Convert to a torch tensor and move to the specified device
self$feat_sc <- torch_tensor(self$feat_sc, dtype = torch_float32)$to(device = self$device)

# Repeat the process for self$feat_sp
self$feat_sp[is.na(self$feat_sp)] <- 0
self$feat_sp <- torch_tensor(self$feat_sp, dtype = torch_float32)$to(device = self$device)

# Check if self$adata_sc exists
if (!is.null(self$adata_sc)) {
  # Assuming self$feat_sc is a torch tensor, use dim() to get dimensions
  self$dim_input <- dim(self$feat_sc)[2]
}

# Set the number of cells and spots based on observations in adata_sc and adata
# Assuming adata_sc and adata are lists or similar objects with n_obs (number of observations) property
self$n_cell <- self$adata_sc$n_obs
self$n_spot <- self$adata$n_obs

GraphST$set("public", "train", function() {
  if (self$datatype %in% c("Stereo", "Slide")) {
    # Encoder_sparse is imported
    self$model <- Encoder_sparse(self$dim_input, self$dim_output, self$graph_neigh)$to(device = self$device)
  } else {
    #Encoder is imported
    self$model <- Encoder(self$dim_input, self$dim_output, self$graph_neigh)$to(device = self$device)
  }
  
  # Using nn_bce_with_logits_loss for the loss function
  self$loss_CSL <- nn_bce_with_logits_loss()
  
  # Setting up the optimizer, assuming self$model has a parameters method similar to PyTorch
  self$optimizer <- optim_adam(self$model$parameters(), lr = self$learning_rate, weight_decay = self$weight_decay)
  
  cat("Begin to train ST data...\n")
  self$model$train()
  
})

GraphSTR$set("public", "train", function() {
  cat("Begin to train ST data...\n")
  
  pblapply(1:self$epochs, function(epoch) {
    self$model$train()
    
    self$features_a <- permutation(self$features)  # Ensure permutation is defined
    
    # Assuming model returns a list with required elements
    model_output <- self$model(self$features, self$features_a, self$adj)
    self$hiden_feat <- model_output[[1]]
    self$emb <- model_output[[2]]
    ret <- model_output[[3]]
    ret_a <- model_output[[4]]
    
    self$loss_sl_1 <- self$loss_CSL(ret, self$label_CSL)
    self$loss_sl_2 <- self$loss_CSL(ret_a, self$label_CSL)
    self$loss_feat <- nnf_mse_loss(self$features, self$emb)
    
    loss <-  self$alpha * self$loss_feat + self$beta * (self$loss_sl_1 + self$loss_sl_2)
    
    self$optimizer$zero_grad()
    loss$backward()
    self$optimizer$step()
  }, .progress = TRUE)  #requires pbapply
  
  cat("Optimization finished for ST data!\n")
})

GraphSTR$set("public", "evaluate", function() {
  torch::no_grad({
    self$model$eval()
    
    if (self$deconvolution) {
      model_output <- self$model(self$features, self$features_a, self$adj)
      self$emb_rec <- model_output[[2]]  # Assuming the second element is the desired output
      
      return(self$emb_rec)
    } else {
      model_output <- self$model(self$features, self$features_a, self$adj)
      self$emb_rec <- model_output[[2]]  # Similarly assuming the second element
      
      if (self$datatype %in% c("Stereo", "Slide")) {
        self$emb_rec <- nnf_normalize(self$emb_rec, p = 2, dim = 2)$detach()$cpu()$numpy()
      } else {
        self$emb_rec <- self$emb_rec$detach()$cpu()$numpy()
      }
      
      # Assuming self$adata is a list or environment and has an element $obsm where 'emb' can be stored
      self$adata$obsm$emb <- self$emb_rec
      
      return(self$adata)
    }
  })
}) 
GraphSTR$set("public", "train_sc", function() {
  # Assuming Encoder_sc is already defined in R
  self$model_sc <- Encoder_sc(self$dim_input, self$dim_output)$to(device = self$device)
  
  # Setting up the optimizer for the single-cell model
  self$optimizer_sc <- optim_adam(self$model_sc$parameters(), lr = self$learning_rate_sc)
  
  cat('Begin to train scRNA data...\n')
  
  # Assuming the epochs are defined; replace tqdm with pbapply for progress bar in R
  pbapply::pblapply(1:self$epochs, function(epoch) {
    self$model_sc$train()
    
    emb <- self$model_sc(self$feat_sc)
    loss <- nnf_mse_loss(emb, self$feat_sc)
    
    self$optimizer_sc$zero_grad()
    loss$backward()
    self$optimizer_sc$step()
  }, .progress = TRUE)  # Optional: Progress bar
  
  cat("Optimization finished for cell representation learning!\n")
})
GraphSTR$set("public", "evaluate_sc", function() {
  torch::no_grad({
    self$model_sc$eval()
    emb_sc <- self$model_sc(self$feat_sc)
    
    return(emb_sc)
  })
})
GraphST$set("public", "train_map", function() {
  emb_sp <- self$train()
  emb_sc <- self$train_sc()
  
  # Assuming adata and adata_sc are list-like objects with an 'obsm' list for embeddings
  self$adata$obsm$emb_sp <- emb_sp$detach()$cpu()$numpy()
  self$adata_sc$obsm$emb_sc <- emb_sc$detach()$cpu()$numpy()
  
  # Normalize features for consistency between ST and scRNA-seq
  emb_sp <- nnf_normalize(emb_sp, p = 2, eps = 1e-12, dim = 2)
  emb_sc <- nnf_normalize(emb_sc, p = 2, eps = 1e-12, dim = 2)
  
  # Initialize the model for mapping, assuming Encoder_map is defined
  self$model_map <- Encoder_map(self$n_cell, self$n_spot)$to(device = self$device)
  
  # Set up the optimizer
  self$optimizer_map <- optim_adam(self$model_map$parameters(), lr = self$learning_rate, weight_decay = self$weight_decay)
  
  cat('Begin to learn mapping matrix...\n')
})

GraphST$set("public", "train_mapping_matrix", function(emb_sp, emb_sc) {
  cat("Begin to learn mapping matrix...\n")
  
  # Training loop with progress bar
  pbapply::pblapply(1:self$epochs, function(epoch) {
    self$model_map$train()
    
    # Invoke the mapping model to get the current mapping matrix
    self$map_matrix <- self$model_map()
    
    # Calculate reconstruction loss and Noise-Contrastive Estimation (NCE) loss
    losses <- self$calculate_loss(emb_sp, emb_sc) # Assuming this function is implemented to return a list of losses
    loss_recon <- losses$loss_recon
    loss_NCE <- losses$loss_NCE
    
    # Combine losses according to specified weights
    total_loss <- self$lamda1 * loss_recon + self$lamda2 * loss_NCE
    
    # Optimization steps
    self$optimizer_map$zero_grad()
    total_loss$backward()
    self$optimizer_map$step()
  }, .progress = TRUE)  # Enables the progress bar for visual feedback
  
  cat("Mapping matrix learning finished!\n")
})
GraphST$set("public", "finalize_mapping", function(emb_sp, emb_sc) {
  torch::no_grad({
    self$model_map$eval()
    
    # In R's torch, tensors are already on the CPU and don't need conversion to numpy arrays
    # Applying softmax to the mapping matrix; assuming self$map_matrix is a torch tensor
    map_matrix <- nnf_softmax(self$map_matrix, dim = 2) # Note: dim=2 in R's 1-based indexing matches dim=1 in Python's 0-based
    
    # Assigning processed embeddings and mapping matrix back to adata structures
    # Assuming self$adata and self$adata_sc are list-like structures that can hold these values
    self$adata$obsm$emb_sp <- emb_sp
    self$adata_sc$obsm$emb_sc <- emb_sc
    self$adata$obsm$map_matrix <- t(map_matrix) # Transpose for spot x cell organization
    
    return(list(self$adata, self$adata_sc))
  })
})
GraphST$set("public", "calculate_loss", function(emb_sp, emb_sc) {
  # Perform softmax normalization by cell on the mapping matrix
  map_probs <- nnf_softmax(self$map_matrix, dim = 2) # Adjusted for R's 1-based indexing
  
  # Calculate predicted spatial spots by matrix multiplication
  self$pred_sp <- torch_matmul(map_probs$t(), emb_sc)
  
  # Reconstruction loss: mean squared error between predicted and actual spatial spots
  loss_recon <- nnf_mse_loss(self$pred_sp, emb_sp, reduction = 'mean')
  
  # Noise-Contrastive Estimation (NCE) loss
  loss_NCE <- self$Noise_Cross_Entropy(self$pred_sp, emb_sp)
  
  return(list(loss_recon = loss_recon, loss_NCE = loss_NCE))
})
GraphST$set("public", "Noise_Cross_Entropy", function(pred_sp, emb_sp) {
  # Calculate cosine similarity between predicted and actual embeddings
  mat <- self$cosine_similarity(pred_sp, emb_sp)
  
  # Sum of exponentiated similarities for normalization, excluding self-similarities
  k <- torch_exp(mat)$sum(dim = 2) - torch_exp(torch_diag(mat))
  
  # Calculate the numerator for positive pairs, using the graph neighborhood matrix
  p <- torch_exp(mat)
  p <- torch_mul(p, self$graph_neigh)$sum(dim = 2)
  
  # Average probability for positive pairs over all possibilities
  ave <- torch_div(p, k)
  loss <- -torch_log(ave)$mean()
  
  return(loss)
})
GraphST$set("public", "cosine_similarity", function(pred_sp, emb_sp) {
  # Matrix multiplication between pred_sp and the transpose of emb_sp
  M <- torch_matmul(pred_sp, emb_sp$t())
  
  # Calculating the norm for pred_sp and emb_sp
  Norm_c <- torch_norm(pred_sp, p = 2, dim = 2)
  Norm_s <- torch_norm(emb_sp, p = 2, dim = 2)
  
  # Reshaping norms for broadcasting and performing outer product
  Norm <- torch_matmul(Norm_c$view(dim(c(pred_sp$size()[1], 1))),
                       Norm_s$view(dim(c(emb_sp$size()[1], 1)))$t()) + -5e-12
  
  # Element-wise division to calculate cosine similarity
  M <- torch_div(M, Norm)
  
  # Handling NaN values by replacing them with a specific value
  if (torch_any(torch_isnan(M))) {
    M <- torch_where(torch_isnan(M), torch_full_like(M, 0.4868), M)
  }
  
  return(M)
})




