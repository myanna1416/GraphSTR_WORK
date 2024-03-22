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
#install.packages("dict")
library(dict)#dict library in python
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

# Example:
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


