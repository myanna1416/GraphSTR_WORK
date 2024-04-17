#READ ME 
#This is the work space for GRAPHSTR 
#Use these imports first when starting up the program 
#This is hwo you would import the functions from other scripts into your main script 


##DUMMY CODE 

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

GraphST$set("public", "train", function() {
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



#DUMMY CODE Zor.R 

      # Check for 'highly_variable' in adata$var
      if (!("highly_variable" %in% names(self$adata$var))) {
        self$adata <- preprocess(self$adata)
      }
      
      # Check for 'adj' in adata$obsm and construct interaction matrices
      if (!("adj" %in% names(self$adata$obsm))) {
        if (self$datatype %in% c("Stereo", "Slide")) {
          self$adata <- construct_interaction_KNN(self$adata)
        } else {
          self$adata <- construct_interaction(self$adata)
        }
      }
      
      # Check for 'label_CSL' in adata$obsm and add contrastive labels
      if (!("label_CSL" %in% names(self$adata$obsm))) {
        self$adata <- add_contrastive_label(self$adata)
      }
      
      # Check for 'feat' in adata$obsm and extract features
      if (!("feat" %in% names(self$adata$obsm))) {
        self$adata <- get_feature(self$adata)
      }
