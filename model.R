#Model 

Discriminator <- nn_module(
  classname = "Discriminator",
  
  initialize = function(n_h) {
    # Initialize the bilinear layer
    self$f_k <- nn_bilinear(in1_features = n_h, in2_features = n_h, out_features = 1)
    
    # Apply custom weights initialization
    self$apply(self$weights_init)
  },
  
  weights_init = function(m) {
    # Xavier uniform initialization for bilinear layers
    if (inherits(m, "nn_bilinear")) {
      nn_init_xavier_uniform_(m$weight)
    }
  },
  
  forward = function(c, h_pl, h_mi, s_bias1 = NULL, s_bias2 = NULL) {
    # Expand c to match the dimensions of h_pl
    c_x <- c$expand_as(h_pl)
    
    # Compute scores using the bilinear layer
    sc_1 <- self$f_k(h_pl, c_x)
    sc_2 <- self$f_k(h_mi, c_x)
    
    # Add optional biases to the scores
    if (!is.null(s_bias1)) {
      sc_1 <- sc_1 + s_bias1
    }
    if (!is.null(s_bias2)) {
      sc_2 <- sc_2 + s_bias2
    }
    
    # Concatenate scores along dimension 2 (1-based indexing in R)
    logits <- torch_cat(list(sc_1, sc_2), dim = 2)
    
    return(logits)
  }
)

AvgReadout <- nn_module(
  classname = "AvgReadout",
  
  initialize = function() {
    # Nothing needed in the initializer for this particular case
  },
  
  forward = function(emb, mask = NULL) {
    if (is.null(mask)) {
      mask <- torch_ones(size = list(nrow(emb), 1))
    }
    
    # Perform matrix multiplication to sum embeddings based on the mask
    vsum <- torch_mm(mask, emb)
    
    # Sum across rows to get total counts per row
    row_sum <- torch_sum(mask, dim = 2)
    
    # Expand row_sum to match dimensions of vsum for element-wise division
    row_sum <- row_sum$expand(c(vsum$size()[2], row_sum$size()[1]))$t()
    
    # Element-wise division to get the average
    global_emb <- vsum / row_sum

    # Normalize the global embeddings
    normalized_emb <- nnf_normalize(global_emb, p = 2, dim = 2)
    
    return(normalized_emb)
  }
)

Encoder <- nn_module(
  classname = "Encoder",
  
  initialize = function(in_features, out_features, graph_neigh, dropout = 0.0, act = nnf_relu) {
    self$in_features <- in_features
    self$out_features <- out_features
    self$graph_neigh <- graph_neigh
    self$dropout <- dropout
    self$act <- act
    
    # Define weights as parameters
    self$weight1 <- nn_parameter(torch_randn(in_features, out_features))
    self$weight2 <- nn_parameter(torch_randn(out_features, in_features))
    
    # Reset parameters using Xavier uniform initialization
    self$reset_parameters()
    
    # Additional modules -- is thie correct?
    # self$disc <- Discriminator$new(out_features)
    self$disc <- Discriminator(out_features)
    
    self$sigm <- nn_sigmoid()
    # self$read <- AvgReadout$new()
    self$read <- AvgReadout()
  },
  
  reset_parameters = function() {
    nn_init_xavier_uniform_(self$weight1)
    nn_init_xavier_uniform_(self$weight2)
  },
  
  forward = function(feat, feat_a, adj) {
    #Encoder_forward(feat, feat_a, adj)
    z <- nnf_dropout(feat, p = self$dropout, training = self$training)
    z <- torch_mm(z, self$weight1)
    z <- torch_mm(adj, z)
    
    hiden_emb <- z

    h <- torch_mm(z, self$weight2)
    h <- torch_mm(adj, h)
    
    emb <- self$act(z)
    
    z_a <- nnf_dropout(feat_a, p = self$dropout, training = self$training)
    z_a <- torch_mm(z_a, self$weight1)
    z_a <- torch_mm(adj, z_a)
    emb_a <- self$act(z_a)
    
    g <- self$read(emb, self$graph_neigh)
    g <- self$sigm(g)
    
    g_a <- self$read(emb_a, self$graph_neigh)
    g_a <- self$sigm(g_a)

    ret <- self$disc(g, emb, emb_a)
    ret_a <- self$disc(g_a, emb_a, emb)
    
    return(list(hiden_emb = hiden_emb, h = h, ret = ret, ret_a = ret_a))
  }
)

Encoder_sparse <- nn_module(
  classname = "Encoder_sparse",
  
  initialize = function(in_features, out_features, graph_neigh, dropout = 0.0, act = nnf_relu) {
    self$in_features <- in_features
    self$out_features <- out_features
    self$graph_neigh <- graph_neigh
    self$dropout <- dropout
    self$act <- act
    
    # Define weights as parameters
    self$weight1 <- nn_parameter(torch_randn(in_features, out_features))
    self$weight2 <- nn_parameter(torch_randn(out_features, in_features))

    # Reset parameters using Xavier uniform initialization
    self$reset_parameters()
    
    # Additional modules
    self$disc <- Discriminator$new(out_features)
    self$sigm <- nn_sigmoid()
    self$read <- AvgReadout$new()
  },
  
  reset_parameters = function() {
    nn_init_xavier_uniform_(self$weight1)
    nn_init_xavier_uniform_(self$weight2)
  },
  
  forward = function(feat, feat_a, adj) {
    cat("Initial feat dimensions: ", dim(feat), "\n")
    z <- nnf_dropout(feat, p = self$dropout, training = self$training)
    z <- torch_mm(z, self$weight1)
    cat("Dimensions after weight1 multiplication: ", dim(z), "\n")
    z <- torch_spmm(adj, z)  # Using sparse matrix multiplication
    cat("Dimensions after adj multiplication: ", dim(z), "\n")
    hiden_emb <- z
    
    h <- torch_mm(z, self$weight2)
    h <- torch_spmm(adj, h)  # Using sparse matrix multiplication again
    
    emb <- self$act(z)
    
    z_a <- nnf_dropout(feat_a, p = self$dropout, training = self$training)
    z_a <- torch_mm(z_a, self$weight1)
    z_a <- torch_spmm(adj, z_a)  # Sparse matrix multiplication for augmented features
    emb_a <- self$act(z_a)
    
    g <- self$read(emb, self$graph_neigh)
    g <- self$sigm(g)
    
    g_a <- self$read(emb_a, self$graph_neigh)
    g_a <- self$sigm(g_a)
    
    ret <- self$disc(g, emb, emb_a)
    ret_a <- self$disc(g_a, emb_a, emb)
    
    return(list(hiden_emb = hiden_emb, h = h, ret = ret, ret_a = ret_a))
  }
)
Encoder_sc <- nn_module(
  classname = "Encoder_sc",
  
  initialize = function(dim_input, dim_output, dropout = 0.0, act = nnf_relu) {
    self$dim_input <- dim_input
    self$dim_output <- dim_output
    self$dropout <- dropout
    self$act <- act
    
    # Setting dimensions for intermediate layers
    self$dim1 <- 256
    self$dim2 <- 64
    self$dim3 <- 32
    
    # Initialize weights for encoder
    self$weight1_en <- nn_parameter(torch_randn(dim_input, self$dim1))
    self$weight2_en <- nn_parameter(torch_randn(self$dim1, self$dim2))
    self$weight3_en <- nn_parameter(torch_randn(self$dim2, self$dim3))
    
    # Initialize weights for decoder
    self$weight1_de <- nn_parameter(torch_randn(self$dim3, self$dim2))
    self$weight2_de <- nn_parameter(torch_randn(self$dim2, self$dim1))
    self$weight3_de <- nn_parameter(torch_randn(self$dim1, dim_input))
    
    # Apply Xavier uniform initialization to all weights
    self$reset_parameters()
  },
  
  reset_parameters = function() {
    nn_init_xavier_uniform_(self$weight1_en)
    nn_init_xavier_uniform_(self$weight2_en)
    nn_init_xavier_uniform_(self$weight3_en)
    
    nn_init_xavier_uniform_(self$weight1_de)
    nn_init_xavier_uniform_(self$weight2_de)
    nn_init_xavier_uniform_(self$weight3_de)
  },
  
  forward = function(x) {
    x <- nnf_dropout(x, p = self$dropout, training = self$training)
    
    # Encoder step
    x <- torch_mm(x, self$weight1_en)
    x <- torch_mm(x, self$weight2_en)
    x <- torch_mm(x, self$weight3_en)
    
    # Decoder step
    x <- torch_mm(x, self$weight1_de)
    x <- torch_mm(x, self$weight2_de)
    x <- torch_mm(x, self$weight3_de)
    
    return(x)
  }
)
Encoder_map <- nn_module(
  classname = "Encoder_map",
  
  initialize = function(n_cell, n_spot) {
    self$n_cell <- n_cell
    self$n_spot <- n_spot
    
    # Initialize the mapping matrix as a parameter
    self$M <- nn_parameter(torch_randn(n_cell, n_spot))
    
    # Apply Xavier uniform initialization to the matrix
    self$reset_parameters()
  },
  
  reset_parameters = function() {
    nn_init_xavier_uniform_(self$M)
  },
  
  forward = function() {
    x <- self$M
    return(x)
  }
)

