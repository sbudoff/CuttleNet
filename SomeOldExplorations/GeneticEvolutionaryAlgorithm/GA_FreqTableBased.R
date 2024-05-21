library(tidyverse)
library(matrixStats)
library(Matrix)

nDELTArow2col2 <- function(expMat){
  # This loss function computes a value between -Inf and Inf where the lower is the better separation
  expMat <- na.omit(expMat)
  n = nrow(expMat)
  classSum2 = sum(rowSums(t(expMat))**2)
  geneSum2 = sum(rowSums(expMat)**2)
  l = n*(geneSum2-classSum2)
}
barcoder <- function(class_mat, verbose = F, binary = F){
  # This Function creates a 'barcode' for each vector of genes defining a class. 
  # This barcode is simply a square matrix of booleans that checks each class against every other class to see if they share a vector of genes.
  # If verbose is true, a visualization of the barcode will be displayed
  class_mat <- na.omit(class_mat)
  n_classes = nrow(t(class_mat))
  n_genes = nrow(class_mat)
  barcode_compare <- data.frame(matrix(ncol = n_classes, nrow = n_classes))
  rownames(barcode_compare) <- colnames(class_mat)
  colnames(barcode_compare) <- colnames(class_mat)
  for (i in 1:n_classes) {
    barcode_i = round(class_mat[,i],1)
    for (j in 1:n_classes) {
      barcode_j = round(class_mat[,j],1)
      if (binary) {
        barcode_compare[i,j] = sum(barcode_i == barcode_j) == nrow(class_mat)
      } else {
        barcode_compare[i,j] = sqrt(sum(((barcode_i-barcode_j)**2)))
      }
    }
  }
  if (verbose) {
    barcode <- barcode_compare %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "class") %>% 
      gather(Class, val, -class) %>% 
      ggplot(aes(Class, class)) + 
      geom_tile(aes(fill = val)) + 
      coord_fixed() + 
      guides(fill = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(barcode)
  }
  sum(rowSums(barcode_compare))
}

RedundantRows <- function(df, precision = 1) {
  df <- round(df,1)
  out <- c(1:nrow(df))
  for (i in 1:nrow(df)) {
    out[i] = length(unique(as.numeric(df[i,])))
  }
  which(out < 1)
}

PlotBubbles <- function(cluster_expression, cluster_expressing, thresh) {
  cluster_by_gene <- cluster_expression %>% 
    rownames_to_column('GENE') %>%
    filter(GENE %in% rownames(thresh)) %>%
    as.data.frame() %>%
    gather(Class, mean, -GENE)
  cluster_by_gene <- cluster_expressing %>% 
    rownames_to_column('GENE') %>%
    filter(GENE %in% rownames(thresh)) %>%
    as.data.frame() %>%
    gather(Class, expressing, -GENE) %>%
    merge(cluster_by_gene, by = c("GENE", "Class"))
  
  # visualize the resulting unique clusters
  visual_cluster <- cluster_by_gene %>% 
    ggplot(aes(Class, GENE)) + 
    geom_point(aes(color = mean, size = expressing)) + 
    coord_fixed() + 
    guides(fill = "none") +
    scale_color_gradient(low = "blue", high = "red") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(visual_cluster)
}

Darwin <- function(data_path, cellClass, key_genes = c(),
                   filtered = F, exponentiated = T, expressingOnly = T,
                   expression_threshold = 0.05, min_expression_var = 0.2,
                   ubiquitous_threshold = 0.5, min_expressing_var = 0.2,
                   min_genes = 9, add_back = 10,
                   n_generations = 100, n_offspring = 100, drop_factor = 10,
                   verbose = T, lossFn = T) {
  # Load the cell class by gene matrices
  load(data_path)
  cluster_expressing <- as.matrix(gene_mats[[cellClass]])
  cluster_expression <- as.matrix(gene_expression_mats[[cellClass]])
  cluster_expressing[is.na(cluster_expressing)] <- 0
  cluster_expression[is.na(cluster_expression)] <- 0
  
  if (filtered) {
    # Remove genes that have no information content
    DropIndex <- which(rowMaxs(as.matrix(cluster_expressing)) < ubiquitous_threshold)
    DropIndex <- unique(c(DropIndex, which(rowMaxs(as.matrix(cluster_expression)) < expression_threshold)))
    DropIndex <- unique(c(DropIndex, which(rowSds(as.matrix(cluster_expressing)) < min_expressing_var)))
    DropIndex <- unique(c(DropIndex, which(rowSds(as.matrix(cluster_expression)) < min_expression_var)))
    # DropIndex <- unique(c(DropIndex, RedundantRows(cluster_expressing[,-1], precision = 1)))
    
    cluster_expression <- cluster_expression[-DropIndex,]
    cluster_expressing <- cluster_expressing[-DropIndex,]
  }
  
  
  # Transform matrices to create fitness matrix
  
  if (exponentiated){
    cluster_expression <- exp(cluster_expression)
    cluster_expressing <- exp(cluster_expressing)
  } else {
    cluster_expression <- cluster_expression
    cluster_expressing <- cluster_expressing
  }
  if (expressingOnly) 
  {
    thresh <- cluster_expressing
  } else {
    thresh <- t(scale(t(cluster_expression)) * scale(t(cluster_expressing)))
  }
  
  # Get number of classes we are dealing with
  n_classes = nrow(t(thresh))
  
  # Get number of genes we are dealing with
  n_genes = nrow(thresh)
  
  # Move the key genes to the head of each matrix so they are never filtered
  key_genes_ind <- which(rownames(thresh) %in% key_genes)
  if (length(key_genes_ind > 0)) {
    n_start <- length(key_genes_ind) + 1
    x <- 1:n_genes
    x <- c(key_genes_ind, x[!x %in% key_genes_ind])
    thresh <- thresh[x,]
  } else {
    n_start <- 1
  }
  # Let's get ready to RUMBAAAAAAL
  thresh_0 = thresh
  gene_list_0 = rownames(thresh)
  generations = 1:n_generations
  genes_dropped = c()
  genes_remaining = c()
  distance_score = 0
  drop_best = NA
  last_best = drop_best
  stuck = 0
  for (gen in generations) {
    n_genes = nrow(thresh) # Recalculate number of genes in this new generation
    n_genes_min = n_genes  # Reset the minimum number of genes an offspring in the new generation has found
    kids = 1:gen*n_offspring
    if (is.na(sum(last_best == drop_best)) )
    {
      stuck = stuck+1
      drop_best = NA
    } else {
      last_best = drop_best
    }
    if (n_genes > min_genes) {
      for (i in kids) {
        n_drop = sample.int(n_genes/drop_factor, 1) # Pick a random number of genes to drop
        drop = sample(n_start:n_genes, n_drop) # Create a random list of genes to drop
        thresh_i = thresh[-drop,] # Drop the genes
        sqCsti = abs(colSums(thresh_i**2))
        sqRsti = abs(rowSums(thresh_i**2))
        if (lossFn) {
          n_genes_i = nrow(thresh_i) # Compute this offspring's number of genes
          distance_score_i = nDELTArow2col2(thresh_i)
          if (distance_score_i < distance_score) { # Check the distance component of fitness, that is this the offspring with the greatest distance between vectors
            # Store fitness values of the best offspring
            n_genes_min = n_genes_i
            drop_best = drop
            distance_score = distance_score_i
            # Run the barcoder checksum on the fittest offspring and see if it gets to make the next generation
            if (length(drop_best) > 0) {
              thresh_i = thresh[-drop_best,]
              print(sprintf("Dropped %s genes in generation %s, %s genes remain. Distance = %s", length(drop_best), gen, nrow(thresh_i), distance_score))
              # checkSum = barcoder(cluster_expressing[-drop_best,-1], verbose = verbose)
              checkSum = barcoder(thresh_i, verbose = verbose)
              if (is.na(checkSum)) {checkSum = n_classes+1}
              if (checkSum > n_classes) {
                thresh = thresh_i
                stuck = 0
              }
              else {
                print("Undifferentiable gene sets identified, generation discarded.")
              }
            } else {
              drop_best = NA
              print(sprintf("No genes dropped in generation %s, %s genes remain", gen, nrow(thresh)))
              stuck = stuck+1
            }
          }
        }
        # else {
        #   if (sum(sqCsti > 0.01) == n_classes) { # Check that each offspring class is defined by at least one gene
        #     n_genes_i = nrow(thresh_i) # Compute this offspring's number of genes
        #     if (n_genes_i < n_genes_min) { # See if this offspring is the most fit, that is has the smallest set of genes
        #       # distance_score_i = sum(as.matrix(dist(thresh_i)))
        #       # distance_score_i <- sum(rowMaxs(thresh_i) /rowMeans(thresh_i))
        #       checkSum = barcoder(thresh_i, verbose = F)
        #       distance_score_i = checkSum
        #       if (distance_score_i > distance_score) { # Check the distance component of fitness, that is this the offspring with the greatest distance between vectors
        #         # Store fitness values of the best offspring
        #         n_genes_min = n_genes_i
        #         drop_best = drop
        #         distance_score = distance_score_i
        #         # Run the barcoder checksum on the fittest offspring and see if it gets to make the next generation
        #         if (length(drop_best) > 0) {
        #           thresh_i = thresh[-drop_best,]
        #           print(sprintf("Dropped %s genes in generation %s, %s genes remain. Distance = %s", length(drop_best), gen, nrow(thresh_i), distance_score))
        #           # checkSum = barcoder(cluster_expressing[-drop_best,-1], verbose = verbose)
        #           checkSum = barcoder(thresh_i, verbose = verbose)
        #           if (checkSum > n_classes) {
        #             thresh = thresh_i
        #           }
        #           else {
        #             print("Undifferentiable gene sets identified, generation discarded.")
        #           }
        #         } else {
        #           drop_best = NA
        #           print(sprintf("No genes dropped in generation %s, %s genes remain", gen, nrow(thresh)))
        #         }
        #       }
        #     }
        #   }
        # }
        
        if (stuck > 10)
        {
          if (nrow(thresh) < min_genes *1.5) 
          {
            genes_dropped = append(genes_dropped, length(drop_best))
            genes_remaining = append(genes_remaining, nrow(thresh_i))
            generations = c(1:length(genes_dropped))
            break
          }
          else {
            # Add random genes back in
            lost_Genes <- setdiff(gene_list_0, rownames(thresh))
            some_Genes <- sample(lost_Genes, add_back)
            gene_list <- c(rownames(thresh), some_Genes)
            thresh <- thresh_0[gene_list,]
            stuck = 0
            print('Random Genes Added back into the pool')
          }
        }
      }
      
      distance_score = 0
    } else {
      genes_dropped = append(genes_dropped, length(drop_best))
      genes_remaining = append(genes_remaining, nrow(thresh_i))
      generations = c(1:length(genes_dropped))
      break
    }
    genes_dropped = append(genes_dropped, length(drop_best))
    genes_remaining = append(genes_remaining, nrow(thresh_i))
  }
  # Visualize the removal of genes over the generations
  gen_plot <- tibble(generations, genes_dropped, genes_remaining) %>%
    ggplot(aes(x = generations)) +
    geom_line(aes(y = genes_dropped), color = 'red') +
    geom_line(aes(y = genes_remaining), color = 'green')
  print(gen_plot)
  # Reorganize the input data based on the thresholding such that all matrices are combined
  if (exponentiated){
    cluster_expression <- log(cluster_expression)
    cluster_expressing <- log(cluster_expressing)
  } 
  
  cluster_expression <- data.frame(cluster_expression)
  cluster_expressing <- data.frame(cluster_expressing)
  
  PlotBubbles(cluster_expression, cluster_expressing, thresh)
  # Print the minimal set of genes for easy copy and paste
  print("The best set of genes this time are:")
  print(noquote(rownames(thresh)))
  # return the Darwinian matrix
  cluster_by_gene
}

setwd('/home/sam/scRNAseq/FreqTables/')
data_path = 'ExpressingExpressionMats.RData'

key_genes <- c('Rbpms', 'Opn1sw','Grm6')#, 'Prkcq', 'Coch', 'Spp1', 'Meis2', 'Mafb', 'Penk', 'Opn4', 'Etv1', 'Zic1')


minimal_unique_clusters <- Darwin(data_path, cellClass = 'RGC', key_genes = key_genes,
                                  filtered = F, exponentiated = F,
                                  min_genes = 24,
                                  n_generations = 100, n_offspring = 1000, drop_factor = 5,
                                  verbose = T, lossFn = T, add_back = 50,)
