library(tidyverse)
library(matrixStats)
library(Matrix)

barcoder <- function(class_mat, verbose = F){
  # This Function creates a 'barcode' for each vector of genes defining a class. 
  # This barcode is simply a square matrix of booleans that checks each class against every other class to see if they share a vector of genes.
  # If verbose is true, a visualization of the barcode will be displayed
  n_classes = nrow(t(class_mat))
  n_genes = nrow(class_mat)
  barcode_compare <- data.frame(matrix(ncol = n_classes, nrow = n_classes))
  rownames(barcode_compare) <- colnames(class_mat)
  colnames(barcode_compare) <- colnames(class_mat)
  for (i in 1:n_classes) {
    barcode_i = class_mat[,i]
    for (j in 1:n_classes) {
      barcode_j = class_mat[,j]
      barcode_compare[i,j] = sum(barcode_i == barcode_j) == nrow(class_mat)
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

Darwin <- function(data_path, key_genes = c(), 
                   expression_threshold = 0.2,
                   expressing_threshold = 0.9,
                   min_expression_var = 0.2,
                   min_expressing_var = 0.2,
                   n_generations = 100, n_offspring = 100, drop_factor = 10,
                   verbose = F) {
  # Load the cell class by gene matrices
  load(data_path)
  # Rename the RGC workspace used in testing
  if (exists("RGC_cluster_expressing")){
    cluster_expressing = RGC_cluster_expressing
    cluster_expression = RGC_cluster_expression
    cluster_SD = RGC_cluster_SD
    cluster_median = RGC_cluster_median
  }
  # Threshold the genes based on joint minimum filters
  thresh <- select(cluster_expressing, -GENE) >= expressing_threshold
  expr_thresh <- select(cluster_expression, -GENE) >= expression_threshold
  thresh <- ((thresh + expr_thresh) == 2) +0
  # Get number of classes we are dealing with
  n_classes = nrow(t(thresh))
  # Remove genes that have no information content
  DropIndex <- which(rowSums(thresh) == 0)
  DropIndex <- c(DropIndex, which(rowSums(thresh) == n_classes))
  # Remove genes based on intra-class vector variance thresholds
  DropIndex <- c(DropIndex, which(rowSds(as.matrix(cluster_expressing[,-1]))<min_expressing_var))
  DropIndex <- c(DropIndex, which(rowSds(as.matrix(cluster_expression[,-1]))<min_expression_var))
  thresh <- thresh[-DropIndex,]
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
  generations = 1:n_generations
  kids = 1:n_offspring
  genes_dropped = c()
  genes_remaining = c()
  for (gen in generations) {
    n_genes = nrow(thresh) # Recalculate number of genes in this new generation
    n_genes_min = n_genes  # Reset the minimum number of genes an offspring in the new generation has found
    for (i in kids) {
      n_drop = sample.int(n_genes/drop_factor, 1) # Pick a random number of genes to drop
      drop = sample(n_start:n_genes, n_drop) # Create a random list of genes to drop
      thresh_i = thresh[-drop,] # Drop the genes
      if (sum(colSums(thresh_i) > 0) == n_classes) { # Check that each offspring class is defined by at least one gene
        n_genes_i = nrow(thresh_i) # Compute this offspring's number of genes
        if (n_genes_i < n_genes_min) { # See if this offspring is the most fit, that is has the smallest set of genes
          # Store fitness values of the best offspring
          n_genes_min = n_genes_i
          drop_best = drop
        }
      }
    }
    # Run the barcoder checksum on the fittest offspring and see if it gets to make the next generation
    if (length(drop_best) > 0) {
      thresh_i = thresh[-drop_best,]
      print(sprintf("Dropped %s genes in generation %s, %s genes remain", length(drop_best), gen, nrow(thresh_i)))
      genes_dropped = append(genes_dropped, length(drop_best))
      genes_remaining = append(genes_remaining, nrow(thresh_i))
      checkSum = barcoder(thresh_i, verbose = verbose)
      if (checkSum == n_classes) {
        thresh = thresh_i
      }
      else {
        print("Undifferentiable gene sets identified, generation discarded.")
      }
    } else {
      drop_best = NA
      print(sprintf("No genes dropped in generation %s, %s genes remain", gen, nrow(thresh)))}
  }
  # Visualize the removal of genes over the generations
  gen_plot <- data_frame(generations, genes_dropped, genes_remaining) %>%
    ggplot(aes(x = generations)) +
      geom_line(aes(y = genes_dropped), color = 'red') +
      geom_line(aes(y = genes_remaining), color = 'green')
  print(gen_plot)
  # Reorganize the input data based on the thresholding such that all matrices are combined
  cluster_by_gene <- cluster_expression %>% 
    filter(GENE %in% rownames(thresh)) %>%
    as.data.frame() %>% 
    gather(Class, mean, -GENE) 
  cluster_by_gene <- cluster_expressing %>% 
    filter(GENE %in% rownames(thresh)) %>%
    as.data.frame() %>% 
    gather(Class, expressing, -GENE) %>%
    merge(cluster_by_gene, by = c("GENE", "Class"))
  cluster_by_gene <- cluster_median %>% 
    filter(GENE %in% rownames(thresh)) %>%
    as.data.frame() %>% 
    gather(Class, median, -GENE) %>%
    merge(cluster_by_gene, by = c("GENE", "Class"))
  cluster_by_gene <- cluster_SD %>% 
    filter(GENE %in% rownames(thresh)) %>%
    as.data.frame() %>% 
    gather(Class, SD, -GENE) %>%
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
  # Print the minimal set of genes for easy copy and paste
  print("The best set of genes this time are:")
  print(noquote(rownames(thresh)))
  # return the Darwinian matrix
  cluster_by_gene
}

setwd("/home/sam/scRNAseq/SCP509/cluster/")
data_path = "RGC_cluster_ExpressionMats.RData"
key_genes <- c('Rbpms', 'Opn1sw','Grm6','Spp1','Tbr1','Pcdh20','Opn4','Il1rapl2','Adcyap1','Mmp17','Col25a1','Gpr88')

key_genes <- c('Spp1','Opn4', 'Igf1', 'Igfbp2', 'Ppp1r1c', 'Parm1', 'Neurod2', 
               'Spon1', 'Chl1', 'Prkcq', 'Foxp2', 'Ctxn3', 'Irx4', 'Penk', 'Gal', 
               'Zeb2', 'Tfap2d', 'Chrna3', 'RP23-407N2.2', 'Dkk3')

key_genes <- c('Spp1','Opn4', 'Igfbp2', 'Tbx20',
               'Foxp2', 'Ctxn3','Penk')

# min_RBPMS <- min(cluster_expressing['Rbpms',-1])
# min_RBPMS_expr <- min(cluster_expression['Rbpms',-1])

minimal_unique_clusters <- Darwin(data_path, key_genes = key_genes, 
                                  expression_threshold = 0.2,
                                  expressing_threshold = 0.7,
                                  min_expression_var = 0.2,
                                  min_expressing_var = 0.2,
                                  n_generations = 100, n_offspring = 100, drop_factor = 15,
                                  verbose = F) 





key_genes_ind <- which(rownames(RGC_cluster_expressing) %in% key_genes)
thresh <- RGC_cluster_expressing[key_genes_ind,-1]
barcoder(thresh > .7, verbose = T)
