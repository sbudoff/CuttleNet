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

Darwin <- function(data_path, key_genes = c(), gene_subset = c(),
                   filtered = T, exponentiated = T, SNR = F, expressingOnly = T,
                   expression_threshold = 0.05, min_expression_var = 0.2,
                   ubiquitous_threshold = 0.5, min_expressing_var = 0.2,
                   min_genes = 9, add_back = 10,
                   n_generations = 100, n_offspring = 100, drop_factor = 10,
                   verbose = T, lossFn = T) {
  # Load the cell class by gene matrices
  load(data_path)
  # Rename the RGC workspace used in testing
  if (exists("RGC_cluster_expressing")){
    cluster_expressing = RGC_cluster_expressing
    cluster_expression = RGC_cluster_expression
    cluster_SD = RGC_cluster_SD
    cluster_median = RGC_cluster_median
  }

  if (length(gene_subset) > min_genes) {
    cluster_expression <- cluster_expression[gene_subset,-1]
    cluster_expressing <- cluster_expressing[gene_subset,-1]
    cluster_median <- cluster_median[gene_subset,-1]
    cluster_SD <- cluster_SD[gene_subset,-1]
  } else {
    cluster_expression <- cluster_expression[,-1]
    cluster_expressing <- cluster_expressing[,-1]
    cluster_median <- cluster_median[,-1]
    cluster_SD <- cluster_SD[,-1]
  }
  
  if (filtered) {
    # Remove genes that have no information content
    DropIndex <- which(rowMaxs(as.matrix(cluster_expressing)) < ubiquitous_threshold)
    DropIndex <- unique(c(DropIndex, which(rowMaxs(as.matrix(cluster_expression)) < expression_threshold)))
    DropIndex <- unique(c(DropIndex, which(rowSds(as.matrix(cluster_expressing)) < min_expressing_var)))
    DropIndex <- unique(c(DropIndex, which(rowSds(as.matrix(cluster_expression)) < min_expression_var)))
    # DropIndex <- unique(c(DropIndex, RedundantRows(cluster_expressing[,-1], precision = 1)))
    
    cluster_expression <- cluster_expression[-DropIndex,]
    cluster_expressing <- cluster_expressing[-DropIndex,]
    cluster_median <- cluster_median[-DropIndex,]
    cluster_SD <- cluster_SD[-DropIndex,]
  }
  
  
  
  # Transform matrices to create fitness matrix
  
  # COV_mat <- cluster_SD[,-1] / cluster_expression[,-1]
  if (exponentiated){
    cluster_expression <- exp(cluster_expression)
    SNR_mat <- exp(cluster_expression / (cluster_SD+.01))
    cluster_expressing <- exp(cluster_expressing)
  } else {
    cluster_expression <- cluster_expression
    cluster_expressing <- cluster_expressing
    SNR_mat <- cluster_expression / (cluster_SD+.01)
  }
  if (expressingOnly) 
    {
    thresh <- cluster_expressing
  } else {
    if (SNR){
      thresh <- t(scale(t(SNR_mat)) * scale(t(cluster_expressing)))
    } else {
      thresh <- t(scale(t(cluster_expression)) * scale(t(cluster_expressing)))
    }
  }
  
  # print(colSums(t(scale(t(cluster_expression)))))
  # # Threshold the genes based on joint minimum filters
  # thresh <- select(cluster_expressing, -GENE) >= expressing_threshold
  # expr_thresh <- select(cluster_expression, -GENE) >= expression_threshold
  # thresh <- ((thresh + expr_thresh) == 2) +0
  
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
  cluster_by_gene <- cluster_median %>% 
    rownames_to_column('GENE') %>%
    filter(GENE %in% rownames(thresh)) %>%
    as.data.frame() %>%
    gather(Class, median, -GENE) %>%
    merge(cluster_by_gene, by = c("GENE", "Class"))
  cluster_by_gene <- cluster_SD %>% 
    rownames_to_column('GENE') %>%
    filter(GENE %in% rownames(thresh)) %>%
    as.data.frame() %>%
    gather(Class, SD, -GENE) %>%
    merge(cluster_by_gene, by = c("GENE", "Class")) %>%
    arrange(GENE)
  
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

setwd('/home/sam/scRNAseq/FreqTables/')
load('ExpressingExpressionMats.RData') 

setwd("/home/sam/scRNAseq/SCP509/cluster/")
data_path = "RGC_cluster_ExpressionMats.RData"
# key_genes <- c('Rbpms', 'Opn1sw','Grm6','Spp1','Tbr1','Pcdh20','Opn4','Il1rapl2','Adcyap1','Mmp17','Col25a1','Gpr88')
# 
# 
# key_genes <- c('Spp1','Opn4', 'Igfbp2', 'Tbx20', 'Igf1', 
#                 'Foxp2', 'Ctxn3','Penk', 'Cck')
# 
# key_genes <- c('Rbpms', 'Opn1sw','Grm6', 'Igf1', 'Meis2')#  'Pvalb', 'Igf1', 'Meis2', 'Spp1')
               

# unique_clusters <- Darwin(data_path, key_genes = key_genes, gene_subset = c(),
#                             filtered = T, exponentiated = F, SNR = F, expressingOnly = F,
#                             expression_threshold = 0, min_expression_var = 0,
#                             ubiquitous_threshold = 0, min_expressing_var = 0.18,
#                             min_genes = 12, add_back = 10,
#                             n_generations = 200, n_offspring = 100, drop_factor = 5,
#                             verbose = T, lossFn = F)


# key_genes <- c('Rbpms', 'Opn1sw','Grm6','Igf1','Zic1','Spp1','Penk','Lmo2','Meis2','Etv1','Opn4')
# gene_subset <- c('Rbpms', 'Opn1sw','Grm6',
#                  'Tbr1','Pcdh20','Il1rapl2','Adcyap1','Mmp17','Col25a1','Gpr88',
#                  'Igfbp2', 'Ppp1r1c', 'Parm1', 'Neurod2','S100b', 
#                  'Spon1', 'Chl1', 'Prkcq', 'Foxp2', 'Ctxn3', 'Irx4', 'Gal', 'Fxyd6',  
#                  'Zeb2', 'Tfap2d', 'Chrna3', 'RP23-407N2.2', 'Dkk3', 'Cck', 'Trp53i11', 'Irx3',
#                  'Cdkn1c', 'Pcdh10', 'Khdrbs2', 'Necab2', 'Ndnf', 'Isl2', 'Pcdh11x', 'Crhbp', 
#                  'Pcdh9', 'Kcnip4', 'Bex1', 'Rprm', 'Ppp1r17','Evc2', 'Synpr', 'Mt1', 'Pcp4l1',
#                  'Ramp3', 'Lpl','Lypd6','Syt2', 'Coch','Eomes', 'Irx6', 'Cd83', 'Gabra2', 
#                  'Reln', 'Slc6a1', 'Sorl1', 'Kcnb2', 'Kcnab3', 'Marcksl1', 'Gprc5b', 'Six6',
#                  'Sh3bgr', 'Atp2b4', 'Fgf1', 'Pou6f2', 'Kitl', 'Clstn2', 'Kcna1', 'Kcnd2', 
#                  'Lmo1', 'Zfh3x', 'Sv2c', 'Pvalb', 'Ntng1', 'Lmo4', 'Kcnip1', 'Isl1', 
#                  'Dlgap1', 'Ctxn3', 'Cpne4', 'Calb1', 'BC048546', 'Zmat4', 'Syt2', 'Stk32a',
#                  'Gnai1', 'Dlg2', 'Cplx2', 'Vamp1', 'Slc24a2', 'Satb1', 'Pcsk2', 'Gabrg3',
#                  'Fgf13', 'Cacng3', 'Alcam', 'Prph', 'Nrxn3', 'Mgat4c', 'Bdnf', 'Mafb',
#                  'Pou4f1', 'Junb', 'Magi1', 'Syndig1l', 'Stk32c', 'Car10', '6330403K07Rik',
#                   'Zfhx3', 'Vgf', 'Nrgn', 'Lgals1', 'Ebf3', 'Diras2', 'Ano3',
#                  '2510009E07Rik', 'Sema5a', 'Myo1b', 'Ly6h', 'Sptssb', 'Ndrg2', 'Lxn',
#                  'Chrnb3')

key_genes <- c('Rbpms', 'Opn1sw','Grm6')#, 'Prkcq', 'Coch', 'Spp1', 'Meis2', 'Mafb', 'Penk', 'Opn4', 'Etv1', 'Zic1')
# gene_subset <- c('Rbpms', 'Opn1sw','Grm6',
#   'Crhbp', 'Cnr1','Scn4b','Eomes','Neurod2','Foxp2','Kctd4','Lypd1','Igfbp2',
#   'Prkcq','Coch', 'Spp1', 'Meis2', 'Penk', 'Lmo2', 'Opn4', 'Cd24a','Etv1',
#   'Zic1', 'Pvalb', 
#   'Igf1', 'Mafb','Pcdh20','Pcdh9','Pcdh7',
#   'Igfbp2', 'Ppp1r1c', 'Parm1', 'Neurod2','S100b', 
#   'Lmo1', 'Zfh3x', 'Sv2c', 'Pvalb', 'Ntng1', 'Lmo4', 'Kcnip1', 'Isl1', 
#   'Tbr1','Pcdh20','Il1rapl2','Adcyap1','Mmp17','Col25a1','Gpr88', 'Tac1',
#   'Gnai1', 'Dlg2', 'Cplx2', 'Vamp1', 'Slc24a2', 'Satb1', 'Pcsk2', 'Gabrg3', 'Jamb'
# )

useful_genes = c("Ano3",      "Anxa2",     "Atp2b4",    "Brinp1",    "Cacng3",    "Calb1",     "Cck",       "Cd24a",     "Cnr1",     
                "Coch",      "Crhbp",     "Crip2",     "Ctxn3",     "Dkk3",      "Dlgap1",    "Ebf3",      "Eomes",     "Etv1",     
                "Fgf1",      "Foxp2",     "Gabra1",    "Gabrb1",    "Gabrg3",    "Gal",       "Gm2115",    "Gprc5b",    "Grin2a",   
                "Gucy1a3",   "Hcn1",      "Igf1",      "Igfbp2",    "Irx6",      "Isl1",      "Kcna1",     "Kcnb2",     "Kcnd2",    
                "Kcnh2",     "Kcnip1",    "Kcnip2",    "Kitl",      "Lgi3",      "Lmo2",      "Lor",       "Lxn",       "Lypd1",    
                "Lypd6",     "Mafb",      "Marcksl1",  "Meis2",     "Mgat4c",    "Mmp9",      "Ndnf",      "Neurod2",   "Nr2f2",    
                "Nrgn",      "Nrxn3",     "Ntng1",     "Opn4",      "Pcdh10",   "Pcdh11x",   "Pcdh17",    "Pcdh7",     "Pcdh9",    
                "Penk",      "Pou4f3",    "Pou6f2",    "Ppargc1a",  "Ppp1r17",   "Prkcq",     "Ramp3",     "Rgs8",      "Rprm",     
                "Sema5a",    "Serpinb1b", "Six3",      "Slc6a1",    "Sphkap",    "Spon1",     "Spp1",      "Sptssb",    "Sv2b",     
                "Sv2c",      "Syndig1l",  "Syt13",     "Syt2",      "Syt5",      "Syt6",      "Tac1",      "Trp53i11",  "Tusc5",    
                "Vgf",       "Zfhx3",     "Zic1",      "Zmat4")

minimal_unique_clusters <- Darwin(data_path, key_genes = key_genes, gene_subset = useful_genes,
                                  filtered = F, exponentiated = F, SNR = F,
                                  min_genes = 24,
                                  n_generations = 75, n_offspring = 1000, drop_factor = 2,
                                  verbose = T, lossFn = F)



RNAscope <- minimal_unique_clusters %>%
  filter(GENE %in% key_genes)

unique(RNAscope$GENE)

write.csv(RNAscope,'RNAscope_12plex_expression_atlas.csv')

# key_genes_ind <- which(rownames(RGC_cluster_expressing) %in% key_genes)
# thresh <- RGC_cluster_expressing[key_genes_ind,-1]
# barcoder(thresh > .7, verbose = T)
# load(data_path)
# SNR_mat <- RGC_cluster_expression[,-1] / (RGC_cluster_SD[,-1]+.01)
# COV_mat <- RGC_cluster_SD[,-1] / RGC_cluster_expression[,-1]
# exp_mean <- exp(RGC_cluster_expression[,-1])
# exp_SNR <- exp(SNR_mat)
# exp_expressing <- exp(RGC_cluster_expressing[,-1])
# fitness <- exp_SNR * exp_expressing
# fitness <- t(scale(t(fitness)))
# distance_composite sum(rowMads(as.matrix(fitness[key_genes,])))

# test <- minimal_unique_clusters %>%
#   group_by(GENE) %>%
#   summarise(score = (max(mean) * max(expressing)) / (mean(mean) * mean(expressing))) %>%
#   arrange(desc(GENE))
# 
# 
# # Load the cell class by gene matrices
# load(data_path)
# 
# 
# test <- RGC_cluster_expressing[c(1:100),-1]
# 
# print(sqrt(sum(test[,1.] - test[,2])**2))
# dist(t(test))
# load(data_path)
# 
# 
# gene_subset <- c('Rbpms', 'Opn1sw', 'Grm6', 'Pcp4', 'Igf1','Zic1','Cartpt','Spp1','Penk','Lmo2','Meis2','Etv1','Opn4')
# test <- na.omit(select(RGC_cluster_expressing[gene_subset,],-GENE))



# barcoder(test, verbose = T, lossFn = F)
# l = nDELTArow2col2(test)
# gene_subset <- c('Rbpms', 'Opn1sw', 'Grm6', 'Pvalb', 'Pcdh9', 'Kcnip4', 'Pcdh7', 'Gnai1', 'Dlg2', 'Cplx2', 'Vamp1', 'Slc24a2', 'Satb1')
# test2 <- select(RGC_cluster_expressing[gene_subset,],-GENE)
# l2 = nDELTArow2col2(test2)
# 
# barcoder(test2, verbose = T, lossFn = F)
# 
# print(l2<l)
