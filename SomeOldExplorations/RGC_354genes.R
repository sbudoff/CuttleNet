library(tidyverse)

df_root = '/home/sam/Classes/Stats/Bayes/assignments/'
raw_seq_path='/home/sam/scRNAseq/SCP509/expression/RGC_Atlas.csv'
raw_clust_path="/home/sam/scRNAseq/SCP509/cluster/RGC_Atlas_coordinates.txt"
save_path = '/home/sam/scRNAseq/RNAscope/RNAscope/RGC_df354.Rdata'

expMatrix <- read.csv(raw_seq_path, head = TRUE) #, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster"))
rownames(expMatrix) <- expMatrix$GENE

# Open and find the cluster specific sample IDs
clusters <- read.csv(raw_clust_path, head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster", "Experiment","SampleID", "BatchID","Background"))
clusters <- clusters %>%
  select(ID, Cluster) %>%
  mutate(ID = gsub('-','.', ID))
RGC_cluster_expression <- select(expMatrix, GENE)
rownames(RGC_cluster_expression) <- RGC_cluster_expression$GENE

# Load list of genes
load(file='/home/sam/scRNAseq/RNAscope/RNAscope/gene_set_big.Rdata')

df <- expMatrix %>%
    filter(GENE %in% genes) %>%
    select(-GENE) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('ID') %>%
    left_join(clusters, by='ID') %>%
    select(-ID) %>%
    mutate(Cluster = as.factor(Cluster))

key_genes = sub("4833423E24Rik" , "four833423E24Rik", genes)
key_genes = sub("6330403K07Rik" , "six330403K07Rik", key_genes)
key_genes = sub("6430548M08Rik" , "six430548M08Rik", key_genes)
key_genes = sub("2510009E07Rik" , "two510009E07Rik", key_genes)
df <- df %>%
  mutate(four833423E24Rik = `4833423E24Rik`,
         six330403K07Rik = `6330403K07Rik`,
         six430548M08Rik = `6430548M08Rik`,
         two510009E07Rik = `2510009E07Rik`
         ) %>%
  select(-`4833423E24Rik`, -`6330403K07Rik`, -`6430548M08Rik`, -`2510009E07Rik`)

col_order <- unique(c("Cluster", "two510009E07Rik", "four833423E24Rik", 
                      "six330403K07Rik", "six430548M08Rik", sort(key_genes)))
# Remove genes not found in dataset
col_order <- col_order[!col_order == col_order[!col_order %in% names(df)]]

df <- df[, col_order]

# Manually enter in the Tran proposed classification scheme
# Proposed marker sets from Tran (2019) Figure 1F
h = list()
h[1] <- data.frame(c(Cluster="1_W3D1.1", Gene1='Serpine2', Gene1_pos=T,Gene2='Amigo2',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[2] <- data.frame(c(Cluster="2_W3D1.2",Gene1='Lypd1',Gene1_pos=T,Gene2='Ntrk1', Gene2_pos=F, Gene3=NA,Gene3_pos=NA))
h[3] <- data.frame(c(Cluster="3_FminiON",Gene1='Foxp2',Gene1_pos=T,Gene2='Irx4',Gene2_pos=T,Gene3=NA, Gene3_pos=NA))
h[4] <- data.frame(c(Cluster="4_FminiOFF",Gene1='Pde1a',Gene1_pos=T,Gene2=NA,Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
h[5] <- data.frame(c(Cluster="5_J-RGC",Gene1='Tbr1',Gene1_pos=T,Gene2='Pcdh20',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
h[6] <- data.frame(c(Cluster="6_W3B",Gene1='Zic1',Gene1_pos=T,Gene2=NA,Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
h[7] <- data.frame(c(Cluster="7_Novel",Gene1='Tbx20',Gene1_pos=T,Gene2='Tagln2',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[8] <- data.frame(c(Cluster="8_Novel",Gene1='Prkcq',Gene1_pos=T,Gene2='Tac1',Gene2_pos=T, Gene3='Spp1',Gene3_pos=F))
h[9] <- data.frame(c(Cluster="9_Tbr1_Novel",Gene1='Slc7a11',Gene1_pos=T,Gene2='Plpp4',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[10] <- data.frame(c(Cluster="10_Novel",Gene1='Gpr88',Gene1_pos=T,Gene2=NA,Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
h[11] <- data.frame(c(Cluster="11_Novel", Gene1='Serpinb1b',Gene1_pos=T,Gene2='Gm17750',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
h[12] <- data.frame(c(Cluster="12_ooDS_NT",Gene1='Mmp17',Gene1_pos=T, Gene2=NA, Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
h[13] <- data.frame(c(Cluster="13_Novel", Gene1='Lypd1',Gene1_pos=T,Gene2='Ntrk1',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[14] <- data.frame(c(Cluster="14_ooDS_Cck",Gene1='Cartpt',Gene1_pos=T,Gene2='Vit',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[15] <- data.frame(c(Cluster="15_Novel",Gene1='Apela',Gene1_pos=T,Gene2=NA,Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
h[16] <- data.frame(c(Cluster="16_ooDS_DV",Gene1='Cartpt', Gene1_pos=T,Gene2='Col25a1',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[17] <- data.frame(c(Cluster="17_Tbr1_S1",Gene1='Tbr1',Gene1_pos=T,Gene2='Irx4',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
h[18] <- data.frame(c(Cluster="18_Novel",Gene1='Pcdh20',Gene1_pos=T,Gene2='4833423E24Rik',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[19] <- data.frame(c(Cluster="19_Novel",Gene1='Penk',Gene1_pos=T,Gene2='Prdm8', Gene2_pos=T,Gene3='Slc24a2', Gene3_pos=T))
h[20] <- data.frame(c(Cluster="20_Novel",Gene1='Penk',Gene1_pos=T, Gene2='Gal',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
h[21] <- data.frame(c(Cluster="21_Tbr1_S2",Gene1='Tbr1',Gene1_pos=T, Gene2='Calca', Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[22] <- data.frame(c(Cluster="22_M5",Gene1='Serpine2',Gene1_pos=T, Gene2='Cdhr1', Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
h[23] <- data.frame(c(Cluster="23_W3D2", Gene1='Prokr1',Gene1_pos=T, Gene2=NA,Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
h[24] <- data.frame(c(Cluster="24_Novel",Gene1='Fam19a4',Gene1_pos=T,Gene2=NA,Gene2_pos=NA,Gene3=NA, Gene3_pos=NA))
h[25] <- data.frame(c(Cluster="25_Novel",Gene1='Slc17a7',Gene1_pos=T, Gene2=NA, Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
h[26] <- data.frame(c(Cluster="26_Novel",Gene1='Penk',Gene1_pos=T,Gene2='Igfbp5',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
h[27] <- data.frame(c(Cluster="27_Novel",Gene1='Prkcq', Gene1_pos=T,Gene2=NA,Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
h[28] <- data.frame(c(Cluster="28_FmidiOFF",Gene1='Foxp2',Gene1_pos=T, Gene2='Cdk15',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[29] <- data.frame(c(Cluster="29_Novel",Gene1='Stxbp6', Gene1_pos=T, Gene2='Prlr', Gene2_pos=T,   Gene3=NA, Gene3_pos=NA))
h[30] <- data.frame(c(Cluster="30_Novel",Gene1='Postn',Gene1_pos=T,Gene2=NA, Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
h[31] <- data.frame(c(Cluster="31_M2",Gene1='Tbx20', Gene1_pos=T,Gene2='Spp1',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[32] <- data.frame(c(Cluster="32_F_Novel",Gene1='Rhox5',Gene1_pos=T, Gene2=NA, Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
h[33] <- data.frame(c(Cluster="33_M1",Gene1='Adcyap1', Gene1_pos=T,Gene2='Opn4',Gene2_pos=T, Gene3='Nmb',  Gene3_pos=F))
h[34] <- data.frame(c(Cluster="34_Novel",Gene1='Tpbg',  Gene1_pos=T,Gene2='Spp1', Gene2_pos=F,Gene3=NA, Gene3_pos=NA))
h[35] <- data.frame(c(Cluster="35_Novel",Gene1='Igfbp4',Gene1_pos=T, Gene2='Chrm2', Gene2_pos=T,Gene3=NA, Gene3_pos=NA))
h[36] <- data.frame(c(Cluster="36_Novel",Gene1='Stxbp6',Gene1_pos=T, Gene2='Coch', Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
h[37] <- data.frame(c(Cluster="37_Novel",Gene1='Ceacam10', Gene1_pos=T, Gene2=NA,Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
h[38] <- data.frame(c(Cluster="38_FmidiON", Gene1='Foxp2',Gene1_pos=T,Gene2='Anxa3', Gene2_pos=T,Gene3=NA, Gene3_pos=NA))
h[39] <- data.frame(c(Cluster="39_Novel", Gene1='Neurod2',Gene1_pos=T, Gene2='S100b',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[40] <- data.frame(c(Cluster="40_M1dup",Gene1='Nmb',Gene1_pos=T, Gene2=NA,Gene2_pos=NA, Gene3=NA, Gene3_pos=NA))
h[41] <- data.frame(c(Cluster="41_AlphaONT", Gene1='Spp1',Gene1_pos=T, Gene2='Kit',Gene2_pos=T, Gene3=NA, Gene3_pos=NA))
h[42] <- data.frame(c(Cluster="42_AlphaOFFS",Gene1='Spp1',Gene1_pos=T,Gene2='Fes',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[43] <- data.frame(c(Cluster="43_AlphaONS", Gene1='Spp1', Gene1_pos=T, Gene2='Il1rapl2', Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
h[44] <- data.frame(c(Cluster="44_Novel",Gene1='Bhlhe22', Gene1_pos=T,Gene2='Fxyd6',Gene2_pos=T,Gene3=NA, Gene3_pos=NA))
h[45] <- data.frame(c(Cluster="45_AlphaOFFT",Gene1='Spp1', Gene1_pos=T,Gene2='Tpbg',Gene2_pos=T,Gene3=NA, Gene3_pos=NA))

# Reorganize the manually entered genes into a dataframe
hypothesis <- t(as.data.frame(h))
rownames(hypothesis) <- NULL
hypothesis <- as.data.frame(hypothesis)
colnames(hypothesis) <- c("Cluster", 'Gene1', 'Gene1_pos', 'Gene2', 'Gene2_pos', 'Gene3', 'Gene3_pos')

# Save the list of frequency dataframes for easy loading elsewhere
save(df, hypothesis, key_genes, file=save_path)
