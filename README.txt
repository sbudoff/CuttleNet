# Order of Scripts to Reproduce Experiments

## Dataset Merging and Preprocessing - R

1) Create Frequency tables

/DatasetMergingProcessing/RGC_expression_mat_grouper.R

/DatasetMergingProcessing/AC_expression_mat_grouper.R

/DatasetMergingProcessing/BP_expression_mat_grouper.R

2) Identify expressing and non-ubiquitous genes

/DatasetMergingProcessing/IdentifyInformativeGenes.R

3) Create Merged expression matrix and metadta lists
/DatasetMergingProcessing/Network_gene_lists.R
## GraSP

### AUCPR

1) Prepare RData file containing AUCPR targets

/GraSP/AUCPR_subseter.R

2) Train and evaluate NNs with AUCPR target datasets

/GraSP/BarcodeNN_AUCPR_Comparison.ipynb

### ML Methods

Useage Note: multiple seeds {18, 36,54,72, 90} were set manually to run R script independently, portions of script commented out if iteration crashed and restarted below a save point. Timing determined by recording savepoint updates in system folder.

/GraSP/traditionalMLcomparisons.R

### GraSP

/GraSP/BarcodeNN_GradientSelection_ErrorBars.ipynb

### Analysis and Figures


#### Figure 1:

1) /Analysis/Network_graphs_rgc.R

2) /Analysis/Fig1b_allComparisons.R

#### Figure 2:

/Analysis/Fig2b_allComparisons.R


#### Figure 3:

1) /CuttleNet/Analysis/Exp1and2Merger.R

2) /CuttleNet/Analysis/Exp123Merger.R

3) /CuttleNet/Analysis/ExpALLMerger.R

4) /CuttleNet/Analysis/Fig3.R

#### Figure 4A:

/CuttleNet/Analysis/cell_counts.R

#### Figure 5A:

/AnalysisAUCPR_supFIg.R

#### Figure 6A:

/Analysis/Network_gene_lists.R

#### Figure 7A, 8A:

/Analysis/grid_search_best_condition_analysis.R

#### Figure 9A:

/Analysis/Experiment1_analysis.R
/Analysis/Experiment3_analysis_figures.R


### Clustering

/CuttleNet/UMAP_tSNE/clustering300genes.R

## CuttleNet

### CuttleNet as implemented for figure

/CuttleNet/BarcodeNN_Class2Subclass_CuttleNet.ipynb

### CuttleNet varient architectures example

/CuttleNet_Inference_XeniumSeg-LongArms.ipynb

## CellTICS

1) /Benchmark/CellTICS_converter.ipynb
2a) /Benchmark/MultiRunCellTICS.py
2b) /Benchmark/MultiRunCellTICS_ablation.py
3a) /Benchmark/CellTICS_result_compiler.R
3b) /Benchmark/CellTICS_result_compiler_ablation.R
