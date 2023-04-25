#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(argparse))
suppressMessages(library(ggraph)) ## even though this script does not use this package directly, it expects this to be loaded in the background for clonalNetwork function in scRepertoire package
suppressMessages(ucsc <- modules::use("/Users/tha8tf/Documents/GitHub/immunereceptoire-clones/ucsc.R"))
#suppressMessages(graphics <- modules::use("graphics.R"))
suppressMessages(scRep <- modules::use("/Users/tha8tf/Documents/GitHub/immunereceptoire-clones/screpertoire_onesample.R"))

D40_COLORS <- c("#FB1C0D", "#0DE400", "#0D00FF", "#E8B4BD", "#FD00EA", "#0DD1FE", "#FF9B0D", "#0D601C", "#C50D69", "#CACA16", "#722A91", "#00DEBF", "#863B00", "#5D7C91", "#FD84D8", "#C100FB", "#8499FC", "#FD6658", "#83D87A", "#968549", "#DEB6FB", "#832E60", "#A8CAB0", "#FE8F95", "#FE1CBB", "#DF7CF8", "#FF0078", "#F9B781", "#4D493B", "#1C5198", "#7C32CE", "#EFBC16", "#7CD2DE", "#B30DA7", "#9FC0F6", "#7A940D", "#9B0000", "#946D9B", "#C8C2D9", "#94605A")


# Useful links
# 1. https://github.com/ncborcherding/scRepertoire
# 2. https://ncborcherding.github.io/vignettes/vignette.html
# 3. https://ncborcherding.github.io/vignettes/beyond_screpertoire.html

set.seed(45)

create_combined_multi_sample_object <- function(args) {
  
  contig_list = purrr::map(args$contigfile,read.table,sep = ",",header = T,stringsAsFactors = F)
  
  combined_sample_object = scRep$create_combined_sample_object(contig_list=contig_list,samples =as.character(args$samples), immunetype= args$immunetype,removeNA=args$removeNA,removeMulti=args$removeMulti,filterMulti=args$filterMulti,tcells=args$tcells)
  
  return (combined_sample_object)
}


get_args <- function(){
  
  parser <- ArgumentParser(description="scRepertoire TCR/BCR sequencing data analysis")
  parser$add_argument(
    "--contigfile",
    help="Path to the cellranger generated filtered_contig_annotation.csv file",
    type="character", nargs="+"
  )
  parser$add_argument(
    "--seuratobject",
    help="Path to the seurat object that contains corresponding gene expression of the TCR/BCR sequencing in --contigfile argument. Note: Cluster analysis must be done either through Seurat or through an external program such that the seurat object contains atleast one set of clusters assigned to the cells.",
    type="character", required="True"
  )
  parser$add_argument(
    "--samples",
    help="name of the metadata header in seuratobject that indicates the clusters of interest. Default: 'seurat_clusters'. Can take more than 1 set of clusters",
    type="character", nargs="+"
  )
  parser$add_argument(
    "--cluster",
    help="name of the metadata header in seuratobject that indicates the clusters of interest. Default: 'seurat_clusters'. Can take more than 1 set of clusters",
    type="character", default="seurat_clusters", nargs="+"
  )
  parser$add_argument(
    "--immunetype",
    help="Is it TCR or BCR data or contains both type of immune receptors? Options are: 'TCR', 'BCR',or 'both'. ",
    required="True", choices=c("TCR","BCR","both"), type="character"
  )
  parser$add_argument(
    "--tcells", 
    help="Enrichment of which type of T cells if --immunetype is 'TCR' or 'both'. Options are 'alpha-beta','gamma-delta', or 'both'" ,
    required="True", choices=c("alpha-beta","gamma-delta","both"), type="character"
  )
  parser$add_argument(
    "--removeMulti",
    help="if set to 'True', this is a stringent filter to remove any cell barcode with more than 2 immune receptor chains. Default: False includes and incorporates  cells with > 2 chains", 
    action="store_true"
  )
  parser$add_argument(
    "--filterMulti",
    help="if set to 'True', it isolates the top 2 expressed chains in cell barcodes with multiple chains. Default: False includes and incorporates cells with > 2 chains", 
    action="store_true"
  )
  parser$add_argument(
    "--removeNA",
    help="if set to 'True', this is a stringent filter to remove any cell barcode with an NA value in at least one of the chains. Default: False includes and incorporates cells with 1 NA value.", 
    action="store_true"
  )
  parser$add_argument(
    "--clone",
    help="Definition of clone: Options are 'gene', 'nt','aa','strict'. 
    'gene' - use the VDJC genes comprising the TCR/Ig
    'nt' - use the nucleotide sequence of the CDR3 region
    'aa' - use the amino acid sequence of the CDR3 region
    'strict' - DEFAULT: use the VDJC genes comprising the TCR/Ig + the nucleotide sequence of the CDR3 region. This is the proper definition of clonotype. 
    note that the clonotype is called using the combination of genes or nt/aa CDR3 sequences for both loci.
    'gene' is most sensitive, 'nt' or 'aa' moderately so, and 'strict' the most specific.", 
    choices=c("gene","nt","aa","strict"), type="character", default="strict", nargs="+"
  )
  
  args = parser$parse_args(commandArgs(trailingOnly = TRUE))
  return (args)
}

# Parse arguments
args <- get_args()

print("Used parameters")
print(args)

print("Loading contigs")
combined_sample_object <- create_combined_multi_sample_object(args)

print("Loading seurat object")
vdjc_seurat_object <- scRep$create_vdjc_seurat_object(args$seuratobject,args$cluster,combined_sample_object)

print("Creating directory to store files")
dir.create("scRepertoire_output_multisample_tcr_multisample_try_5")
setwd("./scRepertoire_output_multisample_tcr_multisample_try_5")

print("Visualizing VDJC Genes")
#visualize_vdjc_genes_multi_sample(args,combined_sample_object)
scRep$visualize_vdjc_genes(combined_sample_object = combined_sample_object,immunetype = args$immunetype)

print("Visualizing Unique Clonotypes")
scRep$visualize_bcr_tcr(args$clone,args$cluster,vdjc_seurat_object,combined_sample_object,args$immunetype)

saveRDS(vdjc_seurat_object,"vdjc_seurat_object.rds")

print("Visualizing Clonal Homeostasis")
vdjc_seurat_object = scRep$visualize_clonalhomeostasis(args$clone,args$cluster,vdjc_seurat_object,combined_sample_object)

print("Visualizing Clonal Overlap")
scRep$visualize_clonaloverlap(args$clone, args$cluster, vdjc_seurat_object)

print("Visualizing Clonal Network")
scRep$visualize_clonalnetwork(args$clone, args$cluster, vdjc_seurat_object, combined_sample_object)

print("Visualizing Chord Diagrams")
scRep$visualize_chord_diagrams(args$clone, args$cluster, vdjc_seurat_object, combined_sample_object)  

print("Visualizing Clonal Diversity")
scRep$visualize_clonaldiversity(args$clone, args$cluster, vdjc_seurat_object)

print("Saving SeuratObject with scRepertoire metadata")
saveRDS(vdjc_seurat_object,"vdjc_seurat_object.rds")

print("Analyses Done!")

