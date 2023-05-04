cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: screpetoire:v1.7.1


inputs:

  seurat_object:
    type: File
    inputBinding:
      prefix: "--seuratobject"
    doc: |
      Path to the seurat object that contains corresponding gene expression of the 
      TCR/BCR sequencing in --contigfile argument. Note: Cluster analysis must be done 
      either through Seurat or through an external program such that the seurat object 
      contains atleast one set of clusters assigned to the cells.

  filtered_contig_annotations_file:
    type: File
    inputBinding:
      prefix: "--contigfile"
    doc: |
      Path to the cellranger generated filtered_contig_annotation.csv file

  sample_names:
    - type: 
      - "null"
      - string
      - string[]
    inputBinding:
      prefix: "--samples"
    doc: |
      name of the samples in seuratobject that indicates the samples to analyze. 
      Can take more than 1 sample names

  cluster_choice:
  	- type: 
      - "null"
      - string
      - string[]
    inputBinding:
      prefix: "--cluster"
    doc: |
      name of the metadata header in seuratobject that indicates the clusters of 
      interest. Default: 'seurat_clusters'. Can take more than 1 set of clusters

  immune_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "TCR"
      - "BCR"
      - "both"
    inputBinding:
      prefix: "--immunetype"
    doc: |
      Is it TCR or BCR data or contains both type of immune receptors? 
      Options are: 'TCR', 'BCR',or 'both'.

  type_of_t_cells:
	type:
    - "null"
    - type: enum
      symbols:
      - "alpha-beta"
      - "gamma-delta"
      - "both"
    inputBinding:
      prefix: "--tcells"
    doc: |
      Enrichment of which type of T cells if --immunetype is 'TCR' or 'both'. 
      Options are 'alpha-beta','gamma-delta', or 'both'"

  remove_barcodes_with_multiple_chains:
    type: boolean?
    inputBinding:
      prefix: "--removeMulti"
    doc: |
      if set to 'True', this is a stringent filter to remove any cell barcode 
      with more than 2 immune receptor chains. Default: False includes and incorporates 
      cells with > 2 chains
      
  isolate_top_chains:
    type: boolean?
    inputBinding:
      prefix: "--filterMulti"
    doc: |
      if set to 'True', it isolates the top 2 expressed chains in cell barcodes with 
      multiple chains. Default: False includes and incorporates cells with > 2 chains

  remove_cells_with_missing_chain_info:
    type: boolean?
    inputBinding:
      prefix: "--removeNA"
    doc: |
      if set to 'True', this is a stringent filter to remove any cell barcode with 
      an NA value in at least one of the chains. Default: False includes and incorporates 
      cells with 1 NA value

  clone_definition:
    type: 
    - "null"
    - type: enum
      symbols:
      - "gene"
      - "nt"
      - "aa"
      - "strict"
    inputBinding:
      prefix: "--clone"
    doc: |
      Definition of clone: Options are 'gene', 'nt','aa','strict'. 
      'gene' - use the VDJC genes comprising the TCR/Ig
      'nt' - use the nucleotide sequence of the CDR3 region
      'aa' - use the amino acid sequence of the CDR3 region
      'strict' - DEFAULT: use the VDJC genes comprising the TCR/Ig + the nucleotide 
      sequence of the CDR3 region. This is the proper definition of clonotype. 
      note that the clonotype is called using the combination of genes or nt/aa CDR3 
      sequences for both loci. 'gene' is most sensitive, 'nt' or 'aa' moderately so, 
      and 'strict' the most specific.


outputs:

  umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_res_*.png"
    doc: |
      Clustered cells UMAP.
      PNG format

  umap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_res_*.pdf"
    doc: |
      Clustered cells UMAP.
      PDF format

  slh_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_slh_res_*.png"
    doc: |
      Silhouette scores. Downsampled to max 500 cells per cluster.
      PNG format

  slh_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_slh_res_*.pdf"
    doc: |
      Silhouette scores. Downsampled to max 500 cells per cluster.
      PDF format

  umap_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_idnt_res_*.png"
    doc: |
      Split by dataset clustered cells UMAP.
      PNG format

  umap_spl_idnt_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_idnt_res_*.pdf"
    doc: |
      Split by dataset clustered cells UMAP.
      PDF format

  cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_idnt_res_*.png"
    doc: |
      Grouped by cluster split by dataset cells composition plot. Downsampled.
      PNG format

  cmp_gr_clst_spl_idnt_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_idnt_res_*.pdf"
    doc: |
      Grouped by cluster split by dataset cells composition plot. Downsampled.
      PDF format

  cmp_gr_idnt_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_clst_res_*.png"
    doc: |
      Grouped by dataset split by cluster cells composition plot. Downsampled.
      PNG format

  cmp_gr_idnt_spl_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_clst_res_*.pdf"
    doc: |
      Grouped by dataset split by cluster cells composition plot. Downsampled.
      PDF format

  umap_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_res_*.png"
    doc: |
      Split by grouping condition clustered cells UMAP.
      PNG format

  umap_spl_cnd_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_res_*.pdf"
    doc: |
      Split by grouping condition clustered cells UMAP.
      PDF format

  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_cnd_res_*.png"
    doc: |
      Grouped by cluster split by condition cells composition plot. Downsampled.
      PNG format

  cmp_gr_clst_spl_cnd_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_cnd_res_*.pdf"
    doc: |
      Grouped by cluster split by condition cells composition plot. Downsampled.
      PDF format

  cmp_gr_cnd_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_clst_res_*.png"
    doc: |
      Grouped by condition split by cluster cells composition plot. Downsampled.
      PNG format

  cmp_gr_cnd_spl_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_clst_res_*.pdf"
    doc: |
      Grouped by condition split by cluster cells composition plot. Downsampled.
      PDF format

  umap_spl_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_ph_res_*.png"
    doc: |
      Split by cell cycle phase clustered cells UMAP.
      PNG format

  umap_spl_ph_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_ph_res_*.pdf"
    doc: |
      Split by cell cycle phase clustered cells UMAP.
      PDF format

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.png"
    doc: |
      Grouped by cell cycle phase split by dataset cells composition plot. Downsampled.
      PNG format

  cmp_gr_ph_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.pdf"
    doc: |
      Grouped by cell cycle phase split by dataset cells composition plot. Downsampled.
      PDF format

  cmp_gr_ph_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_ph_spl_clst_res_*.png"
    doc: |
      Grouped by cell cycle phase split by cluster cells composition plot. Downsampled.
      PNG format

  cmp_gr_ph_spl_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_ph_spl_clst_res_*.pdf"
    doc: |
      Grouped by cell cycle phase split by cluster cells composition plot. Downsampled.
      PDF format

  xpr_avg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_avg_res_*.png"
    doc: |
      Log normalized scaled average gene expression per cluster.
      PNG format

  xpr_avg_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_avg_res_*.pdf"
    doc: |
      Log normalized scaled average gene expression per cluster.
      PDF format

  xpr_per_cell_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_[!sgnl_]*.png"
    doc: |
      Log normalized gene expression on cells UMAP.
      PNG format

  xpr_per_cell_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_[!sgnl_]*.pdf"
    doc: |
      Log normalized gene expression on cells UMAP.
      PDF format

  xpr_per_cell_sgnl_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_*.png"
    doc: |
      Log normalized gene expression density on cells UMAP.
      PNG format

  xpr_per_cell_sgnl_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_*.pdf"
    doc: |
      Log normalized gene expression density on cells UMAP.
      PDF format

  xpr_dnst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_res_*.png"
    doc: |
      Log normalized gene expression density per cluster.
      PNG format

  xpr_dnst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_res_*.pdf"
    doc: |
      Log normalized gene expression density per cluster.
      PDF format

  xpr_htmp_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_htmp_res_*.png"
    doc: |
      Normalized gene expression heatmap grouped by cluster.
      PNG format

  xpr_htmp_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_htmp_res_*.pdf"
    doc: |
      Normalized gene expression heatmap grouped by cluster.
      PDF format

  gene_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_gene_markers.tsv"
    doc: |
      Differentially expressed genes between each pair of clusters for all resolutions.
      TSV format

  ucsc_cb_config_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser"
    doc: |
      Directory with UCSC Cellbrowser configuration data.

  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser/html_data"
    doc: |
      Directory with UCSC Cellbrowser html data.

  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: "*_cellbrowser/html_data/index.html"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser html data.

  seurat_data_rds:
    type: File
    outputBinding:
      glob: "*_data.rds"
    doc: |
      Reduced Seurat data in RDS format

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Reduced Seurat data in h5seurat format

  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_data.h5ad"
    doc: |
      Reduced Seurat data in h5ad format

  seurat_data_scope:
    type: File?
    outputBinding:
      glob: "*_data.loom"
    doc: |
      Reduced Seurat data in SCope compatible loom format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_screpertoire_multisample.R"]

stdout: screpetoire_stdout.log
stderr: screpertoire_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell TCR-seq and BCR-seq Clonal Analysis"
s:name: "Single-cell TCR-seq and BCR-seq Clonal Analysis"
s:alternateName: "Clonal analysis of single-cell TCR-seq and BCR-seq data"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-rna-cluster.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Single-cell TCR-seq and BCR-seq Clonal Analysis

  A toolkit for single-cell immune receptor profiling. scRepertoire was built to 
  process data derived from the 10x Genomics Chromium Immune Profiling for both T-cell 
  receptor (TCR) and immunoglobulin (Ig) enrichment workflows and subsequently interacts 
  with the popular Seurat and SingleCellExperiment R packages. It also allows for general 
  analysis of single-cell clonotype information without the use of expression 
  information. The package functions as a wrapper for Startrac and powerTCR R packages.


s:about: |
  usage: run_screpertoire_multisample.R
        [-h] --seuratobject SEURATOBJECT [--contigfile [CONTIGFILE [CONTIGFILE ...]]]
        [--samples [SAMPLES [SAMPLES ...]]] [--cluster [CLUSTER [CLUSTER ...]]]
        [--immunetype {TCR,BCR,both}]
        [--tcells {alpha-beta,gamma-delta,both}] 
        [--removeMulti] [--filterMulti] [--removeMulti] [--removeNA] 
        [--clone {gene,nt,aa,strict}]
       
  Single-cell TCR-seq and BCR-seq Clonal Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --seuratobject SEURATOBJECT
                          Path to the seurat object that contains corresponding gene 
                          expression of the TCR/BCR sequencing in --contigfile argument. 
                          Note: Cluster analysis must be done either through Seurat or 
                          through an external program such that the seurat object 
                          contains atleast one set of clusters assigned to the cells.
    --contigfile [CONTIGFILE [CONTIGFILE ...]]        
    					  Path to the cellranger generated 
    					  filtered_contig_annotation.csv file
    --samples [SAMPLES [SAMPLES ...]]
                          name of the samples in seuratobject that indicates 
                          the samples to analyze. Can take more than 1 sample names
    --cluster [CLUSTER [CLUSTER ...]]
                          name of the metadata header in seuratobject that indicates
                          the clusters of interest. Default: 'seurat_clusters'. Can 
                          take more than 1 set of clusters
    --immunetype {TCR,BCR,both}
                          Is it TCR or BCR data or contains both type of immune receptors? 
                          Options are: 'TCR', 'BCR',or 'both'.
    --tcells {alpha-beta,gamma-delta,both}
                          Enrichment of which type of T cells if --immunetype is 'TCR' or
                          'both'. Options are 'alpha-beta','gamma-delta', or 'both'. 
    --removeMulti         if set to 'True', this is a stringent filter to remove 
                          any cell barcode with more than 2 immune receptor chains. 
                          Default: False includes and incorporates  cells with > 2 chains.
    --filterMulti         if set to 'True', it isolates the top 2 expressed chains in 
    					  cell barcodes with multiple chains. Default: False includes 
    					  and incorporates cells with > 2 chains. 
    --removeNA            if set to 'True', this is a stringent filter to remove any 
    					  cell barcode with an NA value in at least one of the chains. 
    					  Default: False includes and incorporates cells with 1 NA value.
    --clone {gene,nt,aa,strict}             
    					  Definition of clone: Options are 'gene', 'nt','aa','strict'. 
    					  'gene' - use the VDJC genes comprising the TCR/Ig
    					  'nt' - use the nucleotide sequence of the CDR3 region
    					  'aa' - use the amino acid sequence of the CDR3 region
    					  'strict' - DEFAULT: use the VDJC genes comprising the TCR/Ig + 
    					  the nucleotide sequence of the CDR3 region. This is the proper 
    					  definition of clonotype. Note that the clonotype is called 
    					  using the combination of genes or nt/aa CDR3 sequences for both 
    					  loci.'gene' is most sensitive, 'nt' or 'aa' moderately so, 
    					  and 'strict' the most specific. Default: strict