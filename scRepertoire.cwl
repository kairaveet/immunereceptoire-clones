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

  number_of_unique_clonotypes:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Number_of_Unique_Clonotypes.pdf"
    doc: |
      bar plot of absolute number of clonotypes per cluster/group.

  percentage_of_unique_clonotypes:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Percentage_of_Unique_Clonotypes.pdf"
    doc: |
      bar plot of percentage of clonotypes per cluster/group.

  clonal_chord_diagram:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Clonal_ChordDiagram.pdf"
    doc: |
      chord diagram of number of clones shared between groups.    

  diversity_plot:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Diversity.pdf"
    doc: |
      box plot of five measures of clonal diversity per group.

  clonal_homeostasis_umap:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Clonal_Homeostasis_UMAP.pdf"
    doc: |
      clonal homeostasis category for each cell on UMAP.

  clonal_homeostasis:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Clonal_Homeostasis.pdf"
    doc: |
      stacked bar plot of clonal homeostasis of cells in each group.     

  clonal_network_umap:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Clonal_Network.pdf"
    doc: |
      network drawn on umap showing sharing of clones
      
  jaccard_clonal_overlap:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_jaccard_Clonal_Overlap.pdf"
    doc: |
      jaccard metric of clonal overlap between two groups

  morisita_clonal_overlap:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_morisita_Clonal_Overlap.pdf"
    doc: |
      morisita metric of clonal overlap between two groups

  scaled_clonal_overlap:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_overlap_Clonal_Overlap.pdf"
    doc: |
      scaled clonal overlap between two groups 

  absolute_clonal_overlap:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_raw_Clonal_Overlap.pdf"
    doc: |
       absolute number of shared clones between two groups     

  vdjc_scaled_frequencies:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_visualize_frequencies_scaled.pdf"
    doc: |
       bar plot showing percentage of cells with each vdjc gene in each group  

  vdjc_absolute_freqeuncies:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_visualize_frequencies.pdf"
    doc: |
      bar plot showing number of cells with each vdjc gene in each group 

  number_of_unique_clonotypes_txt:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Number_of_Unique_Clonotypes.txt"
    doc: |
      absolute number of clonotypes per cluster in txt.

  percentage_of_unique_clonotypes_txt:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Percentage_of_Unique_Clonotypes.txt"
    doc: |
      percentage of clonotypes per cluster/group in txt.

  diversity_txt:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Diversity.txt"
    doc: |
      five measures of clonal diversity per group in txt.

  clonal_homeostasis_txt:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_Clonal_Homeostasis_UMAP.txt"
    doc: |
      clonal homeostasis of cells in each group in txt.  
      
  jaccard_clonal_overlap_txt:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_jaccard_Clonal_Overlap.txt"
    doc: |
       jaccard clonal overlap score between two groups in txt.  

  morisita_clonal_overlap_txt:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_morisita_Clonal_Overlap.txt"
    doc: |
       morisita clonal overlap score between two groups in txt.   

  scaled_clonal_overlap_txt:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_overlap_Clonal_Overlap.txt"
    doc: |
       scaled clonal overlap score between two groups in txt.   

  absolute_clonal_overlap_txt:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_raw_Clonal_Overlap.txt"
    doc: |
       absolute/raw clonal overlap score between two groups in txt.  

  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: "*_clonal_ellbrowser/html_data"
    doc: |
      Directory with UCSC Cellbrowser html data.

  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: "*_clonal_cellbrowser/html_data/index.html"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser html data.

  vdjc_seurat_object:
    type: File
    outputBinding:
      glob: "vdjc_seurat_object.rds"
    doc: |
      Reduced Seurat data in RDS format


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