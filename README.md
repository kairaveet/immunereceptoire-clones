# immunereceptoire-clones
scRepertoire analysis for sciDAP

Command to run the above files: 

Rscript /Users/tha8tf/Documents/GitHub/immunereceptoire-clones/run_screpertoire_multisample.R --contigfile "/Users/tha8tf/Downloads/filtered_contig_annotations_e3.csv" "/Users/tha8tf/Downloads/filtered_contig_annotations_e2.csv" "/Users/tha8tf/Downloads/filtered_contig_annotations_e0.csv" "/Users/tha8tf/Downloads/filtered_contig_annotations_e1.csv" --seuratobject "/Users/tha8tf/Downloads/sc_data_e0toe3.rds" --immunetype "TCR" --tcells "alpha-beta" --clone "strict" --cluster "seurat_clusters" "orig.ident" --samples "1" "2" "3" "4"
