import("scRepertoire", attach=FALSE) 
import("Seurat", attach=FALSE)
import("RColorBrewer", attach=FALSE)
import("ggraph", attach=FALSE)
import("circlize", attach=FALSE)

export(
  "create_combined_sample_object",
  "visualize_vdjc_genes",
  "create_vdjc_seurat_object",
  "check_seurat_object",
  "visualize_bcr_tcr",
  "visualize_clonalhomeostasis",
  "visualize_clonaloverlap",
  "visualize_clonalnetwork",
  "visualize_chord_diagrams",
  "visualize_clonaldiversity"
)

suppressMessages(ucsc <- modules::use("/Users/tha8tf/Documents/GitHub/immunereceptoire-clones/ucsc.R"))
#suppressMessages(graphics <- modules::use("graphics.R"))

D40_COLORS <- c("#FB1C0D", "#0DE400", "#0D00FF", "#E8B4BD", "#FD00EA", "#0DD1FE", "#FF9B0D", "#0D601C", "#C50D69", "#CACA16", "#722A91", "#00DEBF", "#863B00", "#5D7C91", "#FD84D8", "#C100FB", "#8499FC", "#FD6658", "#83D87A", "#968549", "#DEB6FB", "#832E60", "#A8CAB0", "#FE8F95", "#FE1CBB", "#DF7CF8", "#FF0078", "#F9B781", "#4D493B", "#1C5198", "#7C32CE", "#EFBC16", "#7CD2DE", "#B30DA7", "#9FC0F6", "#7A940D", "#9B0000", "#946D9B", "#C8C2D9", "#94605A")


# Useful links
# 1. https://github.com/ncborcherding/scRepertoire
# 2. https://ncborcherding.github.io/vignettes/vignette.html
# 3. https://ncborcherding.github.io/vignettes/beyond_screpertoire.html

set.seed(45)

create_combined_sample_object <- function(contig_list,samples ="S1", immunetype,removeNA,removeMulti,filterMulti,tcells) {
  
  ## if there are two types of t cells, how do I want to deal with this list? 
  
  if (immunetype == "TCR"){
    if (tcells == "alpha-beta") {
      combined_sample_object <- scRepertoire::combineTCR(contig_list,samples = samples,cells = "T-AB", removeNA = removeNA, removeMulti = removeMulti, filterMulti = filterMulti)
    } else {
      combined_sample_object <- scRepertoire::combineTCR(contig_list,samples = samples,cells = "T-GD", removeNA = removeNA, removeMulti = removeMulti, filterMulti = filterMulti)
    }
    
  } else{
    combined_sample_object <- scRepertoire::combineBCR(contig_list,samples = samples, removeNA = removeNA, removeMulti = removeMulti)
  } 
  
  # else (args$immunetype == "both"){
  #   print("two types of immune data")
  # }
  
  return (combined_sample_object)
}


visualize_vdjc_genes <- function(combined_sample_object,immunetype){
  
  gene_vector = c("V","D","J","C")
  
  if (immunetype == "TCR"){
    chain_vec = c("TRA","TRB")
  } else{
    chain_vec = c("IGH","IGL")
  }
  
  for (g in gene_vector){
    for (c in chain_vec){
      
      tryCatch(
        {grDevices::pdf(base::paste0(c,"_",g,"_gene_visualize_frequencies.pdf"),width = 11, height=11)
          p <- scRepertoire::vizGenes(combined_sample_object, gene = g, chain = c, plot = "bar", order = "variance", scale = FALSE)
          base::print(p)
          grDevices::dev.off()
          
          grDevices::pdf(base::paste0(c,"_",g,"_gene_visualize_frequencies_scaled.pdf"),width = 11, height=11)
          p <- scRepertoire::vizGenes(combined_sample_object, gene = g, chain = c, plot = "bar", order = "variance", scale = TRUE)
          base::print(p)
          grDevices::dev.off()
        }, error = function(e){
          base::print(paste("No", g,"genes detected", "in", c, "chain",e))
        }
      )
    }
  }
}


create_vdjc_seurat_object <- function(seuratobject,cluster,combined_sample_object){
  
  seurat_object = check_seurat_object(seuratobject,cluster)
  
  ## add a check here to see if the seurat object already has S1 in its prefix or suffix so it doesnt do this unnecessarily.
  
  ## then add a check if the barcodes in general match between seurat object and the combined object
  
  # seurat_object = Seurat::RenameCells(seurat_object,add.cell.id = "S1")
  
  seurat_object = Seurat::RenameCells(object = seurat_object,new.names = base::paste0(seurat_object@meta.data[["orig.ident"]],"_", Seurat::Cells(x = seurat_object)))
  
  vdjc_seurat_object = scRepertoire::combineExpression(combined_sample_object,seurat_object,addLabel = T)
  
  return(vdjc_seurat_object)
  
}

check_seurat_object <- function(seuratobject,cluster){
  
  so <- readRDS(seuratobject)
  
  for (i in 1:base::length(cluster)){
    metadata_header <- cluster[i]
    
    tryCatch(
      { check_if_factor = TRUE 
      if (base::class(so@meta.data[[metadata_header]]) != "factor"){
        check_if_factor = FALSE
        base::print(paste("Performing analyses for" , metadata_header, "that has the following unique identities:", unique(base::sort(so@meta.data[[metadata_header]]))))
      }
      if (check_if_factor == TRUE){
        base::print(paste("Performing analyses for" , metadata_header, "that has the following unique identities:", base::levels(so@meta.data[[metadata_header]])))
      }
      }, error = function(e){
        base::print(paste("Error with", metadata_header,"in seurat object with error",e)) ## should include a code to ignore if the mentioned metadata header is missing and run the analyses on the remaining headers. 
      }
    )
  }
  return(so) 
}


visualize_bcr_tcr <- function(clone, cluster, vdjc_seurat_object,combined_sample_object,immunetype){
  
  if (immunetype == "TCR"){
    chain_vec = c("TRA","TRB")
  } else{
    chain_vec = c("IGH","IGL")
  }
  
  for (clonecall_type in clone){
    for (chain in chain_vec){
      for (c in cluster){
        
        tryCatch({
          grDevices::pdf(base::paste0(c,"_",clonecall_type,"_Number_of_Unique_Clonotypes.pdf"),width = 11, height=11)
          p <- scRepertoire::quantContig(vdjc_seurat_object, cloneCall=clonecall_type, scale = F, split.by = c)  + ggplot2::labs(fill=c) + ggplot2::xlab(c) + ggplot2::scale_fill_manual(values=D40_COLORS) #+ ggplot2::scale_color_manual(values = metadata_color_list[c]) 
          base::print(p)
          grDevices::dev.off()
          
          quantContig_output <- scRepertoire::quantContig(vdjc_seurat_object, cloneCall=clonecall_type, scale = F, split.by = c, exportTable = T)
          base::colnames(quantContig_output)[2] = c
          utils::write.table(quantContig_output,base::paste0(c,"_",clonecall_type,"_Number_of_Unique_Clonotypes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
          
          grDevices::pdf(base::paste0(c,"_",clonecall_type,"_Percentage_of_Unique_Clonotypes.pdf"),width = 11, height=11)
          p <- scRepertoire::quantContig(vdjc_seurat_object, cloneCall=clonecall_type, scale = T, split.by = c)  + ggplot2::labs(fill=c) + ggplot2::xlab(c) + ggplot2::scale_fill_manual(values=D40_COLORS) #+ ggplot2::scale_fill_manual(values = metadata_color_list[c])
          base::print(p)
          grDevices::dev.off()
          
          quantContig_output <- scRepertoire::quantContig(vdjc_seurat_object, cloneCall=clonecall_type, scale = T, split.by = c, exportTable = T)
          base::colnames(quantContig_output)[2] = c
          utils::write.table(quantContig_output,base::paste0(c,"_",clonecall_type,"_Percentage_of_Unique_Clonotypes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
          
          grDevices::pdf(base::paste0(c,"_",clonecall_type,"_",chain,"_Number_of_Unique_Clonotypes.pdf"),width = 11, height=11)
          p <- scRepertoire::quantContig(vdjc_seurat_object, cloneCall=clonecall_type, scale = F, split.by = c,chain = chain)  + ggplot2::labs(fill=c) + ggplot2::xlab(c) + ggplot2::scale_fill_manual(values=D40_COLORS) #+ ggplot2::scale_fill_manual(values = metadata_color_list[c])
          base::print(p)
          grDevices::dev.off()
          
          quantContig_output <- scRepertoire::quantContig(vdjc_seurat_object, cloneCall=clonecall_type, scale = F, split.by = c, exportTable = T,chain = chain)
          base::colnames(quantContig_output)[2] = c
          utils::write.table(quantContig_output,base::paste0(c,"_",clonecall_type,"_",chain,"_Number_of_Unique_Clonotypes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
          
          grDevices::pdf(base::paste0(c,"_",clonecall_type,"_",chain, "_Percentage_of_Unique_Clonotypes.pdf"),width = 11, height=11)
          p <- scRepertoire::quantContig(vdjc_seurat_object, cloneCall=clonecall_type, scale = T, split.by = c,chain = chain)  + ggplot2::labs(fill=c) + ggplot2::xlab(c) + ggplot2::scale_fill_manual(values=D40_COLORS) #+ ggplot2::scale_fill_manual(values = metadata_color_list[c])
          base::print(p)
          grDevices::dev.off()
          
          quantContig_output <- scRepertoire::quantContig(vdjc_seurat_object, cloneCall=clonecall_type, scale = T, split.by = c, exportTable = T,chain = chain)
          base::colnames(quantContig_output)[2] = c
          utils::write.table(quantContig_output,base::paste0(c,"_",clonecall_type,"_",chain,"_Percentage_of_Unique_Clonotypes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
          
          ### run the below only if clonecall includes "nt" or "aa" ### 
          # grDevices::pdf(base::paste0(clonecall_type,"_",chain, "_Length_of_Unique_Clonotypes.pdf"),width = 11, height=11)
          # p = lengthContig(combined_sample_object, cloneCall=clonecall_type, chain = chain) 
          # base::print(p)
          # grDevices::dev.off()
          # 
          # grDevices::pdf(base::paste0(clonecall_type, "_Length_of_Unique_Clonotypes.pdf"),width = 11, height=11)
          # p = lengthContig(combined_sample_object, cloneCall=clonecall_type) 
          # base::print(p)
          # grDevices::dev.off()
          
          
        }, error = function(e){
          base::print(paste("An error in determination of Number/Percent of Clonotypes for", clonecall_type, "and",chain, "using", c,e))
        })
        
        
      }
    }
  }
  
}

visualize_clonalhomeostasis <- function(clone,cluster,vdjc_seurat_object,combined_sample_object){
  
  for (clonecall_type in clone){
    for (c in cluster){
      tryCatch({
        
        grDevices::pdf(base::paste0(clonecall_type,"_",c, "_Clonal_Homeostasis.pdf"),width = 11, height=11)
        p <- scRepertoire::clonalHomeostasis(vdjc_seurat_object, cloneCall=clonecall_type, split.by = c)  + ggplot2::labs(fill=c) + ggplot2::xlab(c) #+ ggplot2::scale_color_manual(values = metadata_color_list[c])
        base::print(p)
        grDevices::dev.off()
        
        clonal_homeostasis_output <- scRepertoire::clonalHomeostasis(vdjc_seurat_object, cloneCall=clonecall_type, split.by = c, exportTable = T)
        utils::write.table(clonal_homeostasis_output,base::paste0(clonecall_type,"_",c,"_Clonal_Homeostasis.txt"), sep = "\t", col.names = T, row.names = T, quote = F)
        
        grDevices::pdf(base::paste0(clonecall_type,"_",c, "_Clonal_Homeostasis_UMAP.pdf"),width = 11, height=11)
        temp_seurat_object <- scRepertoire::combineExpression(combined_sample_object, vdjc_seurat_object, cloneCall=clonecall_type)
        p = Seurat::DimPlot(temp_seurat_object, group.by = "cloneType") + ggplot2::theme(plot.title = ggplot2::element_blank()) 
        base::print(p)
        grDevices::dev.off()
        
        print("Exporting Clonal Homeostasis to UCSC Cellbrowser")
        ucsc$export_cellbrowser(
          seurat_data=temp_seurat_object,
          assay="RNA",
          slot="counts",
          short_label="RNA",
          is_nested=TRUE,
          meta_fields = "cloneType", meta_fields_names = "ClonalHomeostasis",
          rootname=paste("test_clonal", "_cellbrowser/rna", sep=""),
        )
        
        
        rm(temp_seurat_object)
        
      }, error = function(e){
        base::print(paste("Error in Clonal Homeostasis for clone call type",clonecall_type, "using",c,e))
      }
      )
    }
  }
  
}


visualize_clonaloverlap <- function(clone, cluster, vdjc_seurat_object){
  
  method_vec = c("overlap", "morisita", "jaccard", "raw")
  for (clonecall_type in clone){
    for (c in cluster){
      for (m in method_vec){
        tryCatch({
          
          grDevices::pdf(base::paste0(clonecall_type,"_",c,"_",m,"_Clonal_Overlap.pdf"),width = 11, height=11)
          p <- scRepertoire::clonalOverlap(vdjc_seurat_object, cloneCall=clonecall_type, split.by = c,method = m)  + ggplot2::labs(fill=c) + ggplot2::xlab(c) #+ ggplot2::scale_fill_manual(values = metadata_color_list[c])
          base::print(p)
          grDevices::dev.off()
          
          clonal_overlap_output <- scRepertoire::clonalOverlap(vdjc_seurat_object, cloneCall=clonecall_type, split.by = c, exportTable = T,method = m)
          base::colnames(clonal_overlap_output)[base::ncol(clonal_overlap_output)] = c
          utils::write.table(clonal_overlap_output,base::paste0(clonecall_type,"_",c,"_",m,"_Clonal_Overlap.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
          
          
        }, error = function(e){
          base::print(paste("Error in Clonal Overlap for clone call type",clonecall_type, "using",c, "and method",m,e))
        }
        )
      }
    }
  }
  
}


visualize_clonalnetwork <- function(clone, cluster, vdjc_seurat_object, combined_sample_object){
  
  for (clonecall_type in clone){
    for (c in cluster){
      tryCatch({
        
        temp_seurat_object <- scRepertoire::combineExpression(combined_sample_object, vdjc_seurat_object, cloneCall=clonecall_type)
        grDevices::pdf(base::paste0(clonecall_type,"_",c,"_Clonal_Network.pdf"),width = 11, height=11)
        p <- scRepertoire::clonalNetwork(temp_seurat_object, reduction = "rnaumap",identity = c, cloneCall = clonecall_type,filter.clones = NULL,filter.identity = NULL) + ggplot2::scale_color_manual(values=D40_COLORS) + ggplot2::scale_fill_manual(values=D40_COLORS) 
        base::print(p)
        grDevices::dev.off()
        
        clonal_network_output <- scRepertoire::clonalNetwork(temp_seurat_object, reduction = "rnaumap",identity = c, cloneCall = clonecall_type,exportTable = T)
        utils::write.table(clonal_network_output,base::paste0(clonecall_type,"_",c,"_Clonal_Overlap.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
        
        rm(temp_seurat_object)  
      }, error = function(e){
        base::print(paste("Error in Clonal Network for clone call type",clonecall_type, "using",c,e))
      }
      )
    }
  }
  
}


visualize_chord_diagrams <- function(clone, cluster,vdjc_seurat_object, combined_sample_object){
  
  for (clonecall_type in clone){
    for (c in cluster){
      tryCatch({
        
        temp_seurat_object <- scRepertoire::combineExpression(combined_sample_object, vdjc_seurat_object, cloneCall=clonecall_type)
        grDevices::pdf(base::paste0(clonecall_type,"_",c,"_Clonal_ChordDiagram.pdf"),width = 11, height=11)
        circles <- scRepertoire::getCirclize(temp_seurat_object,group.by = c)
        grid_cols = D40_COLORS[1:base::length(unique(temp_seurat_object@meta.data[[c]]))]
        names(grid_cols) = base::sort(unique(temp_seurat_object@meta.data[[c]]))
        p <- circlize::chordDiagram(circles,self.link = 1,grid.col=grid_cols)
        base::print(p)
        grDevices::dev.off()
        
      }, error = function(e){
        base::print(paste("Error in Chord Diagram for clone call type",clonecall_type, "using",c,e))
      }
      )
    }
  }
  
}

visualize_clonaldiversity <- function(clone, cluster,vdjc_seurat_object){
  
  for (clonecall_type in clone){
    for (c in cluster){
      
      tryCatch({
        
        grDevices::pdf(base::paste0(clonecall_type,"_",c,"_Clonal_Diversity.pdf"),width = 11, height=11)
        p <- scRepertoire::clonalDiversity(vdjc_seurat_object, cloneCall=clonecall_type, split.by = c)  + ggplot2::labs(fill=c) + ggplot2::xlab(c) #+ ggplot2::scale_fill_manual(values = metadata_color_list[c])
        base::print(p)
        grDevices::dev.off()
        
        clonal_diversity_output <- scRepertoire::clonalDiversity(vdjc_seurat_object, cloneCall=clonecall_type, split.by = c,exportTable = T)
        utils::write.table(clonal_diversity_output,base::paste0(clonecall_type,"_",c,"_Clonal_Diversity.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
        
      }, error = function(e){
        base::print(paste("Error in Clonal Diversity for clone call type",clonecall_type, "using",c,e))
      }
      )
    }
  }
  
}
