import("dplyr", attach=FALSE)
import("Seurat", attach=FALSE)
import("data.table", attach=FALSE)
import("reticulate", attach=FALSE)
import("tidyselect", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)

export(
    "export_cellbrowser"
)


get_matrix <- function(object, slot){
    if (slot == "counts") {
        counts <- SeuratObject::GetAssayData(object=object, slot="counts")
    } else if (slot == "scale.data") {
        counts <- SeuratObject::GetAssayData(object=object, slot="scale.data")
    } else if (slot == "data") {
        counts <- SeuratObject::GetAssayData(object=object)
    } else {
        base::print(base::paste("Slot", slot, "not found. Quit."))
        base::quit(save="no", status=1, runLast=FALSE)
    }
}

write_matrix <- function(inMat, outFname, sliceSize=1000){
    fnames <- c()
    mat <- inMat
    geneCount <- base::nrow(mat)
    startIdx <- 1
    while (startIdx < geneCount){
        endIdx <- min(startIdx + sliceSize - 1, geneCount)
        matSlice <- mat[startIdx:endIdx,]
        denseSlice <- base::as.matrix(matSlice)
        dt <- data.table::data.table(denseSlice)
        dt <- base::cbind(gene=base::rownames(matSlice), dt)
        writeHeader <- FALSE
        if (startIdx == 1) { 
            writeHeader <- TRUE
        }
        sliceFname <- base::paste0("temp", startIdx, ".txt")
        data.table::fwrite(dt, sep="\t", file=sliceFname, quote=FALSE, col.names=writeHeader)
        fnames <- base::append(fnames, sliceFname)
        startIdx <- startIdx + sliceSize
    }
    base::system(base::paste("cat", base::paste(fnames, collapse=" "), "| gzip >", outFname, sep=" "))
    base::unlink(fnames)
}

cb_build <- function(object, slot, short_label, rootname, label_field, show_labels, features, markers, meta_fields, meta_fields_names, is_nested, dot_radius, dot_alpha, color_data) {
    if (!base::dir.exists(rootname)) {
        base::dir.create(rootname, recursive=TRUE)
    }
    idents <- SeuratObject::Idents(object)
    meta <- object@meta.data
    cell_order <- base::colnames(object)
    counts <- get_matrix(object, slot)
    genes <- base::rownames(x=object)
    dr <- object@reductions

    gzPath <- base::file.path(rootname, "exprMatrix.tsv.gz")
    if ((((base::ncol(counts)/1000)*(base::nrow(counts)/1000))>2000) && methods::is(counts, "sparseMatrix")){
        write_matrix(counts, gzPath);
    } else {
        mat <- base::as.matrix(counts)
        df <- base::as.data.frame(mat, check.names=FALSE)
        df <- base::data.frame(gene=genes, df, check.names=FALSE)
        z <- base::gzfile(gzPath, "w")
        utils::write.table(x=df, sep="\t", file=z, quote=FALSE, row.names=FALSE)
        base::close(con=z)
    }
    embeddings = names(dr)
    embeddings_conf <- c()
    for (embedding in embeddings) {
        emb <- dr[[embedding]]
        df <- emb@cell.embeddings
        if (base::ncol(df) > 2){
            df <- df[, 1:2]
        }
        base::colnames(df) <- c("x", "y")
        df <- base::data.frame(cellId=base::rownames(df), df, check.names=FALSE)
        utils::write.table(
            df[cell_order,],
            sep="\t",
            file=base::file.path(rootname, base::sprintf("%s.coords.tsv", embedding)),
            quote=FALSE,
            row.names=FALSE
        )
        embeddings_conf <- c(
            embeddings_conf,
            base::sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}', embedding)
        )
    }
    df <- base::data.frame(row.names=cell_order, check.names=FALSE)
    enum_fields <- c()
    for (i in 1:length(meta_fields)){
        field <- meta_fields[i]
        name <- meta_fields_names[i]
        if (field %in% base::colnames(meta)){
            df[[name]] <- meta[[field]]
            if (!is.numeric(df[[name]])) {
                enum_fields <- c(enum_fields, name)
            }
        }
    }
    df <- base::data.frame(Cell=base::rownames(df), df, check.names=FALSE)
    utils::write.table(
        base::as.matrix(df[cell_order,]),
        sep="\t",
        file=base::file.path(rootname, "meta.tsv"),
        quote=FALSE,
        row.names=FALSE
    )

    local_config <- '
name="%s"
shortLabel="%s"
geneLabel="Feature"
exprMatrix="exprMatrix.tsv.gz"
meta="meta.tsv"
radius=%i
alpha=%f
geneIdType="auto"
clusterField="%s"
labelField="%s"
enumFields=%s
showLabels=%s
coords=%s'

    local_config <- base::sprintf(
        local_config,
        short_label,
        short_label,
        dot_radius,
        dot_alpha,
        label_field,
        label_field,
        base::paste0("[", base::paste(base::paste0('"', enum_fields, '"'), collapse=", "), "]"),
        base::ifelse(show_labels, "True", "False"),
        base::paste0("[", base::paste(embeddings_conf, collapse = ",\n"), "]")
    )

    if (!is.null(features)){
        utils::write.table(
            features,
            sep="\t",
            file=base::file.path(rootname, "quickGenes.tsv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE
        )
        local_config <- base::paste(
            local_config,
            'quickGenesFile="quickGenes.tsv"',
            sep="\n"
        )
    }

    if (!is.null(markers)){
        utils::write.table(
            markers,
            sep="\t",
            file=base::file.path(rootname, "markers.tsv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=TRUE
        )
        local_config <- base::paste(
            local_config,
            'markers=[{"file":"markers.tsv", "shortLabel":"Cluster-specific markers"}]',
            sep="\n"
        )
    }

    if (!is.null(color_data)){
        utils::write.table(
            color_data,
            sep=",",
            file=base::file.path(rootname, "colors.csv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE
        )
        local_config <- base::paste(
            local_config,
            'colors="colors.csv"',
            sep="\n"
        )
    }

    html_data_dir <- base::file.path(rootname, "html_data")
    if (is_nested){
        local_config <- base::paste(local_config, 'dataRoot="../"', sep="\n")
        desc <- 'title="%s"\nabstract=""\nmethods=""\nbiorxiv_url=""\ncustom={}'
        desc <- base::sprintf(desc, short_label)
        base::cat(
            desc,
            file=base::file.path(rootname, "desc.conf")
        )
        base::cat(
            'shortLabel="Multiple datasets"',
            file=base::file.path(rootname, "../cellbrowser.conf")
        )
        html_data_dir <- base::file.path(rootname, "../html_data")
    }

    base::cat(
        local_config,
        file=base::file.path(rootname, "cellbrowser.conf")
    )
    #cb <- reticulate::import(module = "cellbrowser")
    cb <- reticulate::import_from_path(module = "cellbrowser",path="/Users/tha8tf/Library/r-miniconda/envs/r-reticulate")
    cb$cellbrowser$build(rootname, html_data_dir)
}

get_color_data <- function(seurat_data, metadata_column, palette_colors){
    categories <- base::levels(seurat_data@meta.data[[metadata_column]])
    if (is.null(categories)) {                                                 # not a factor
        categories <- base::sort(                                              # alphabetically sorted
            base::unique(
                base::as.vector(
                    as.character(seurat_data@meta.data[[metadata_column]])
                )
            )
        )
    }
    color_data <- base::data.frame(
        category=categories,
        color=palette_colors[1:length(categories)],
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    return (color_data)
}


export_cellbrowser <- function(seurat_data, assay, slot, rootname, label_field=NULL, features=NULL, markers=NULL, dot_radius=3, dot_alpha=0.5, palette_colors=NULL, short_label="cellbrowser", is_nested=FALSE, meta_fields=NULL, meta_fields_names=NULL){
    base::tryCatch(
        expr = {
            backup_assay <- SeuratObject::DefaultAssay(seurat_data)
            SeuratObject::DefaultAssay(seurat_data) <- assay                 # safety measure
            SeuratObject::Idents(seurat_data) <- "new.ident"                 # safety measure

            color_data <- base::data.frame(                                 # here we will collect all colors if palette_colors was provided
                category=character(),
                color=character(),
                check.names=FALSE,
                stringsAsFactors=FALSE
            )

            if (is.null(meta_fields) || is.null(meta_fields_names)){
                meta_fields <-       c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "S.Score", "G2M.Score",    "Phase", "rna_doublets", "atac_doublets", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment",       "nucleosome_signal", "frip", "blacklist_fraction")
                meta_fields_names <- c("RNA UMIs",   "Genes",        "Mitochondrial %", "Novelty score",            "S score", "G to M score", "Phase", "RNA doublets", "ATAC doublets", "ATAC UMIs",   "Peaks",         "TSS enrichment score", "Nucleosome signal", "FRiP", "Bl. regions")
            }

            if (length(base::unique(base::as.vector(as.character(seurat_data@meta.data$new.ident)))) > 1){
                meta_fields <- base::append(meta_fields, "new.ident", 0)
                meta_fields_names <- base::append(meta_fields_names, "Dataset", 0)
                if (!is.null(palette_colors)){
                    color_data <- base::rbind(
                        color_data,
                        get_color_data(seurat_data, "new.ident", palette_colors)
                    )
                }
            }

            if (
                all(base::as.vector(as.character(seurat_data@meta.data$new.ident)) != base::as.vector(as.character(seurat_data@meta.data$condition))) &&
                length(base::unique(base::as.vector(as.character(seurat_data@meta.data$condition)))) > 1
            ){
                meta_fields <- base::append(meta_fields, "condition", 1)
                meta_fields_names <- base::append(meta_fields_names, "Condition", 1)
                if (!is.null(palette_colors)){
                    color_data <- base::rbind(
                        color_data,
                        get_color_data(seurat_data, "condition", palette_colors)
                    )
                }
            }

            clustering_fields <- base::grep("_res\\.", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            clustering_fields_names <- base::as.vector(
                base::unlist(                                                                              # when nothing found it returns named list that should be unlisted
                    base::sapply(
                        clustering_fields,
                        function(line) {
                            split_line <- base::unlist(base::strsplit(line, split="_res\\."))
                            base::paste("Clustering (", split_line[1], " ", split_line[2], ")", sep="")
                        }
                    )
                )
            )
            meta_fields <- base::append(meta_fields, clustering_fields)
            meta_fields_names <- base::append(meta_fields_names, clustering_fields_names)
            if (!is.null(palette_colors) && length(clustering_fields) > 0){                                # need to check if clustering_fields is not empty
                for (i in 1:length(clustering_fields)){
                    color_data <- base::rbind(
                        color_data,
                        get_color_data(seurat_data, clustering_fields[i], palette_colors)
                    )
                }
            }

            custom_fields <- base::grep("^custom_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            custom_fields_names <- base::gsub("custom_", "Custom ", custom_fields)
            meta_fields <- base::append(meta_fields, custom_fields)
            meta_fields_names <- base::append(meta_fields_names, custom_fields_names)
            if (!is.null(palette_colors) && length(custom_fields) > 0){                                    # need to check if custom_fields is not empty
                for (i in 1:length(custom_fields)){
                    color_data <- base::rbind(
                        color_data,
                        get_color_data(seurat_data, custom_fields[i], palette_colors)
                    )
                }
            }

            quartile_fields <- base::grep("^quartile_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            quartile_fields_names <- base::gsub("quartile_", "Quartile ", quartile_fields)
            meta_fields <- base::append(meta_fields, quartile_fields)
            meta_fields_names <- base::append(meta_fields_names, quartile_fields_names)
            if (!is.null(palette_colors) && length(quartile_fields) > 0){                                  # need to check if quartile_fields is not empty
                for (i in 1:length(quartile_fields)){
                    color_data <- base::rbind(
                        color_data,
                        get_color_data(seurat_data, quartile_fields[i], palette_colors)
                    )
                }
            }

            show_labels <- TRUE                                                      # by default we try to show labels, but we hide them if cluster_field wasn't provided
            if (is.null(label_field) || !(label_field %in% (meta_fields_names))){    # either not provided or not correct label_field
                show_labels <- FALSE                                                 # hide labels
                markers <- NULL                                                      # we can't use markers if label_field wasn't provided
                for (i in 1:length(meta_fields)){                                    # searching the first field from the meta.data that is not unique for all cells
                    current_field <- meta_fields[i]
                    if (current_field %in% base::colnames(seurat_data@meta.data) && length(base::unique(base::as.vector(as.character(seurat_data@meta.data[, current_field])))) > 1){
                        label_field <- meta_fields_names[i]
                        break
                    }
                }
            }

            if (!is.null(markers) && base::nrow(markers) == 0){                      # in case markers somehow was an empty dataframe
                markers <- NULL
            }

            if (is.null(features) && !is.null(markers)){                             # setting default features from markers only if they are provided with correct label_field
                features <- markers %>%
                            dplyr::group_by(cluster) %>%
                            dplyr::top_n(                                            # can't use dplyr::distinct for grouped data
                                n=tidyselect::all_of(floor(60/length(base::unique(markers$cluster)))),
                                wt=avg_log2FC
                            ) %>%
                            dplyr::pull(feature)                                     # can't use dplyr::distinct for vector
                features <- base::unique(features)                                   # some of the genes may be identical for several clusters so we want to remove them
            }

            color_data <- color_data %>%                                             # to make sure we don't have duplicates
                          dplyr::distinct(category, .keep_all=TRUE)

            cb_build(
                seurat_data,
                slot=slot,
                short_label=short_label,
                rootname=rootname,
                label_field=label_field,
                show_labels=show_labels,
                features=features,
                markers=markers,
                meta_fields=meta_fields,
                meta_fields_names=meta_fields_names,
                is_nested=is_nested,
                dot_radius=dot_radius,
                dot_alpha=dot_alpha,
                color_data=if(base::nrow(color_data) > 0) color_data else NULL
            )
            base::print(base::paste("Exporting UCSC Cellbrowser data to", rootname, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export UCSC Cellbrowser data with error -", e))
        },
        finally = {
            SeuratObject::DefaultAssay(seurat_data) <- backup_assay
        }
    )
}