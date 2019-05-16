#' Plot consensus matrix as a heatmap
#' 
#' The consensus matrix is a NxN 
#' matrix, where N is the number of cells.
#' It represents similarity between the cells based 
#' on the averaging of clustering results from all 
#' combinations of clustering parameters. Similarity 0 
#' (blue) means that the two cells are always assigned to different clusters. 
#' In contrast, similarity 1 (red) means that the two cells are always assigned 
#' to the same cluster. The consensus matrix is clustered by hierarchical 
#' clustering and has a diagonal-block structure. Intuitively, the perfect 
#' clustering is achieved when all diagonal blocks are completely red 
#' and all off-diagonal elements are completely blue.
#' 
#' @name moSC3_plot_consensus
#' @aliases moSC3_plot_consensus, moSC3_plot_consensus,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
moSC3_plot_consensus.SingleCellExperiment <- function(object, k, show_pdata) {
    if (is.null(metadata(object)$moSC3$consensus)) {
        warning(paste0("Please run moSC3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$moSC3$consensus[[as.character(k)]]$hc
    consensus <- metadata(object)$moSC3$consensus[[as.character(k)]]$consensus
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- make_col_ann_for_heatmaps(object, show_pdata)
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            # make same names for the annotation table
            rownames(ann) <- colnames(consensus)
        }
    }
    do.call(pheatmap::pheatmap, c(list(consensus, cluster_rows = hc, cluster_cols = hc, cutree_rows = k, 
        cutree_cols = k, show_rownames = FALSE, show_colnames = FALSE), list(annotation_col = ann)[add_ann_col]))
}

#' @rdname moSC3_plot_consensus
#' @aliases moSC3_plot_consensus
setMethod("moSC3_plot_consensus", signature(object = "SingleCellExperiment"), moSC3_plot_consensus.SingleCellExperiment)

#' Plot consensus matrix as a heatmap
#' 
#' The consensus matrix is a NxN 
#' matrix, where N is the number of cells.
#' It represents similarity between the cells based 
#' on the averaging of clustering results from all 
#' combinations of clustering parameters. Similarity 0 
#' (blue) means that the two cells are always assigned to different clusters. 
#' In contrast, similarity 1 (red) means that the two cells are always assigned 
#' to the same cluster. The consensus matrix is clustered by hierarchical 
#' clustering and has a diagonal-block structure. Intuitively, the perfect 
#' clustering is achieved when all diagonal blocks are completely red 
#' and all off-diagonal elements are completely blue.
#' 
#' @name moSC3_plot_consensus_omics
#' @aliases moSC3_plot_consensus_omics
#' 
#' @param object a list of objects of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
#' @export

moSC3_plot_consensus_omics.list <- function(object, k, show_pdata) {
  if (is.null(metadata(object[[1]])$moSC3$omics_cons)) {
    warning(paste0("Please run moSC3_consensus() first!"))
    return(object)
  }
  #get list with the smallest number of elements
  a = unlist(lapply(object, function(x)  x@int_colData@nrows))
  ind = which.min(a)
  obj = object[[ind]]
  hc <- metadata(obj)$moSC3$consensus[[as.character(k)]]$hc
  consensus <- metadata(obj)$moSC3$omics_cons
  
  add_ann_col <- FALSE
  ann <- NULL
  if (!is.null(show_pdata)) {
    ann <- make_col_ann_for_heatmaps(obj, show_pdata)
    if (!is.null(ann)) {
      add_ann_col <- TRUE
      # make same names for the annotation table
      rownames(ann) <- colnames(consensus)
    }
  }
  do.call(pheatmap::pheatmap, c(list(consensus, cluster_rows = hc, cluster_cols = hc, cutree_rows = k, 
                                     cutree_cols = k, show_rownames = FALSE, show_colnames = FALSE), list(annotation_col = ann)[add_ann_col]))
}

#' @rdname moSC3_plot_consensus_omics
#' @aliases moSC3_plot_consensus_omics
setMethod("moSC3_plot_consensus_omics", signature(object = "list"), moSC3_plot_consensus_omics.list)
#' Plot silhouette indexes of the cells
#' 
#' A silhouette is a quantitative measure of the diagonality of the consensus 
#' matrix. An average silhouette width (shown at the bottom left of the silhouette 
#' plot) varies from 0 to 1, where 1 represents a perfectly block-diagonal 
#' consensus matrix and 0 represents a situation where there is no 
#' block-diagonal structure. The best clustering is achieved when the average 
#' silhouette width is close to 1.
#' 
#' @name moSC3_plot_silhouette
#' @aliases moSC3_plot_silhouette, moSC3_plot_silhouette,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
moSC3_plot_silhouette.SingleCellExperiment <- function(object, k) {
    if (is.null(metadata(object)$moSC3$consensus)) {
        warning(paste0("Please run moSC3_consensus() first!"))
        return(object)
    }
    silh <- metadata(object)$moSC3$consensus[[as.character(k)]]$silhouette
    plot(silh, col = "black")
}

#' @rdname moSC3_plot_silhouette
#' @aliases moSC3_plot_silhouette
setMethod("moSC3_plot_silhouette", signature(object = "SingleCellExperiment"), moSC3_plot_silhouette.SingleCellExperiment)

#' Plot expression matrix used for moSC3 clustering as a heatmap
#' 
#' The expression panel represents the original input expression matrix 
#' (cells in columns and genes in rows) after the gene filter. 
#' Genes are clustered by kmeans with k = 100 (dendrogram on the left) and 
#' the heatmap represents the expression levels of the gene cluster centers 
#' after log2-scaling.
#' 
#' @name moSC3_plot_expression
#' @aliases moSC3_plot_expression, moSC3_plot_expression,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
moSC3_plot_expression.SingleCellExperiment <- function(object, k, show_pdata) {
    if (is.null(metadata(object)$moSC3$consensus)) {
        warning(paste0("Please run moSC3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$moSC3$consensus[[as.character(k)]]$hc
    dataset <- get_processed_dataset(object)
    if (!is.null(metadata(object)$moSC3$svm_train_inds)) {
        dataset <- dataset[, metadata(object)$moSC3$svm_train_inds]
    }
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- make_col_ann_for_heatmaps(object, show_pdata)
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            # make same names for the annotation table
            rownames(ann) <- colnames(dataset)
        }
    }
    
    if(nrow(dataset) > 100) {
        do.call(pheatmap::pheatmap, c(list(dataset, cluster_cols = hc, kmeans_k = 100, cutree_cols = k, 
                                           show_rownames = FALSE, show_colnames = FALSE), list(annotation_col = ann)[add_ann_col]))
    } else {
        do.call(pheatmap::pheatmap, c(list(dataset, cluster_cols = hc, cutree_cols = k, 
                                           show_rownames = FALSE, show_colnames = FALSE), list(annotation_col = ann)[add_ann_col]))
    }
}

#' @rdname moSC3_plot_expression
#' @aliases moSC3_plot_expression
setMethod("moSC3_plot_expression", signature(object = "SingleCellExperiment"), moSC3_plot_expression.SingleCellExperiment)

#' Plot expression of DE genes of the clusters identified by \code{moSC3} as a heatmap
#' 
#' \code{moSC3} plots gene expression profiles of the 50 genes with the lowest p-values. 
#' 
#' @name moSC3_plot_de_genes
#' @aliases moSC3_plot_de_genes, moSC3_plot_de_genes,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param p.val significance threshold used for the DE genes
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
moSC3_plot_de_genes.SingleCellExperiment <- function(object, k, p.val, show_pdata) {
    if (is.null(metadata(object)$moSC3$consensus)) {
        warning(paste0("Please run moSC3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$moSC3$consensus[[as.character(k)]]$hc
    dataset <- get_processed_dataset(object)
    if (!is.null(metadata(object)$moSC3$svm_train_inds)) {
        dataset <- dataset[, metadata(object)$moSC3$svm_train_inds]
    }
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- make_col_ann_for_heatmaps(object, show_pdata)
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            # make same names for the annotation table
            rownames(ann) <- colnames(dataset)
        }
    }
    
    de_genes <- organise_de_genes(object, k, p.val)
    de_genes <- head(de_genes, 50)
    # remove Inf when the p-value is actually 0 (less than the accuracy limit)
    de_genes[de_genes < 1e-17] <- 1e-17
    row_ann <- data.frame(log10_padj = -log10(de_genes))
    rownames(row_ann) <- names(de_genes)
    
    do.call(pheatmap::pheatmap, c(list(dataset[names(de_genes), , drop = FALSE], show_colnames = FALSE, 
        cluster_rows = FALSE, cluster_cols = hc, cutree_cols = k, annotation_row = row_ann, cellheight = 10), 
        list(annotation_col = ann)[add_ann_col]))
}

#' @rdname moSC3_plot_de_genes
#' @aliases moSC3_plot_de_genes
setMethod("moSC3_plot_de_genes", signature(object = "SingleCellExperiment"), moSC3_plot_de_genes.SingleCellExperiment)

#' Plot expression of marker genes identified by \code{moSC3} as a heatmap.
#' 
#' By default the genes with the area under the ROC curve (AUROC) > 0.85 
#' and with the p-value < 0.01 are selected and the top 10 marker 
#' genes of each cluster are visualized in this heatmap.
#' 
#' @name moSC3_plot_markers
#' @aliases moSC3_plot_markers, moSC3_plot_markers,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param auroc area under the ROC curve
#' @param p.val significance threshold used for the DE genes
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
moSC3_plot_markers.SingleCellExperiment <- function(object, k, auroc, p.val, show_pdata) {
    if (is.null(metadata(object)$moSC3$consensus)) {
        warning(paste0("Please run moSC3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$moSC3$consensus[[as.character(k)]]$hc
    dataset <- get_processed_dataset(object)
    if (!is.null(metadata(object)$moSC3$svm_train_inds)) {
        dataset <- dataset[, metadata(object)$moSC3$svm_train_inds]
    }
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- make_col_ann_for_heatmaps(object, show_pdata)
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            # make same names for the annotation table
            rownames(ann) <- colnames(dataset)
        }
    }
    
    # get all marker genes
    markers <- organise_marker_genes(object, k, p.val, auroc)
    
    if(!is.null(markers)) {
        # get top 10 marker genes of each cluster
        markers <- markers_for_heatmap(markers)
        
        row.ann <- data.frame(Cluster = factor(markers[, 1], levels = unique(markers[, 1])))
        rownames(row.ann) <- markers$feature_symbol
        
        do.call(pheatmap::pheatmap, c(list(dataset[markers$feature_symbol, , drop = FALSE], show_colnames = FALSE, 
            cluster_rows = FALSE, cluster_cols = hc, cutree_cols = k, annotation_row = row.ann, annotation_names_row = FALSE, 
            gaps_row = which(diff(markers[, 1]) != 0), cellheight = 10), list(annotation_col = ann)[add_ann_col]))
    } else {
        message("No markers have been found, try to lower significance thresholds!")
    }
}

#' @rdname moSC3_plot_markers
#' @aliases moSC3_plot_markers
setMethod("moSC3_plot_markers", signature(object = "SingleCellExperiment"), moSC3_plot_markers.SingleCellExperiment)

#' Plot stability of the clusters
#' 
#' Stability index shows how stable each cluster is accross the selected 
#' range of ks. The stability index varies between 0 and 1, where 1 means that 
#' the same cluster appears in every solution for different k.
#' 
#' @name moSC3_plot_cluster_stability
#' @aliases moSC3_plot_cluster_stability, moSC3_plot_cluster_stability,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' 
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw labs ylim
moSC3_plot_cluster_stability.SingleCellExperiment <- function(object, k) {
    if (is.null(metadata(object)$moSC3$consensus)) {
        warning(paste0("Please run moSC3_consensus() first!"))
        return(object)
    }
    # calculate stability of the clusters check if there are more than 1 k value in ks range
    stability <- NULL
    stability <- calculate_stability(metadata(object)$moSC3$consensus, k)
    
    d <- data.frame(Cluster = factor(1:length(stability)), Stability = stability)
    ggplot(d, aes(x = d$Cluster, y = d$Stability)) + geom_bar(stat = "identity") + ylim(0, 1) + 
        labs(x = "Cluster", y = "Stability Index") + theme_bw()
}

#' @rdname moSC3_plot_cluster_stability
#' @aliases moSC3_plot_cluster_stability
setMethod("moSC3_plot_cluster_stability", signature(object = "SingleCellExperiment"), moSC3_plot_cluster_stability.SingleCellExperiment)
