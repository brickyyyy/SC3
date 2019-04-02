# setClass("SingleCellExperimentCollection",
#          contains="SingleCellExperiment",
#          validity = function(object) {
#            msg <- NULL
#            if (!all(sapply(object, is, "SingleCellExperiment")))
#              msg <- c(msg, "members must all be 'SingleCellExperiment' classes")
#            tryCatch({
#              if (anyDuplicated(names(object)))
#                msg <- c(msg, "each setName must be distinct")
#            }, error=function(err) {
#              msg <<- c(msg, conditionMessage(err))
#            })
#            if (!is.null(msg))
#              msg
#            else
#              TRUE
#          })
# 
# setMethod("SingleCellExperimentCollection",
#           signature=signature(
#             object="list",
#             ks = NULL,
#             gene_filter = TRUE, pct_dropout_min = 10, pct_dropout_max = 90,
#             d_region_min = 0.04, d_region_max = 0.07, svm_num_cells = NULL,
#             svm_train_inds = NULL, svm_max = 5000, n_cores = NULL, kmeans_nstart = NULL,
#             kmeans_iter_max = 1e+09, k_estimator = FALSE, biology = FALSE, rand_seed = 1),
#           function(object, ..., ks, gene_filter, pct_dropout_min, pct_dropout_max,
#                    d_region_min, d_region_max, svm_num_cells, svm_train_inds, svm_max, n_cores, kmeans_nstart,
#                    kmeans_iter_max, k_estimator, biology, rand_seed) {
#             new("SingleCellExperimentCollection", object)
#           })
# 
# setGeneric("SingleCellExperimentCollection",
#            function(object)
#              standardGeneric("SingleCellExperimentCollection"))