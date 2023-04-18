#'Running the liana tool to retrieve LR interactions
#'
#'This function runs the liana_wrap function to retrieve LR interactions for the given seurat object using the
#'method from CellPhoneDB and the liana Consensus database.
#'
#'
#'@param seurat_object seurat object with expression data
#'@param condition_ident name of the meta data field containing condition information to split the object
#'@param annotation_ident name of the meta data field containing the cluster annotation
#'@param organism either human or mouse
#'@param out_path output path to save the results to
#'@import Seurat liana dplyr stringr
#'@return list of LR interaction tables by condition
#' @export
run_LR_analysis <- function(seurat_object, condition_ident, annotation_ident, organism, out_path) {

  Idents(seurat_object) <- condition_ident
  res_list = list()
  for (ident in levels(seurat_object)) {
    ident_rename <- str_replace_all(ident, "[,;.:-]", "_")
    sub_obj <- subset(x = seurat_object, idents = ident)

    if(organism=="human"){
      liana_cpdb_res <- liana_wrap(sub_obj, resource = c("Consensus"), idents_col = annotation_ident, method = "cellphonedb", expr_prop = 0, permutation.params =
      list(
        nperms = 1000
      ))
    } else if (organism=="mouse"){
      liana_cpdb_res <- liana_wrap(sub_obj, resource = c("MouseConsensus"), idents_col = annotation_ident, method = "cellphonedb", expr_prop = 0, permutation.params =
      list(
        nperms = 1000
      ))
    } else {
      print("Unsopported organism defined. Currently only human and mouse data can be analyzed.")
      stop("exit execution of LR analysis")
    }

    res_table <- liana_cpdb_res %>%
      filter(pvalue < 0.05)

    res_table <- select(res_table, source, target, ligand, receptor, lr.mean)

    colnames(res_table)[colnames(res_table) == "ligand"] <- "gene_A"
    colnames(res_table)[colnames(res_table) == "receptor"] <- "gene_B"
    colnames(res_table)[colnames(res_table) == "lr.mean"] <- "MeanLR"

    res_table["type_gene_A"] <- "Ligand"
    res_table["type_gene_B"] <- "Receptor"

    res_table <- select(res_table, source, target, gene_A, gene_B, type_gene_A, type_gene_B, MeanLR)

    write.csv(res_table, file = paste0(out_path, "/", ident_rename, "_LR_table.csv"), row.names = FALSE)
    res_list[[ident]] <- res_table
  }
  return(res_list)
}


#'Running the LR to TF analysis from the LR2TF package
#'
#'This function runs a TF activity analysis with the dorothea tool and combines active TFs with Ligand and Receptors
#'based on expression and a-priori knowledge from the Omnipath database
#'
#'@param seurat_object seurat object with expression data
#'@param parameters list of arguments that are needed to run the TF analysis
#'@param LR_results tables with LR interactions by consitions
#'@import Seurat rcompanion LR2TF
#'@return list of LRTF interaction tables by condition
#' @export
run_TF_analysis <- function(seurat_object, parameters, LR_results = NA) {

  if (is.null(parameters[["comparison_list"]])) {
    comparison_list = NA
  }else {
    comparison_list = parameters[["comparison_list"]]
  }
  if (is.null(parameters[["log2fc"]])) {
    log2fc = 1.0
  } else {
    log2fc = parameters[["log2fc"]]
  }
  if (is.null(parameters[["pval"]])) {
    pval = 0.05
  } else {
    pval = parameters[["pval"]]
  }


  all_results = LR2TF::dorothea_tf_prediction(seurat_object,
                                              parameters[["out_path"]],
                                              parameters[["confidence_level"]],
                                              parameters[["organism"]],
                                              parameters[["condition_annotation"]],
                                              parameters[["celltype_annotation"]],
                                              comparison_list,
                                              log2fc,
                                              pval)

  Idents(seurat_object) <- parameters[["condition_annotation"]]
  res_list <- list()
  for (con in levels(seurat_object)) {
    con_rename <- str_replace_all(con, "[,;.:-]", "_")
    if (parameters[["organism"]] == "human") {
      if (is.null(all_results[[con_rename]][["condition"]])) {
        result_table <- LR2TF::generate_CrossTalkeR_input_significant_table(all_results[[con_rename]][["cluster"]],
                                                                            parameters[["confidence_level"]],
                                                                            all_results[[paste0(con_rename, "_average_expression")]])
      } else {
        result_table <- LR2TF::generate_CrossTalkeR_input_significant_table(all_results[[con_rename]][["condition"]],
                                                                            parameters[["confidence_level"]],
                                                                            all_results[[paste0(con_rename, "_average_expression")]])
      }
    } else if (parameters[["organism"]] == "mouse") {
      if (is.null(all_results[[con_rename]][["condition"]])) {
        result_table <- LR2TF::generate_CrossTalkeR_input_mouse_significant_table(all_results[[con_rename]][["cluster"]],
                                                                                  parameters[["confidence_level"]],
                                                                                  all_results[[paste0(con_rename, "_average_expression")]])
      } else {
        result_table <- LR2TF::generate_CrossTalkeR_input_mouse_significant_table(all_results[[con_rename]][["condition"]],
                                                                                  parameters[["confidence_level"]],
                                                                                  all_results[[paste0(con_rename, "_average_expression")]])
      }
    } else {
      print("Unsopported organism defined. Currently only human and mouse data can be analyzed.")
      stop("exit execution of TF analysis")
    }


    if (!anyNA(LR_results[[con]])) {
      res_list[[con]] <- LR2TF::combine_LR_and_TF(result_table,
                                                  LR_results[[con]],
                                                  parameters[["out_path"]],
                                                  con_rename,
                                                  add_node_type = FALSE)
    } else {
      res_list[[con]] <- result_table
    }
  }
  return(res_list)
}


#'Running the LR and LR to TF analysis
#'
#' Wrapper to run the LR analysis with liana and the TF analysis with the LR2TF package
#'
#'@param seurat_object seurat object with expression data
#'@param parameters list of arguments that are needed to run the LR and TF analyses
#'@import Seurat rcompanion LR2TF
#'@return list of LRTF interaction tables by condition that can be used to run CrossTalkeR
#' @export
run_LRTF_analysis <- function(seurat_object, parameters){
  LR_tables <- run_LR_analysis(seurat_object, parameters[["condition_annotation"]], parameters[["celltype_annotation"]], parameters[["organism"]], parameters[["out_path"]])
  CTR_input <- run_TF_analysis(seurat_object, parameters, LR_tables)

  return(CTR_input)
}
