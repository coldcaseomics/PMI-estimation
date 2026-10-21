# Source 01 and 02 beforehand. Do not commit local sample-level objects.
# Start with integer(0) to avoid expensive unintentional reruns.
folds_to_run <- integer(0)
if (length(folds_to_run) == 0L) {
  message("No folds selected. Set folds_to_run <- 1:5 explicitly to run.")
} else {
  dir.create("results/generated",recursive=TRUE,showWarnings=FALSE)
  for (fold in folds_to_run) {
    result <- run_pmi_fold(fold,meta_model,rse_final,fold_id)
    # Aggregate summaries only; no sample-level metadata stored in GitHub.
    summary <- data.frame(
      fold=result$fold,train=result$n_train,test=result$n_test,
      ref=result$reference,expression_genes=result$n_expression_genes,
      r_gt_02=result$n_r02,FDR_lt_05=result$n_fdr,
      candidate_genes=result$n_candidates)
    utils::write.csv(summary,
      sprintf("results/generated/fold_%d_summary.csv",fold),row.names=FALSE)
    utils::write.csv(result$panels,
      sprintf("results/generated/fold_%d_panels.csv",fold),row.names=FALSE)
  }
}
