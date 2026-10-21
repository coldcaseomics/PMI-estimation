# Independently reconstructed fold workflow. Not the original saved RStudio script.
# R packages: edgeR, glmnet, SummarizedExperiment.
pmi_formula <- PMI ~ Age + I(Age^2) + RIN + I(RIN^2) + Sex +
  Manner_Of_Death + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5
covariate_formula <- ~ Age + I(Age^2) + RIN + I(RIN^2) + Sex +
  Manner_Of_Death + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5
panel_sizes <- c(5,10,20,30,50,70,100,150,200,300,400,500,750,1000)

normalize_fold <- function(counts_train, counts_test) {
  stopifnot(requireNamespace("edgeR", quietly=TRUE))
  keep <- rowSums(edgeR::cpm(counts_train) > 1) >= ceiling(ncol(counts_train)*0.20)
  tr <- counts_train[keep,,drop=FALSE]
  te <- counts_test[keep,,drop=FALSE]
  uq <- apply(edgeR::cpm(tr),2,stats::quantile,probs=.75)
  ref_i <- which.min(abs(uq-mean(uq)))
  dtr <- edgeR::calcNormFactors(edgeR::DGEList(counts=tr),
                               method="TMM", refColumn=as.integer(ref_i))
  xtr <- edgeR::cpm(dtr,log=TRUE,prior.count=1)

  # Assumption, not author-supplied test normalization code: pairwise
  # test TMM to TRAIN reference, scaled to train reference norm factor.
  ref_count <- tr[,ref_i,drop=FALSE]
  ref_lib <- colSums(tr)[ref_i]
  ref_nf <- dtr$samples$norm.factors[ref_i]
  test_nf <- vapply(seq_len(ncol(te)),function(i){
    paired <- cbind(ref_count,te[,i,drop=FALSE])
    dp <- edgeR::calcNormFactors(edgeR::DGEList(counts=paired),
                                 method="TMM",refColumn=1L)
    dp$samples$norm.factors[2]/dp$samples$norm.factors[1]*ref_nf
  }, numeric(1))
  dte <- edgeR::DGEList(counts=te)
  dte$samples$norm.factors <- test_nf
  xte <- edgeR::cpm(dte,log=TRUE,prior.count=1)
  stopifnot(identical(rownames(xtr),rownames(xte)))
  list(train=xtr,test=xte,reference=colnames(tr)[ref_i],n_genes=nrow(tr))
}

screen_fold <- function(train_logcpm, meta_train) {
  pmi_fit <- stats::lm(pmi_formula,data=meta_train)
  pmi_r <- stats::residuals(pmi_fit)
  Xcov <- stats::model.matrix(covariate_formula,data=meta_train)
  gene_r <- stats::lm.fit(Xcov,t(train_logcpm))$residuals
  stopifnot(length(pmi_r)==nrow(gene_r))
  r <- as.numeric(stats::cor(pmi_r,gene_r,method="pearson"))
  df <- stats::df.residual(pmi_fit)-1
  p <- 2*stats::pt(-abs(r*sqrt(df/(1-r^2))),df=df)
  fdr <- stats::p.adjust(p,method="BH")
  ok <- is.finite(r) & is.finite(fdr) & (fdr < .05) & (abs(r) > .20)
  out <- data.frame(gene=rownames(train_logcpm),r=r,p=p,fdr=fdr,
                    candidate=ok,stringsAsFactors=FALSE)
  list(table=out,candidates=out[ok,,drop=FALSE],
       n_r02=sum(is.finite(r)&abs(r)>.20),n_fdr=sum(fdr<.05,na.rm=TRUE))
}

bootstrap_rank <- function(train_logcpm,meta_train,candidates,
                           nboot=100,seed=123) {
  stopifnot(requireNamespace("glmnet",quietly=TRUE))
  gene <- candidates$gene
  if (!length(gene)) return(candidates)
  xg <- t(train_logcpm[gene,,drop=FALSE])
  xc <- stats::model.matrix(covariate_formula,data=meta_train)[,-1,drop=FALSE]
  x <- cbind(xg,xc)
  penalties <- c(rep(1,ncol(xg)),rep(0,ncol(xc)))
  y <- meta_train$PMI
  counts <- integer(length(gene))
  set.seed(seed)
  for(b in seq_len(nboot)) {
    idx <- sample.int(length(y),replace=TRUE)
    cv <- glmnet::cv.glmnet(x[idx,,drop=FALSE],y[idx],alpha=.5,
                           penalty.factor=penalties,nfolds=5,
                           standardize=TRUE)
    co <- as.matrix(stats::coef(cv,s="lambda.min"))
    counts <- counts + as.integer(co[gene,1] != 0)
  }
  ranked <- candidates
  ranked$selection_frequency <- counts
  ranked <- ranked[order(-ranked$selection_frequency,ranked$p,
                         -abs(ranked$r)),,drop=FALSE]
  rownames(ranked) <- NULL
  ranked
}

evaluate_panels <- function(train_logcpm,test_logcpm,meta_train,meta_test,
                            ranked,sizes=panel_sizes,seed=123) {
  stopifnot(requireNamespace("glmnet",quietly=TRUE))
  out <- data.frame(panel_size=sizes,MAE=NA_real_,RMSE=NA_real_,
                    status="insufficient candidates")
  for (i in seq_along(sizes)) {
    k <- sizes[i]
    if (nrow(ranked)<k) next
    genes <- ranked$gene[seq_len(k)]
    xtr <- t(train_logcpm[genes,,drop=FALSE])
    xte <- t(test_logcpm[genes,,drop=FALSE])
    set.seed(seed+i)
    fit <- glmnet::cv.glmnet(xtr,meta_train$PMI,alpha=.5,nfolds=5,
                             standardize=TRUE)
    pred <- as.numeric(stats::predict(fit,newx=xte,s="lambda.min"))
    err <- pred-meta_test$PMI
    out$MAE[i] <- mean(abs(err))
    out$RMSE[i] <- sqrt(mean(err^2))
    out$status[i] <- "evaluated"
  }
  out
}

run_pmi_fold <- function(fold, meta_model, rse_final, fold_id,
                         evaluate=TRUE,seed=123) {
  train_idx <- which(fold_id != fold)
  test_idx <- which(fold_id == fold)
  meta_train <- meta_model[train_idx,,drop=FALSE]
  meta_test <- meta_model[test_idx,,drop=FALSE]
  ct <- SummarizedExperiment::assay(rse_final)
  norm <- normalize_fold(ct[,train_idx,drop=FALSE],ct[,test_idx,drop=FALSE])
  screen <- screen_fold(norm$train,meta_train)
  n_candidates <- nrow(screen$candidates)
  cat(sprintf("Fold %d: %d candidates, %d abs(r)>.2; ref %s\n",
              fold,n_candidates,screen$n_r02,norm$reference))
  ranking <- if(n_candidates>=1) bootstrap_rank(
    norm$train,meta_train,screen$candidates,seed=seed+fold) else screen$candidates
  panels <- if(evaluate) evaluate_panels(
    norm$train,norm$test,meta_train,meta_test,ranking,seed=seed+fold) else NULL
  list(fold=fold,n_train=length(train_idx),n_test=length(test_idx),
       reference=norm$reference,n_expression_genes=norm$n_genes,
       n_r02=screen$n_r02,n_fdr=screen$n_fdr,n_candidates=n_candidates,
       gene_screen=screen$table,ranking=ranking,panels=panels)
}
