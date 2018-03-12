## ----step0, echo = TRUE--------------------------------------------------
library(DWFM)

covmats_etc <- vector(mode="list", length=2)
covmats_etc[[1]] <- step0(X[[1]], 2)
covmats_etc[[2]] <- step0(X[[2]], 2)


## ----Procrustes rotation,echo=TRUE---------------------------------------
lambdahat_proc <- vector(mode="list", length=2)
lambdahat_proc[[1]] <- covmats_etc[[1]]$lambdahat
lambdahat_proc[[2]] <- covmats_etc[[2]]$lambdahat

for(j in 1:length(covmats_etc)){
  for(i in 1:dim(covmats_etc[[j]]$lambdahat)[3]){
    lambdahat_proc[[j]][,,i] <- vegan::procrustes(Lambda_vals[[j]][,,i], covmats_etc[[j]]$lambdahat[,,i], scale=FALSE, symmetric=FALSE)$Yrot
  }
}

## ----sparse lambdahat, echo = TRUE---------------------------------------
sparse_params <- vector(mode="list", length=2)
for(i in 1:length(sparse_params)){
  sparse_params[[i]] <- sparse_lambdahat_cardinality(covmats_etc[[i]])
}

spca_cutoff = 0.01
sps_lambdahat_proc = lambdahat_proc
for(j in 1:length(lambdahat_proc)){
  sps_lambdahat_proc$C[[j]] <- trunc_C(lambdahat_proc[[j]], sparse_params[[j]]$col_card)
  sps_lambdahat_proc$M[[j]] <- trunc_M(lambdahat_proc[[j]])
  sps_lambdahat_proc$SPCA[[j]] <- array(data=NA, dim = dim(sps_lambdahat_proc$C[[j]]))
  for(i in 1:dim(lambdahat_proc[[j]])[3]){
    temp <- elasticnet::spca(covmats_etc[[j]]$covmat[[i]],K = sparse_params[[j]]$factsize_est, type = "Gram", sparse = "varnum", para =sparse_params[[j]]$col_card )
    temp$loadings[abs(temp$loadings)<spca_cutoff] <-0
    sps_lambdahat_proc$SPCA[[j]][,,i] <- temp$loadings%*%diag(sqrt(temp$var.all*temp$pev), ncol=ncol(temp$loadings))
  }
}


## ----run algorithm, echo = TRUE------------------------------------------
test_algo <- vector(mode="list", length=2)
test_algo[[1]] <- estimateBkQ(sps_lambdahat_proc$C[[1]], constraint = 1, qadj = "n")
test_algo[[2]] <- estimateBkQ(sps_lambdahat_proc$C[[2]], constraint = 3, norm_groups = list(c(1,2), c(3,4) ), qadj = "n")

#ok, good.
test_algo[[1]]$small_load
test_algo[[2]]$small_load

## ----generate reconstructed covariance matrices and pseudo-LFC for plotting, echo = TRUE----

reconstructed <- vector(mode="list", length=2)
reconstructed[[1]] = recon_covlh(test_algo[[1]])
reconstructed[[2]] = recon_covlh(test_algo[[2]])

contrast_mat <- logfc <- vector(mode="list", length=2)
contrast_mat[[1]] <- matrix(data = c(1,2,1,3,2,3), ncol = 2, byrow=TRUE)
contrast_mat[[2]] <- matrix(data = c(1,2,3,4,1,3,2,4), ncol=2, byrow = TRUE)
logfc[[1]] <- logFC_mat(test_algo[[1]]$Bhat,contrasts=contrast_mat[[1]], condit_names= c("c1", "c2", "c3"))
logfc[[2]] <- logFC_mat(test_algo[[2]]$Bhat,contrasts=contrast_mat[[2]], condit_names= c("c1", "c2", "c3", "c4"))

## ---- echo = TRUE--------------------------------------------------------
colvals = c()
colvals[[1]] <- color_break_vals(reconstructed[[1]]$recon_cov)
colvals[[2]] <- color_break_vals(reconstructed[[2]]$recon_cov)

colvals[[1]][8:9] <- colvals[[1]][8:9]/8
colvals[[2]][8:9] <- colvals[[2]][8:9]/8

## ---- echo = TRUE--------------------------------------------------------

fancy <- vector(mode = "list", length = 2)
for(i in 1:length(fancy)){
  fancy[[i]] = fancyplots_setup(logfc[[i]], reconstructed[[i]]$recon_cov, subset_order = 1:ncol(reconstructed[[i]]$recon_cov[[1]]), subset_name = "dta order", bluered_palette,val_scale = colvals[[i]] )
}

## ---- echo = TRUE--------------------------------------------------------
pdf("fancyplots.pdf", width = 22 , height = 20)
layout_3x1(base_dim = 5, fancy[[1]]$charts, orient=1, contrast_mat[[1]])
layout_3x1(base_dim = 5, fancy[[1]]$charts, orient=2, contrast_mat[[1]])
layout_3x1(base_dim = 5, fancy[[1]]$labels, orient=2, contrast_mat[[1]])
layout_2x2(base_dim = 5, fancy[[2]]$charts, orient = 1, contrast_mat[[2]])
layout_2x2(base_dim = 5, fancy[[2]]$charts, orient = 2, contrast_mat[[2]])
layout_4x1(base_dim = 5, fancy[[2]]$charts, contrast_mat[[2]])
dev.off()


## ---- echo = TRUE--------------------------------------------------------
sqrt(crossprod(as.vector(reconstructed[[1]]$recon_lh) - as.vector(Lambda_vals[[1]])))/sqrt(crossprod(as.vector(Lambda_vals[[1]])))
sqrt(crossprod(as.vector(reconstructed[[2]]$recon_lh) - as.vector(Lambda_vals[[2]])))/sqrt(crossprod(as.vector(Lambda_vals[[2]])))

