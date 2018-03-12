
#K = number of conditions (conditsize)
#n = vector, with number of samples per condition (length(n) = K)
#p = number of variables
#m = number of factors (factsize)

#' Calculate the angle between two vectors
#'
#' This function calculates the angle, in degrees, between two vectors.
#'
#' Uses \code{acos()}, which will sometimes return \code{NaN} if the vectors are very close to 180 degrees apart.
#' @param x1 A vector or matrix with one column
#' @param x2 same as x1, same length
#' @return The angle between \code{x1} and \code{x2}.
#' @examples
#' angle(c(0,1), c(1,0))
#' angle(c(1.5, 1.5), c(-1,-1))
#' @export
angle <-function(x1,x2){
    return((180/pi)*acos(crossprod(x1, x2)/(sqrt(crossprod(x1)) * sqrt(crossprod(x2)))))
}

#' Preliminary computations for \code{estimateBkQ}
#'
#' Returns a list of covariance matrices, non-sparse lambdahats and other useful preliminary items.
#'
#' @param dta a list of datasets you wish to use in the covma_fa algorithm
#' @param m how many columns do you want in your scaled eigen decomposition?
#' @return A list containing covariance matrices, non-sparse lambdahats (from the scaled eigen-decomposition),
#' a matrix with percent variance explained by the first \code{m} eigen vectors
#' and results from the eigen decomposition of the covariance of each dataset.
#' @examples
#' #X[[1]] contains simulated data from a 3x1 experimental design scenario
#' #X[[2]] contains simulated data from a 2x2 experimental design scenario
#'
#' covmats_etc <- vector(mode="list", length=2)
#' covmats_etc[[1]] <- step0(X[[1]], 2)
#' covmats_etc[[2]] <- step0(X[[2]], 2)
#' @export
step0 <- function(dta, m = 10){
  covmat <- eigen_res <- c()
  percent_var <- matrix(data=NA, nrow = length(dta), ncol = m)
  lambdahat <- array(data=0, dim=c(nrow(dta[[1]]), m, length(dta)))

  for(i in 1:length(dta)){
    covmat[[i]] = stats::cov(t(dta[[i]]))
    eigen_res[[i]] = eigen(covmat[[i]])
    for(j in 1:m){
      percent_var[i,j] <- eigen_res[[i]]$values[j]/sum(eigen_res[[i]]$values)
    }
    lambdahat[,,i] <- eigen_res[[i]]$vectors[,1:m]%*%diag(sqrt(eigen_res[[i]]$values[1:m]), nrow=m)
    rownames(covmat[[i]]) = colnames(covmat[[i]]) = rownames(dta[[1]])
  }
  rownames(lambdahat) = rownames(dta[[1]])

  return(list("covmat" = covmat,"lambdahat"=lambdahat,"percent_var"= percent_var, "eigen_res" = eigen_res))
}

#' Estimate the number and cardinality of latent factors
#'
#' Estimates the number and cardinality of latent factors by running the leading eigenvector community detection algorithm.
#'
#' Transforms each covariance matrix into a graph, and uses \code{igraph::leading.eigenvector.community()}, to detect distinct communities.
#' The smallest number of communities in one covariance matrix is returned as the estimated number of factors.
#' The sizes of each community are returned as the target column cardinalities.
#'
#' @param step0_obj output from \code{step0()}
#' @return A list with the estimated number of factors (factsize_est) and the number of non-zero elements in each factor (col_card)
#' @examples
#' #generate the covariance matrices from example datasets
#' covmats_etc <- vector(mode="list", length=2)
#' covmats_etc[[1]] <- step0(X[[1]], 2)
#' covmats_etc[[2]] <- step0(X[[2]], 2)
#' #calculate the number of factors and number of non-sparse elements for each factor
#' sparse_params <- vector(mode="list", length=2)
#' for(i in 1:length(sparse_params)){
#'  sparse_params[[i]] <- sparse_lambdahat_cardinality(covmats_etc[[i]])
#'  }
#' @export
sparse_lambdahat_cardinality <-function(step0_obj){

  K = dim(step0_obj$lambdahat)[3]

  cov_graph = c()
  le = c()
  comm_numb = c()
  CommSizes <- matrix(data=0, nrow = dim(step0_obj$lambdahat)[1], ncol= K)

  #find the number and sizes of communities from leading eigenvector community detection
  for(i in 1:K){
    cov_graph[[i]] = igraph::graph.adjacency(adjmatrix= step0_obj$covmat[[i]], mode="undirected", weighted =TRUE ,diag =TRUE)
    le[[i]] = as.list(tryCatch(igraph::leading.eigenvector.community(cov_graph[[i]], steps = -1, weights = abs(igraph::E(cov_graph[[i]])$weight),
                                                                     start = NULL, options =list(maxiter=7000),
                                                                     callback = NULL, extra = NULL, env = parent.frame()),
                               error=function(c){NA}))
    comm_numb[i] <- length(le[[i]])
    CommSizes[1:comm_numb[i],i] <- sort(igraph::sizes(le[[i]]), decreasing=TRUE)
  }
  factsize_est = min(comm_numb)
  nzero_elts = sort(igraph::sizes(le[[which.min(comm_numb)]]), decreasing=TRUE)
  return(list("factsize_est"=factsize_est,"col_card"=nzero_elts))
}

#' Adjusts the orientation of a set of matrices for maximum similarity.
#'
#' If matrices have a single column, it checks whether multiplying non-reference conditions by \code{-1} decreases the angle between the reference condition and another condition.
#' If matrices have more than a single column, it performs a procrustes rotation for maximum similarity.
#'
#' Uses \code{angle()} which uses \code{acos()}, and therefore may return an error if \code{m=1} and \code{lambdahat[,,k]} and \code{lambdahat[,,j]} are very close to 180 degrees different.
#' This is relatively unlikely to happen.
#'
#' @param lambdahat An array of dimension \code{(p,m,K)}
#' @param ref_condit Which of the \code{K} matrices (slices of the array) should be used as reference? (ie: what should everything else be rotated to?)
#' @param m how many factors in the dataset?
#' @return An array matching the dimensions of \code{lambdahat}
#' @examples \dontrun{
#' test_lh <- adjust_orient(covmats_etc$lambdahat, ref_condit = 1, m=2)
#' }
#' @export
adjust_orient <- function(lambdahat, ref_condit = 1,m){
  lambdahat_proc <- lambdahat
  if(m==1){
    K = ncol(lambdahat)
    for(i in (1:K)[-ref_condit]){
      a1 = angle(lambdahat[,ref_condit],lambdahat[,i])
      a2 = angle(lambdahat[,ref_condit],(-1)*lambdahat[,i])
      if(a2 < a1){
        lambdahat_proc[,i] <- (-1)*lambdahat_proc[,i]
      }
    }
  }
  if(m >1){
    K = dim(lambdahat)[3]
    for(i in (1:K)[-ref_condit]){
        lambdahat_proc[,,i] <-vegan::procrustes(lambdahat[,,ref_condit], lambdahat[,,i], scale=FALSE, symmetric=FALSE)$Yrot
    }
  }
  return(lambdahat_proc)
}


#' Truncates the columns of an array by cardinality.
#'
#' Sets to zero the elements of an array so that each column has a specified number of non-zero elements.
#'
#' @param lambdahat An array
#' @param col_card A vector of length \code{dim(lambdahat)[2]}
#' @return Array lambdahat where \code{lambdahat[,i,k]} has \code{col_card[i]} non-zero elements for every \code{k}.
#' @examples
#' sparse_params <- covmats_etc <- vector(mode="list", length=2)
#' for(i in 1:length(sparse_params)){
#'  covmats_etc[[i]] <- step0(X[[i]],2)
#'  sparse_params[[i]] <- sparse_lambdahat_cardinality(covmats_etc[[i]])
#' }
#' #we need to rotate all lambdahats to the true Lambda values before sparsifying
#' lambdahat_proc <- vector(mode="list", length=2)
#' lambdahat_proc[[1]] <- covmats_etc[[1]]$lambdahat
#' lambdahat_proc[[2]] <- covmats_etc[[2]]$lambdahat
#' for(j in 1:length(covmats_etc)){
#'  for(i in 1:dim(covmats_etc[[j]]$lambdahat)[3]){
#'   lambdahat_proc[[j]][,,i] <- vegan::procrustes(Lambda_vals[[j]][,,i],
#'   covmats_etc[[j]]$lambdahat[,,i], scale=FALSE, symmetric=FALSE)$Yrot
#'  }
#' }
#' #sparsify
#' spca_cutoff = 0.01
#' sps_lambdahat_proc = lambdahat_proc
#' for(j in 1:length(lambdahat_proc)){
#' sps_lambdahat_proc$C[[j]] <- trunc_C(lambdahat_proc[[j]], sparse_params[[j]]$col_card)
#' }
#' @export
trunc_C <- function(lambdahat, col_card){
  for(i in 1:dim(lambdahat)[3]){
    for(j in 1:dim(lambdahat)[2]){
      ids <- order(abs(lambdahat[,j,i]))[1:(dim(lambdahat)[1]-col_card[j])]
      lambdahat[ids,j,i] <- 0
    }
  }
  return(lambdahat)
}


#' Truncates an array by magnitude
#'
#' Sets to zero any entries in an array having absolute value below some cutoff.
#'
#' @param lambdahat An array
#' @param cutoff A cutoff value
#' @return lambdahat where all values between \code{-cutoff} and \code{cutoff} are set to 0
#' @examples
#' covmats_etc <- vector(mode="list", length=2)
#' for(i in 1:length(covmats_etc)){
#'  covmats_etc[[i]] <- step0(X[[i]],2)
#' }
#' #we need to rotate all lambdahats to the ture Lambda values before sparsifying
#' lambdahat_proc <- vector(mode="list", length=2)
#' lambdahat_proc[[1]] <- covmats_etc[[1]]$lambdahat
#' lambdahat_proc[[2]] <- covmats_etc[[2]]$lambdahat
#' for(j in 1:length(covmats_etc)){
#'  for(i in 1:dim(covmats_etc[[j]]$lambdahat)[3]){
#'   lambdahat_proc[[j]][,,i] <- vegan::procrustes(Lambda_vals[[j]][,,i],
#'   covmats_etc[[j]]$lambdahat[,,i], scale=FALSE, symmetric=FALSE)$Yrot
#'  }
#' }
#' #sparsify
#' sps_lambdahat_proc = lambdahat_proc
#' for(j in 1:length(lambdahat_proc)){
#' sps_lambdahat_proc$M[[j]] <- trunc_M(lambdahat_proc[[j]])
#' }
#' @export
trunc_M <- function(lambdahat, cutoff = 0.1){
  lambdahat[abs(lambdahat) < cutoff] <- 0
  return(lambdahat)
}


#' Estimates B^k and Q for covma_fa algorithm from corresponding paper
#'
#' Returns a list which includes intermediate and final estimates for \code{Q} and \code{Bhat^k}
#'
#' @param lambdahat a sparse or non sparse array with one lambdahat for each condition
#' @param constraint what normalization constraint do you want? choices are 1 and 3
#' @param norm_groups if constraint==3, you need to specify how you want the B^k values to be normalized
#' @param qadj Do you want to include the s_ij scaling step?
#' @param small_load_cutoff Do you want to remove variables which load <= some small value on all columns of Q? Defaults to 0.
#' @return B^k, Q, a list of the variables removed due to small loadings, and intermediate values of B^k and Q if applicable.
#' Bhat and Qhat are the final values, Bhat_v1 and Qhat_v1 are intermediate values.
#' @examples
#' sparse_params <- covmats_etc <- vector(mode="list", length=2)
#' for(i in 1:length(sparse_params)){
#'  covmats_etc[[i]] <- step0(X[[i]],2)
#'  sparse_params[[i]] <- sparse_lambdahat_cardinality(covmats_etc[[i]])
#' }
#' #we need to rotate all lambdahats to the ture Lambda values before sparsifying
#' lambdahat_proc <- vector(mode="list", length=2)
#' lambdahat_proc[[1]] <- covmats_etc[[1]]$lambdahat
#' lambdahat_proc[[2]] <- covmats_etc[[2]]$lambdahat
#' for(j in 1:length(covmats_etc)){
#'  for(i in 1:dim(covmats_etc[[j]]$lambdahat)[3]){
#'   lambdahat_proc[[j]][,,i] <- vegan::procrustes(Lambda_vals[[j]][,,i],
#'   covmats_etc[[j]]$lambdahat[,,i], scale=FALSE, symmetric=FALSE)$Yrot
#'  }
#' }
#' #sparsify
#' spca_cutoff = 0.01
#' sps_lambdahat_proc = lambdahat_proc
#' for(j in 1:length(lambdahat_proc)){
#' sps_lambdahat_proc$C[[j]] <- trunc_C(lambdahat_proc[[j]], sparse_params[[j]]$col_card)
#' }
#' test_algo <- vector(mode="list", length=2)
#' test_algo[[1]] <- estimateBkQ(sps_lambdahat_proc$C[[1]], constraint = 1, qadj = "n")
#' test_algo[[2]] <- estimateBkQ(sps_lambdahat_proc$C[[2]], constraint = 3,
#' norm_groups = list(c(1,2), c(3,4) ), qadj = "n")
#' @export
estimateBkQ <- function(lambdahat, constraint, norm_groups = NULL, qadj, small_load_cutoff=0){
  factsize = dim(lambdahat)[2]
  conditsize = dim(lambdahat)[3]

  if(constraint==1){
    design_levels=1
    norm_groups = list(1:dim(lambdahat)[3])
  }
  if(constraint==3){
    design_levels=2
  }

  Qhat_v1 <- estimate_Q(lambdahat, design_levels, qadj)

  Bhat_v1 <- estimate_B(Qhat_v1, lambdahat, norm_groups, factsize, design_levels)

  negB <- which(rowSums(Bhat_v1<0)!=0) #do we have any negative values in Bhat?

  if(length(negB)>0){
    lambdahat <- update_lambda(lambdahat, negB)
    Qhat_v2 <- estimate_Q(lambdahat, design_levels, qadj)
    Bhat_v2 <- estimate_B(Qhat_v2, lambdahat, norm_groups, factsize, design_levels)
    rownames(Bhat_v2) <- rownames(Qhat_v2) <- rownames(lambdahat)

    t_id = which(rowSums(abs(Qhat_v2) <= small_load_cutoff) == factsize) #check for small load vars

    return(list("Qhat" = Qhat_v2[-t_id,], "Bhat"=Bhat_v2[-t_id,],
                "Bhat_v1" = Bhat_v1, "Qhat_v1" = Qhat_v1, "small_load" = rownames(Qhat_v2)[t_id]))
  }

  if(length(negB)==0){ #no negative values to worry about!
    rownames(Bhat_v1) <- rownames(Qhat_v1) <- rownames(lambdahat)
    t_id = which(rowSums(abs(Qhat_v1) <= small_load_cutoff) == factsize) #check for small load vars
    if(length(t_id) >0 ){
      return(list("Qhat" = Qhat_v1[-t_id,], "Bhat"=Bhat_v1[-t_id,],"small_load" = rownames(Qhat_v1)[t_id]))
    } else{
      return(list("Qhat" = Qhat_v1, "Bhat"=Bhat_v1,"small_load" = "nothing too small"))
    }
  }

}


#' Estimates Qhat via covmat_fa algorithm
#'
#' @param lambdahat An array of dimensions \code{(p,m,K)}
#' @param design_levels parameter indicating normalization structure.
#' @param qadj "y" or "n", indicating whether \code{Qhat} should be adjusted w/ the scaling factor {s_ij}
#' @return estimated Qhat
#' @export
estimate_Q <-function(lambdahat, design_levels, qadj){
  #re-estimate Qhat
  Qhat <- (1/design_levels)*apply(lambdahat, c(1,2), sum)

  if(qadj=="y"){
    #rescale as necessary
    scale_mult = matrix(data=dim(lambdahat)[3], nrow = nrow(Qhat), ncol=ncol(Qhat))

    for(i in 1:ncol(scale_mult)){
      scale_mult[,i] <- scale_mult[,i]/rowSums(lambdahat[,i,]!=0)
    }
    scale_mult[scale_mult=="Inf"] <- 0

    for(i in 1:ncol(Qhat)){
      Qhat[,i] <- Qhat[,i]*scale_mult[,i]
    }
  }


 return(Qhat)
}

#' Estimates Bhat via covmat_fa algorithm
#'
#' @param Qhat Estimated Qhat from application of covma_fa algorithm
#' @param lambdahat An array of dimensions \code{(p,m,K)}
#' @param norm_groups Which of the \code{K} conditions should be normalized together?
#' @param factsize Number of latent factors (same as \code{m})
#' @param design_levels 1 or 2, depending on the experimental design. see vignette for details
#' @return estimated Bhat
#' @export
estimate_B <-function(Qhat, lambdahat, norm_groups, factsize, design_levels){
  B_summands = array(data=0, dim=c(nrow(Qhat), dim(lambdahat)[2],dim(lambdahat)[3]))

  for(i in 1:dim(B_summands)[1]){
    for(j in 1:dim(B_summands)[3]){
      for(z in 1:dim(B_summands)[2]){
        if(Qhat[i,z]!=0){
          B_summands[i,z,j] = lambdahat[i,z,j]/Qhat[i,z]
        }
      }
    }
  }

  Bhat = (1/factsize)*apply(B_summands, c(1,3), sum)
  #normalize taking the exp design into account
  for(i in 1:length(norm_groups)){
    nz1 = which(rowSums(Bhat[,norm_groups[[i]]]==0)!=ncol(Bhat)/design_levels)
    Bhat[nz1,norm_groups[[i]]] = Bhat[nz1,norm_groups[[i]]]/rowSums(Bhat[nz1,norm_groups[[i]]])
  }

  return(Bhat)
}

#' Updates lambdahat to remove entries that lead to negative values in Bhat
#'
#' @param lambdahat An array of dimensions \code{(p,m,K)}
#' @param negB parameter indicating normalization structure.
#' @return updated lambdahat, where noisy entries are removed
#' @export
update_lambda <-function(lambdahat, negB){
  factsize = ncol(lambdahat)
  for(i in 1:length(negB)){
    for(z in 1:factsize){
      neg = sum(lambdahat[negB[i],z,]<0)
      pos = sum(lambdahat[negB[i],z,]>0)

      if(neg==pos){
        lambdahat[negB[i],z,] <- 0
      }

      if(neg>pos){#set positive to 0, fact1
        id = which(lambdahat[negB[i],z,]>0)
        lambdahat[negB[i],z,id] <- 0
      }
      if(pos>neg){#set negative to 0, fact1
        id = which(lambdahat[negB[i],z,]<0)
        lambdahat[negB[i],z,id] <- 0
      }

    }
  }
  return(lambdahat)
}

#' Calculates the pseudo-log fold change between all 2-way comparisons of interest
#'
#'
#' Calculates \code{log2(Bhat[i,j]/Bhat[i,k])} where \code{j,k} are specified via \code{contrasts} matrix
#' @param Bhat Bhat matrix from the output of \code{estimateBkQ()}
#' @param contrasts a matrix with 2 columns and n rows, specifying the \code{n} total 2-way comparisons of interest.
#' @param condit_names Names of each of the conditions/datasets of interest.
#' @return A matrix with one column for each row of \code{contrasts}. Note that \code{Nan}, \code{Inf}, and \code{-Inf} values are set to \code{0}
#' @examples
#' sparse_params <- covmats_etc <- vector(mode="list", length=2)
#' for(i in 1:length(sparse_params)){
#'  covmats_etc[[i]] <- step0(X[[i]],2)
#'  sparse_params[[i]] <- sparse_lambdahat_cardinality(covmats_etc[[i]])
#' }
#' #we need to rotate all lambdahats to the ture Lambda values before sparsifying
#' lambdahat_proc <- vector(mode="list", length=2)
#' lambdahat_proc[[1]] <- covmats_etc[[1]]$lambdahat
#' lambdahat_proc[[2]] <- covmats_etc[[2]]$lambdahat
#' for(j in 1:length(covmats_etc)){
#'  for(i in 1:dim(covmats_etc[[j]]$lambdahat)[3]){
#'   lambdahat_proc[[j]][,,i] <- vegan::procrustes(Lambda_vals[[j]][,,i],
#'   covmats_etc[[j]]$lambdahat[,,i], scale=FALSE, symmetric=FALSE)$Yrot
#'  }
#' }
#' #sparsify
#' spca_cutoff = 0.01
#' sps_lambdahat_proc = lambdahat_proc
#' for(j in 1:length(lambdahat_proc)){
#' sps_lambdahat_proc$C[[j]] <- trunc_C(lambdahat_proc[[j]], sparse_params[[j]]$col_card)
#' }
#' #run the algorithm
#' test_algo <- vector(mode="list", length=2)
#' test_algo[[1]] <- estimateBkQ(sps_lambdahat_proc$C[[1]], constraint = 1, qadj = "n")
#' test_algo[[2]] <- estimateBkQ(sps_lambdahat_proc$C[[2]], constraint = 3,
#' norm_groups = list(c(1,2), c(3,4) ), qadj = "n")
#' #create reconstructed covariance and lambdahat matrices
#' reconstructed_output <- vector(mode="list", length=2)
#' for(i in 1:length(reconstructed_output)){
#'  reconstructed_output[[i]] <- recon_covlh(test_algo[[i]])
#' }
#' #calculate psuedo-LFC values
#' contrast_mat <- logfc <- vector(mode="list", length=2)
#' contrast_mat[[1]] <- matrix(data = c(1,2,1,3,2,3), ncol = 2, byrow=TRUE)
#' contrast_mat[[2]] <- matrix(data = c(1,2,3,4,1,3,2,4), ncol=2, byrow = TRUE)
#' logfc[[1]] <- logFC_mat(test_algo[[1]]$Bhat,contrasts=contrast_mat[[1]],
#' condit_names= c("c1", "c2", "c3"))
#' logfc[[2]] <- logFC_mat(test_algo[[2]]$Bhat,contrasts=contrast_mat[[2]],
#' condit_names= c("c1", "c2", "c3", "c4"))
#' @export
logFC_mat <- function(Bhat, contrasts, condit_names){
  logFC_mat <- matrix(data=NA, nrow = nrow(Bhat), ncol = nrow(contrasts))
  cn = c()
  for(i in 1:ncol(logFC_mat)){
    logFC_mat[,i] = suppressWarnings(log2(Bhat[,contrasts[i,1]]/Bhat[,contrasts[i,2]]))
    cn[i] = paste(condit_names[contrasts[i,1]], condit_names[contrasts[i,2]], sep="_")
  }
  colnames(logFC_mat) = cn
  logFC_mat[logFC_mat %in% c("NaN", "Inf", "-Inf")] <- 0
  return(logFC_mat)
}

#' Calculates reconstructed lambdahat and covariance matrices based on the output of \code{estimateBkQ()}.
#'
#' @param method_res The output from \code{estimateBkQ()}.
#' @param row_subset Optional subset of rows to use in reconstruction. If \code{NULL}, all rows are used.
#' @return A list containing  reconstructed covariance matrices (recon_cov), and reconstructed lambdahats (recon_lh).
#' @examples
#' sparse_params <- covmats_etc <- vector(mode="list", length=2)
#' for(i in 1:length(sparse_params)){
#'  covmats_etc[[i]] <- step0(X[[i]],2)
#'  sparse_params[[i]] <- sparse_lambdahat_cardinality(covmats_etc[[i]])
#' }
#' #we need to rotate all lambdahats to the ture Lambda values before sparsifying
#' lambdahat_proc <- vector(mode="list", length=2)
#' lambdahat_proc[[1]] <- covmats_etc[[1]]$lambdahat
#' lambdahat_proc[[2]] <- covmats_etc[[2]]$lambdahat
#' for(j in 1:length(covmats_etc)){
#'  for(i in 1:dim(covmats_etc[[j]]$lambdahat)[3]){
#'   lambdahat_proc[[j]][,,i] <- vegan::procrustes(Lambda_vals[[j]][,,i],
#'   covmats_etc[[j]]$lambdahat[,,i], scale=FALSE, symmetric=FALSE)$Yrot
#'  }
#' }
#' #sparsify
#' spca_cutoff = 0.01
#' sps_lambdahat_proc = lambdahat_proc
#' for(j in 1:length(lambdahat_proc)){
#' sps_lambdahat_proc$C[[j]] <- trunc_C(lambdahat_proc[[j]], sparse_params[[j]]$col_card)
#' }
#' test_algo <- vector(mode="list", length=2)
#' test_algo[[1]] <- estimateBkQ(sps_lambdahat_proc$C[[1]], constraint = 1, qadj = "n")
#' test_algo[[2]] <- estimateBkQ(sps_lambdahat_proc$C[[2]], constraint = 3,
#' norm_groups = list(c(1,2), c(3,4) ), qadj = "n")
#' reconstructed_output <- vector(mode="list", length=2)
#' for(i in 1:length(reconstructed_output)){
#'  reconstructed_output[[i]] <- recon_covlh(test_algo[[i]])
#' }
#' @export
recon_covlh <- function(method_res, row_subset=NULL){
  if(is.null(row_subset)){
    row_subset = nrow(method_res$Bhat)
  }
  K= ncol(method_res$Bhat)

  recon_lh = array(data=NA, dim=c(nrow(method_res$Qhat),ncol(method_res$Qhat), K) )
  recon_cov= c()

  for(i in 1:K){
    recon_lh[,,i] <- diag(method_res$Bhat[,i])%*%method_res$Qhat

    recon_cov[[i]] = recon_lh[,,i]%*%t(recon_lh[,,i])

    colnames(recon_cov[[i]]) =  rownames(recon_cov[[i]]) = rownames(recon_lh) = rownames(method_res$Qhat)

  }

  names(recon_cov) = dimnames(recon_lh)[[3]]
  return(list("recon_cov" =recon_cov , "recon_lh" = recon_lh))
}


