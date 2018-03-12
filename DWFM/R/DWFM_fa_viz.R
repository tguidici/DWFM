

#' A vector with cutpoints to be used in conjunction with fancyplots
#'
#' Takes as inputs all covariance matrices to be plotted via fancyplots, and returns a sequence of cutpoints so that there are approximately the same number of data points between each cut.
#'
#' The extreme cutpoints are located at the min and max of all covariance matrices, with intermediate cutpoints spaced approximately at
#' \code{quantile(, probs = seq(0,1,length=9))} over the negative and positive covariance values separately.
#'
#'
#' @param covmats A list of covariance matrices that will be plotted via fancyplots
#' @return a vector of 16 cutpoints.
#' @examples \dontrun{
#'  colvals <- color_break_vals(reconstructed$recon_cov)
#'  }
#' @export
color_break_vals <- function(covmats){

  cov_lims_v2 = matrix(data=NA, nrow = length(covmats), ncol = 2)
  for(i in 1:nrow(cov_lims_v2)){
    cov_lims_v2[i,] = range(covmats[[i]])
  }

  cov_percentiles_p =cov_percentiles_n = matrix(data=NA, nrow = nrow(cov_lims_v2), ncol = 9)
  for(i in 1:nrow(cov_percentiles_p)){
    cov_percentiles_p[i,] = stats::quantile(covmats[[i]][covmats[[i]]>0], probs = seq(0,1,length=9))
    cov_percentiles_n[i,] = stats::quantile(covmats[[i]][covmats[[i]]<0], probs = seq(0,1,length=9))
  }

 return(c(min(cov_lims_v2[,1]), apply(cov_percentiles_n[,2:8], 2, stats::median), apply(cov_percentiles_p[,2:8], 2, stats::median),max(cov_lims_v2[,2])))
}


#' Horizontal bar chart of pseudo_LFC values
#'
#' Produces a horizontally oriented bar chart of pseudo_LFC values, suitable for plotting with fancyplots.
#'
#' @param dta A vector or matrix with one column, containing data to be plotted
#' @param minbar most negative value across all barcharts that will appear in fancyplots
#' @param maxbar most positive value across all barcharts that will appear in fancyplots
#' @return a ggplot object suitable for plotting with fancyplots
#' @export
horiz_bar<-function(dta, minbar, maxbar){

  dta = as.data.frame(cbind(1:length(dta),dta))
  colnames(dta)=c("x", "y")

  horzplot <- ggplot2::ggplot(data=dta, ggplot2::aes(x,y))+
    ggplot2::geom_bar(stat="identity", position="identity", color="lightgray", fill="lightgrey")+
    ggplot2::labs(list(title=NULL, x=NULL,y=NULL))+
    ggplot2::scale_y_continuous(limits=c(minbar,maxbar)) +
    ggplot2::theme(legend.position="none",
          axis.text.x=ggplot2::element_blank(),
          axis.text.y=ggplot2::element_blank(),
          panel.background=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          panel.grid.major.x=ggplot2::element_blank(),
          panel.grid.minor.x=ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line(colour="gray95", size=0.25),
          plot.margin=ggplot2::unit(rep(0,4), "lines"))
  return(horzplot)
}

#' Labels for a horizontal bar chart of pseudo-Log fold change values
#'
#' Produces a horizontally oriented, empty bar chart of psuedo-LFC values, suitable for plotting with fancyplots.
#' Data is not plotted, but plot label and axis are plotted.
#'
#' @param dta A vector or matrix with one column, containing data to be plotted
#' @param minbar most negative value across all barcharts that will appear in fancyplots
#' @param maxbar most positive value across all barcharts that will appear in fancyplots
#' @param dta_name dataset name, used for labeling
#' @return a ggplot object suitable for plotting with fancyplots
#' @export
horiz_bar_label<-function(dta, dta_name, minbar, maxbar){

  dta = as.data.frame(cbind(1:length(dta),dta))
  colnames(dta)=c("x", "y")

  horzbar_lab <- ggplot2::ggplot(data=dta, ggplot2::aes(x,y))+
    ggplot2::geom_blank()+
    ggplot2::labs(list(title=dta_name, x=NULL,y=NULL))+
    ggplot2::scale_y_continuous(limits=c(minbar,maxbar)) +
    ggplot2::theme(legend.position="none",
          axis.text.x=ggplot2::element_blank(),
          panel.background=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          panel.grid.major.x=ggplot2::element_blank(),
          panel.grid.minor.x=ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line(colour="gray95", size=0.25),
          plot.margin=ggplot2::unit(rep(0,4), "lines"),
          panel.border = ggplot2::element_rect(colour = "gray", fill=NA, size=1))

  return(horzbar_lab)
}

#' Vertical bar chart of pseudo_LFC values
#'
#' Produces a vertically oriented bar chart of pseudo_LFC values, suitable for plotting with fancyplots.
#'
#' @param dta A vector or matrix with one column, containing data to be plotted
#' @param minbar most negative value across all barcharts that will appear in fancyplots
#' @param maxbar most positive value across all barcharts that will appear in fancyplots
#' @return a ggplot object suitable for plotting with fancyplots
#' @export
vert_bar <- function(dta,minbar, maxbar){
  dta = as.data.frame(cbind(length(dta):1,dta))
  colnames(dta)=c("x", "y")

  vertbar <- ggplot2::ggplot(data=dta, ggplot2::aes(x,y))+
    ggplot2::geom_bar(stat="identity", position="identity", color="lightgray", fill="lightgrey")+
    ggplot2::coord_flip() +
    ggplot2::labs(list(title=NULL, x=NULL,y=NULL))+
    ggplot2::scale_y_continuous(limits=c(minbar,maxbar)) +
    ggplot2::theme(legend.position="none",
          axis.text.x=ggplot2::element_blank(),
          axis.text.y=ggplot2::element_blank(),
          panel.background=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          panel.grid.major.y=ggplot2::element_blank(),
          panel.grid.minor.y=ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_line(colour="gray95", size=0.25),
          plot.margin=ggplot2::unit(rep(0,4), "lines"))
  return(vertbar)
}

#' Labels for a vertical bar chart of pseudo_LFC values
#'
#' Produces a vertically oriented, empty bar chart of pseudo_LFC values, suitable for plotting with fancyplots.
#' Data is not plotted, but plot label and axis are plotted.
#'
#' @param dta A vector or matrix with one column, containing data to be plotted
#' @param minbar most negative value across all barcharts that will appear in fancyplots
#' @param maxbar most positive value across all barcharts that will appear in fancyplots
#' @param dta_name dataset name, used for labeling
#' @return a ggplot object suitable for plotting with fancyplots
#' @export
vert_bar_label <-function(dta, dta_name, minbar, maxbar){
  #could just make sure data has a name; also could input c(minbar, maxbar)
  dta = as.data.frame(cbind(length(dta):1,dta))
  colnames(dta)=c("x", "y")

  vertbar_lab <- ggplot2::ggplot(data=dta, ggplot2::aes(x,y)) +
    ggplot2::geom_blank() + ggplot2::coord_flip() +
    ggplot2::labs(list(title=dta_name, x=NULL,y=NULL))+
    ggplot2::scale_y_continuous(limits=c(minbar,maxbar)) +
    ggplot2::theme(legend.position="none",
          panel.background=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          axis.title.y=ggplot2::element_blank(),
          panel.grid.major.y=ggplot2::element_blank(),
          panel.grid.minor.y=ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_line(colour="gray95", size=0.25),
          plot.margin=ggplot2::unit(rep(0,4), "lines"),
          panel.border = ggplot2::element_rect(colour = "gray", fill=NA, size=1))+
    ggplot2::scale_x_discrete(breaks=c(-1,1), labels=c("", ""))

  return(vertbar_lab)
}

#' A heatmap of a reconstructed covariance matrix
#'
#' Creates a heatmap via \code{ggplot()} and \code{cut()}.
#'
#' The continuous valued input matrix \code{dta} is converted to a matrix of factors via \code{cut()}.
#' This matrix is then plotted as a heatmap where the factor levels correspond to colors in \code{color_pal}.
#'
#' @param dta a reconstructed covariance matrix
#' @param color_pal (optional) a color palette to use in the heatmap
#' @param val_scale a vector of cutpoints to partition data in reconstructed covariance matrices, \code{length(val_scale) = length(color_pal) + 1}
#' @return a ggplot object
#' @export
ggplot_heatmap_v3 <- function(dta,color_pal = bluered_palette, val_scale){

  rownames(dta) = nrow(dta):1 #gotta reverse the order of rows, but not columns.
  colnames(dta) = 1:nrow(dta)
  melt_dta = reshape2::melt(dta)
  melt_dta[,"value"] = cut(dta, breaks=val_scale, include.lowest=TRUE)
  melt_dta[,"value"] = factor(melt_dta[,"value"], levels = rev(levels(melt_dta[,"value"])))
  colnames(melt_dta) = c("X1","X2", "value")

  heatmap <- ggplot2::ggplot(melt_dta, ggplot2::aes(y=X1, x=X2)) + ggplot2::geom_tile(ggplot2::aes(fill=value)) +
    ggplot2::scale_fill_manual(values = rev(color_pal), drop=FALSE) +
    ggplot2::labs(list(title=NULL, x=NULL,y=NULL))+
    ggplot2::theme(axis.text.y=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          legend.position="none",
          panel.background=ggplot2::element_blank(),
          plot.margin=ggplot2::unit(rep(0,4), "lines"))

  return(heatmap)
}

#' A heatmap of a reconstructed covariance matrix, with color scale
#'
#' Creates a heatmap with color scale, via \code{ggplot()} and \code{cut()}.
#'
#' The continuous valued input matrix \code{dta} is converted to a matrix of factors via \code{cut()}.
#' This matrix is then plotted as a heatmap where the factor levels correspond to colors in \code{color_pal}.
#' The color scale is plotted directly on top of the heatmap.
#' This labeled version is meant to be used for reference or further work in image processing software.
#'
#' @param dta a reconstructed covariance matrix
#' @param color_pal (optional) a color palette to use in the heatmap
#' @param val_scale a vector of cutpoints to partition data in reconstructed covariance matrices, \code{length(val_scale) = length(color_pal) + 1}
#' @param dta_name dataset name, used for labeling
#' @return a ggplot object
#' @export
ggplot_heatmap_label_v3 <- function(dta, dta_name, color_pal, val_scale){

  rownames(dta) = nrow(dta):1 #gotta reverse the order of rows, but not columns.
  colnames(dta) = 1:nrow(dta)
  melt_dta = reshape2::melt(dta)
  melt_dta[,"value"] = cut(dta, breaks=val_scale, include.lowest=TRUE)
  melt_dta[,"value"] = factor(melt_dta[,"value"], levels = rev(levels(melt_dta[,"value"])))
  melt_dta = cbind(melt_dta, "augh"=factor(levels(melt_dta[,"value"])[8],levels = rev(levels(melt_dta[,"value"])) ))
  colnames(melt_dta) = c("X1","X2", "value")

  heatmap_label <- ggplot2::ggplot(melt_dta, ggplot2::aes(y=X1, x=X2)) + ggplot2::geom_tile(ggplot2::aes(fill=value)) +
    ggplot2::scale_fill_manual(values = rev(color_pal), drop=FALSE) +
    ggplot2::labs(list(title=paste("reconstruct cov", dta_name), x=NULL,y=NULL))+
    ggplot2::theme(axis.text.y=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          legend.position=c(0.5, 0.5),
          legend.key.size = ggplot2::unit(0.75, "cm"),
          panel.background=ggplot2::element_blank(),
          plot.margin=ggplot2::unit(rep(0,4), "lines"),
          panel.border = ggplot2::element_rect(colour = "gray", fill=NA, size=1))


  return(heatmap_label)

}

#' Helper function that setups the underlying grid layout for use in fancyplots
#'
#' @param col_width A vector of column widths in inches
#' @param row_height A vector of row heights in inches
#' @export
VP_setup <- function(col_width, row_height){
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = length(row_height), ncol = length(col_width),
                                             widths=grid::unit(col_width, "inches"),
                                             heights = grid::unit(row_height, "inches"))))
}

#' Make all necessary heatmaps and bar charts for fancy plots
#'
#' Produces all necessary plots for things.
#'
#' @param pseudo_LFC output from \code{logFC()} - a matrix of pseudo_LFC values to be plotted between heatmaps
#' @param CovMat_est reconstructed covariance matrices based on estimates B^k and Q
#' @param subset_order a vector giving what variables should be plotted and in the order that they should be plotted
#' @param subset_name (optional) name to be plotted with labeled plots.
#' @param color_pal (optional) color palette. defaults to \code{bluered_palette}
#' @param val_scale a vector of cutpoints to partition data in reconstructed covariance matrices, \code{length(val_scale) = length(color_pal) + 1}
#' @return a list of lists, suitable for use with any of the layout options.
#' \code{charts} contains all necessary heatmaps and barcharts, \code{labels} contains labeled versions on the same, for internal reference
#' @examples
#' #generate covariance matrices and estimate sparse parameters
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
#' #setup for fancy plots
#' colvals = c()
#' colvals[[1]] <- color_break_vals(reconstructed_output[[1]]$recon_cov)
#' colvals[[2]] <- color_break_vals(reconstructed_output[[2]]$recon_cov)
#' fancy <- vector(mode = "list", length = 2)
#' for(i in 1:length(fancy)){
#'   fancy[[i]] = fancyplots_setup(logfc[[i]], reconstructed_output[[i]]$recon_cov,
#'   subset_order = 1:ncol(reconstructed_output[[i]]$recon_cov[[1]]), subset_name = "dta order",
#'   bluered_palette,val_scale = colvals[[i]] )
#' }
#' @export
fancyplots_setup <- function(pseudo_LFC, CovMat_est, subset_order, subset_name = NULL, color_pal = bluered_palette, val_scale ){ #subset_order should be the subset of vars we are interested in, in the correct order
  ####prep
  ratio_bar_v = ratio_bar_lab_v= ratio_bar_h = ratio_bar_lab_h= cov_maps = cov_maps_lab = c()
  bar_lims = c(min(pseudo_LFC[subset_order,]), max(pseudo_LFC[subset_order,]))


  for(i in 1:length(CovMat_est)){
    cov_maps[[i]] = ggplot_heatmap_v3(CovMat_est[[i]][subset_order,subset_order],color_pal, val_scale)
    cov_maps_lab[[i]] = ggplot_heatmap_label_v3(CovMat_est[[i]][subset_order,subset_order], paste(names(CovMat_est)[i], subset_name), color_pal, val_scale)

  }
  for(i in 1:ncol(pseudo_LFC)){
    ratio_bar_v[[i]] =vert_bar(-1*pseudo_LFC[subset_order,i], minbar=-1*bar_lims[2], maxbar=-1*bar_lims[1])
    ratio_bar_lab_v[[i]] = vert_bar_label(-1*pseudo_LFC[subset_order,i], colnames(pseudo_LFC)[i], minbar=-1*bar_lims[2], maxbar=-1*bar_lims[1])

    ratio_bar_h[[i]] = horiz_bar(pseudo_LFC[subset_order,i], minbar=bar_lims[1], maxbar=bar_lims[2])
    ratio_bar_lab_h[[i]]  = horiz_bar_label(pseudo_LFC[subset_order,i], colnames(pseudo_LFC)[i], minbar=bar_lims[1], maxbar=bar_lims[2])
  }

  return(list("charts" = list("ratio_bar_v" = ratio_bar_v,  "ratio_bar_h" = ratio_bar_h, "cov_maps" = cov_maps, "bar_lims" = bar_lims),
              "labels" = list("ratio_bar_v" = ratio_bar_lab_v,"ratio_bar_h" = ratio_bar_lab_h,"cov_maps" = cov_maps_lab, "bar_lims" = bar_lims)))
}


#' Plots fancyplots for 3x1 experimental design
#'
#' Plots pseudoLFC bar charts and reconstructed covariance matrices for data having a 3x1 experimental design setup.
#'
#' @param base_dim (optional) an integer value giving the base size for a cell in layout grid.
#' @param fancyplots_setup_obj_l2 Output from fancyplots_setup. Should be \code{fancyplots_setup_obj$charts} to plot the actual data, or \code{fancyplots_setup_obj$labels} to plot labels
#' @param orient \code{orient=1} has one heatmap in the first row and two on the second row; \code{orient=2} has two heatmaps in the first row and one on the second row
#' @param contrast_mat The contrast matrix used with \code{logFC_mat()}
#' @return fancyplots for 3x1 experimental design
#' @examples
#' #generate covariance matrices and estimate sparse parameters
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
#' #setup for fancy plots
#' colvals = c()
#' colvals[[1]] <- color_break_vals(reconstructed_output[[1]]$recon_cov)
#' colvals[[2]] <- color_break_vals(reconstructed_output[[2]]$recon_cov)
#' fancy <- vector(mode = "list", length = 2)
#' for(i in 1:length(fancy)){
#'   fancy[[i]] = fancyplots_setup(logfc[[i]], reconstructed_output[[i]]$recon_cov,
#'   subset_order = 1:ncol(reconstructed_output[[i]]$recon_cov[[1]]), subset_name = "dta order",
#'   bluered_palette,val_scale = colvals[[i]] )
#' }
#' #plot fancy plots
#' layout_3x1(base_dim = 5, fancy[[1]]$charts, orient=1, contrast_mat[[1]])
#' layout_3x1(base_dim = 5, fancy[[1]]$labels, orient=2, contrast_mat[[1]])
#' @export
layout_3x1 <- function(base_dim = 5, fancyplots_setup_obj_l2, orient, contrast_mat){
  layout_dim = c()
  layout_dim$bdim = base_dim
  layout_dim$barmult = 0.4

  barchart_order = c(which(contrast_mat[,1] == 1 & contrast_mat[,2]==2), which(contrast_mat[,1] == 1 & contrast_mat[,2]==3),
                     which(contrast_mat[,1] == 2 & contrast_mat[,2]==3) )


  grid::grid.newpage()
  vplayout <- function(x, y)
    grid::viewport(layout.pos.row = x, layout.pos.col = y)

  layout_dim$col_width = c(layout_dim$bdim,layout_dim$bdim,layout_dim$barmult*layout_dim$bdim, layout_dim$bdim)
  layout_dim$row_height = c(layout_dim$bdim, layout_dim$barmult*layout_dim$bdim, layout_dim$bdim)

  VP_setup(layout_dim$col_width, layout_dim$row_height)


  if(orient==1){
    print(fancyplots_setup_obj_l2$cov_maps[[1]], vp = vplayout(1, 2))
    print(fancyplots_setup_obj_l2$cov_maps[[2]], vp = vplayout(3, 1))
    print(fancyplots_setup_obj_l2$cov_maps[[3]], vp = vplayout(3, 4))

    print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[1]]], vp = vplayout(2, 1))
    print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[2]]], vp = vplayout(2, 4))
    print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[3]]], vp = vplayout(3, 3))

  }
  if(orient==2){
    print(fancyplots_setup_obj_l2$cov_maps[[1]], vp = vplayout(1, 1))
    print(fancyplots_setup_obj_l2$cov_maps[[2]], vp = vplayout(1, 4))
    print(fancyplots_setup_obj_l2$cov_maps[[3]], vp = vplayout(3, 2))

    print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[1]]], vp = vplayout(1, 3))
    print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[2]]], vp = vplayout(2, 1))
    print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[3]]], vp = vplayout(2, 4))

  }
}

#' Plots fancyplots for 4x1 experimental design
#'
#' Plots pseudoLFC bar charts and reconstructed covariance matrices for data having a 4x1 experimental design setup.
#'
#' @param base_dim (optional) an integer value giving the base size for a cell in layout grid.
#' @param fancyplots_setup_obj_l2 Output from fancyplots_setup. Should be \code{fancyplots_setup_obj$charts} to plot the actual data, or \code{fancyplots_setup_obj$labels} to plot labels
#' @param contrast_mat The contrast matrix used with \code{logFC_mat()}
#' @return fancyplots for 4x1 experimental design
#' @examples \dontrun{
#'  layout_4x1(base_dim = 5, fancy[[2]]$charts, contrast_mat[[2]])
#' }
#' @export
layout_4x1 <- function(base_dim = 5,fancyplots_setup_obj_l2, contrast_mat){
  layout_dim = c()
  layout_dim$bdim = base_dim
  layout_dim$barmult = 0.4

  barchart_order = c(which(contrast_mat[,1] == 1 & contrast_mat[,2]==2), which(contrast_mat[,1] == 1 & contrast_mat[,2]==3),
                     which(contrast_mat[,1] == 2 & contrast_mat[,2]==3), which(contrast_mat[,1] == 2 & contrast_mat[,2]==4),
                     which(contrast_mat[,1] == 1 & contrast_mat[,2]==4), which(contrast_mat[,1] == 3 & contrast_mat[,2]==4))

  grid::grid.newpage()
  vplayout <- function(x, y)
    grid::viewport(layout.pos.row = x, layout.pos.col = y)

  layout_dim$col_width = c(layout_dim$bdim, layout_dim$bdim, layout_dim$barmult*layout_dim$bdim, layout_dim$bdim)
  layout_dim$row_height = c(layout_dim$bdim, layout_dim$barmult*layout_dim$bdim, layout_dim$bdim,
                            layout_dim$barmult*layout_dim$bdim, layout_dim$bdim)

  VP_setup(layout_dim$col_width, layout_dim$row_height)

  print(fancyplots_setup_obj_l2$cov_maps[[1]], vp = vplayout(1, 2))
  print(fancyplots_setup_obj_l2$cov_maps[[2]], vp = vplayout(3, 1))
  print(fancyplots_setup_obj_l2$cov_maps[[3]], vp = vplayout(3, 4))
  print(fancyplots_setup_obj_l2$cov_maps[[4]], vp = vplayout(5, 2))

  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[1]]], vp = vplayout(2, 1))
  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[2]]], vp = vplayout(2, 4))

  print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[3]]], vp = vplayout(3, 3))

  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[4]]], vp = vplayout(4, 1))
  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[5]]], vp = vplayout(4, 2))
  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[6]]], vp = vplayout(4, 4))

}

#' Plots fancyplots for 2x2 experimental design
#'
#' Plots pseudoLFC bar charts and reconstructed covariance matrices for data having a 2x2 experimental design setup.
#'
#' @param base_dim (optional) an integer value giving the base size for a cell in layout grid.
#' @param fancyplots_setup_obj_l2 Output from fancyplots_setup. Should be \code{fancyplots_setup_obj$charts} to plot the actual data, or \code{fancyplots_setup_obj$labels} to plot labels
#' @param orient \code{orient=1} arranges the heatmaps in a square; \code{orient=2} arranges the heatmaps in a V
#' @param contrast_mat The contrast matrix used with \code{logFC_mat()}
#' @return fancyplots for 2x2 experimental design
#' @examples
#' #generate covariance matrices and estimate sparse parameters
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
#' #setup for fancy plots
#' colvals = c()
#' colvals[[1]] <- color_break_vals(reconstructed_output[[1]]$recon_cov)
#' colvals[[2]] <- color_break_vals(reconstructed_output[[2]]$recon_cov)
#' fancy <- vector(mode = "list", length = 2)
#' for(i in 1:length(fancy)){
#'   fancy[[i]] = fancyplots_setup(logfc[[i]], reconstructed_output[[i]]$recon_cov,
#'   subset_order = 1:ncol(reconstructed_output[[i]]$recon_cov[[1]]), subset_name = "dta order",
#'   bluered_palette,val_scale = colvals[[i]] )
#' }
#' #plot fancy plots
#' layout_2x2(base_dim = 5, fancy[[2]]$charts, orient = 1, contrast_mat[[2]])
#' layout_2x2(base_dim = 5, fancy[[2]]$labels, orient = 2, contrast_mat[[2]])
#' @export
layout_2x2<- function(base_dim = 5, fancyplots_setup_obj_l2, orient, contrast_mat){
  layout_dim = c()
  layout_dim$bdim = base_dim
  layout_dim$barmult = 0.4

  barchart_order = c(which(contrast_mat[,1] == 1 & contrast_mat[,2]==2), which(contrast_mat[,1] == 1 & contrast_mat[,2]==3),
                     which(contrast_mat[,1] == 2 & contrast_mat[,2]==4), which(contrast_mat[,1] == 3 & contrast_mat[,2]==4) )


  grid::grid.newpage()
  vplayout <- function(x, y)
    grid::viewport(layout.pos.row = x, layout.pos.col = y)

  if(orient==1){
    layout_dim$col_width = c(layout_dim$bdim,layout_dim$barmult*layout_dim$bdim, layout_dim$bdim)
    layout_dim$row_height = c(layout_dim$bdim, layout_dim$barmult*layout_dim$bdim, layout_dim$bdim)

    VP_setup(layout_dim$col_width, layout_dim$row_height)

    print(fancyplots_setup_obj_l2$cov_maps[[1]], vp = vplayout(1, 1))
    print(fancyplots_setup_obj_l2$cov_maps[[2]], vp = vplayout(1, 3))
    print(fancyplots_setup_obj_l2$cov_maps[[3]], vp = vplayout(3, 1))
    print(fancyplots_setup_obj_l2$cov_maps[[4]], vp = vplayout(3, 3))


    print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[1]]], vp = vplayout(1, 2))
    print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[2]]], vp = vplayout(2, 1))
    print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[3]]], vp = vplayout(2, 3))
    print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[4]]], vp = vplayout(3, 2))

  }

  if(orient == 2){

    layout_dim$col_width = c(layout_dim$bdim,layout_dim$bdim,layout_dim$barmult*layout_dim$bdim, layout_dim$bdim,layout_dim$bdim)
    layout_dim$row_height = c(layout_dim$bdim, layout_dim$barmult*layout_dim$bdim, layout_dim$bdim)

    VP_setup(layout_dim$col_width, layout_dim$row_height)

    print(fancyplots_setup_obj_l2$cov_maps[[1]], vp = vplayout(1, 1))
    print(fancyplots_setup_obj_l2$cov_maps[[2]], vp = vplayout(1, 5))
    print(fancyplots_setup_obj_l2$cov_maps[[3]], vp = vplayout(3, 2))
    print(fancyplots_setup_obj_l2$cov_maps[[4]], vp = vplayout(3, 4))

    print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[1]]], vp = vplayout(1, 3))
    print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[2]]], vp = vplayout(2, 1))
    print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[3]]], vp = vplayout(2, 5))
    print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[4]]], vp = vplayout(3, 3))
  }
}


#' Plots fancyplots for 3x2 experimental design
#'
#' Plots pseudoLFC bar charts and reconstructed covariance matrices for data having a 3x2 experimental design setup.
#'
#' @param base_dim (optional) an integer value giving the base size for a cell in layout grid.
#' @param fancyplots_setup_obj_l2 Output from fancyplots_setup. Should be \code{fancyplots_setup_obj$charts} to plot the actual data, or \code{fancyplots_setup_obj$labels} to plot labels
#' @param contrast_mat The contrast matrix used with \code{logFC_mat()}
#' @return fancyplots for 3x2 experimental design
#' @examples  \dontrun{ layout_3x2(base_dim = 6,fancy[[2]]$charts, contrast_mat[[2]] )}
#' @export
layout_3x2 <- function(base_dim=3,fancyplots_setup_obj_l2, contrast_mat){
  layout_dim = c()
  layout_dim$bdim = base_dim
  layout_dim$barmult = 0.4

  barchart_order = c(which(contrast_mat[,1] == 1 & contrast_mat[,2]==2), which(contrast_mat[,1] == 1 & contrast_mat[,2]==3),
                     which(contrast_mat[,1] == 2 & contrast_mat[,2]==4), which(contrast_mat[,1] == 3 & contrast_mat[,2]==4),
                     which(contrast_mat[,1] == 1 & contrast_mat[,2]==5), which(contrast_mat[,1] == 3 & contrast_mat[,2]==5),
                     which(contrast_mat[,1] == 4 & contrast_mat[,2]==6), which(contrast_mat[,1] == 2 & contrast_mat[,2]==6),
                     which(contrast_mat[,1] == 5 & contrast_mat[,2]==6))


  grid::grid.newpage()
  vplayout <- function(x, y)
    grid::viewport(layout.pos.row = x, layout.pos.col = y)

  layout_dim$col_width = c(layout_dim$bdim,layout_dim$bdim,layout_dim$bdim, layout_dim$barmult*layout_dim$bdim,
                           layout_dim$bdim, layout_dim$bdim,layout_dim$bdim)
  layout_dim$row_height = c(layout_dim$bdim, layout_dim$barmult*layout_dim$bdim, layout_dim$bdim,
                            layout_dim$barmult*layout_dim$bdim, layout_dim$bdim)

  VP_setup(layout_dim$col_width, layout_dim$row_height)

  print(fancyplots_setup_obj_l2$cov_maps[[1]], vp = vplayout(1, 1))
  print(fancyplots_setup_obj_l2$cov_maps[[2]], vp = vplayout(1, 7))
  print(fancyplots_setup_obj_l2$cov_maps[[3]], vp = vplayout(3, 2))
  print(fancyplots_setup_obj_l2$cov_maps[[4]], vp = vplayout(3, 6))
  print(fancyplots_setup_obj_l2$cov_maps[[5]], vp = vplayout(5, 3))
  print(fancyplots_setup_obj_l2$cov_maps[[6]], vp = vplayout(5, 5))

  print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[1]]], vp = vplayout(1, 4)) #center
  #first row of bars
  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[2]]], vp = vplayout(2, 1))
  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[3]]], vp = vplayout(2, 7))

  print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[4]]], vp = vplayout(3, 4)) #center
  #second row of bars
  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[5]]], vp = vplayout(4, 1))
  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[6]]], vp = vplayout(4, 7))
  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[7]]], vp = vplayout(4, 2))
  print(fancyplots_setup_obj_l2$ratio_bar_h[[barchart_order[8]]], vp = vplayout(4, 6))

  print(fancyplots_setup_obj_l2$ratio_bar_v[[barchart_order[9]]], vp = vplayout(5, 4))#center

}




#' Plots fancyplots for 2x1 experimental design
#'
#' Plots pseudoLFC bar charts and reconstructed covariance matrices for data having a 2x1 experimental design setup.
#'
#' @param base_dim (optional) an integer value giving the base size for a cell in layout grid.
#' @param fancyplots_setup_obj_l2 Output from fancyplots_setup. Should be \code{fancyplots_setup_obj$charts} to plot the actual data, or \code{fancyplots_setup_obj$labels} to plot labels
#' @param orient \code{orient=1} has both heatmaps in one row; \code{orient=2} has both heatmaps in one column.
#' @return fancyplots for 2x1 experimental design
#' @examples  \dontrun{layout_2x1(fancy[[2]]$charts, orient = 1)}
#' @export
layout_2x1 <- function(base_dim=5,fancyplots_setup_obj_l2, orient){
  layout_dim = c()
  layout_dim$bdim = base_dim
  layout_dim$barmult = 0.4

  grid::grid.newpage()
  vplayout <- function(x, y)
    grid::viewport(layout.pos.row = x, layout.pos.col = y)

  if(orient=="1"){
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 3, widths=grid::unit(c(layout_dim$bdim,layout_dim$barmult*layout_dim$bdim, layout_dim$bdim), "inches"),
                                               heights = grid::unit(c(layout_dim$bdim), "inches"))))
    print(fancyplots_setup_obj_l2$cov_maps[[1]], vp = vplayout(1, 1))
    print(fancyplots_setup_obj_l2$ratio_bar_v[[1]], vp = vplayout(1, 2))
    print(fancyplots_setup_obj_l2$cov_maps[[2]], vp = vplayout(1, 3))
  }

  if(orient=="2"){
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 3, ncol = 1, heights=grid::unit(c(layout_dim$bdim,layout_dim$barmult*layout_dim$bdim,layout_dim$bdim), "inches"),
                                               widths = grid::unit(c(layout_dim$bdim), "inches"))))
    print(fancyplots_setup_obj_l2$cov_maps[[1]], vp = vplayout(1, 1))
    print(fancyplots_setup_obj_l2$ratio_bar_h[[1]], vp = vplayout(2, 1))
    print(fancyplots_setup_obj_l2$cov_maps[[2]], vp = vplayout(3, 1))
  }
}
