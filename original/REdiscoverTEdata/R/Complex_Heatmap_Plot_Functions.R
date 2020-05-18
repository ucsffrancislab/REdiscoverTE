#' Create subfigures for Figure 3
#'
#' @param data 
#' @param star 
#'
#' @return
#' @importFrom ComplexHeatmap Heatmap
#' @export
#'
#' @examples
plotDEHeat = function(data,star){
  min_DE = min(as.vector(data))
  max_DE = max(as.vector(data))
  range_DE = seq(min(min_DE,0-max_DE),max(abs(min_DE),abs(max_DE)),length.out=100)
  col.pal_DE = colorRamp2(range_DE,colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
  
  DE_h = Heatmap(data,
                 col = col.pal_DE,
                 row_title = "",
                 row_title_side = "left",
                 column_title = "Differential Expression of TEs",
                 column_title_gp = gpar(fontsize = 18),
                 column_title_side = "top",
                 cluster_columns = T,
                 cluster_rows = F,
                 row_title_gp = gpar(fontsize = 18),
                 row_names_gp = gpar(fontsize=18),
                 row_names_side = "left",
                 column_names_gp = gpar(fontsize=18),
                 cell_fun = function(j,i,x,y,width,height,fill){
                   grid.text(star[i,j],x,y,gp=gpar(fontsize=15))
                 },
                 rect_gp = gpar(col = "white", lwd = 2),
                 heatmap_legend_param = list(title = "logFC",legend_direction="horizontal"))
  DE_h
}


#' Create subfigures for Figure 3
#'
#' @param data 
#'
#' @return
#' @importFrom ComplexHeatmap Heatmap
#' @export
#'
#' @examples
plotBetaHeat = function(data){
  min_bVals = min(as.vector(data))
  max_bVals = max(as.vector(data))
  range_bVals = seq(min(min_bVals,0-max_bVals),max(abs(min_bVals),abs(max_bVals)),length.out=100)
  col.pal_bVals = colorRamp2(range_bVals,colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
  
  delta_beta_h = Heatmap(data,
                         col = col.pal_bVals,
                         row_title = "repName",
                         row_title_side = "left",
                         show_row_names = F,
                         column_title = paste("Δβ(Tumor-Normal)\n (500bp ± 5' TEs)"),
                         column_title_side = "top",
                         column_title_gp = gpar(fontsize = 18),
                         row_names_gp = gpar(fontsize=18),
                         column_names_gp = gpar(fontsize=18),
                         cluster_columns = F,
                         rect_gp = gpar(col = "white", lwd = 2),
                         heatmap_legend_param = list(title = "Δβ(T-N)",
                                                     legend_direction="horizontal"))
}


#' Create subfigures for Figure 3
#'
#' @param data 
#' @param star 
#'
#' @return
#' @importFrom ComplexHeatmap Heatmap
#' @export
#'
#' @examples
plotCorHeat = function(data,star){
  min_cor = min(as.vector(data))
  max_cor = max(as.vector(data))
  range_cor = seq(min(min_cor,0-max_cor),max(abs(min_cor),abs(max_cor)),length.out=100)
  col.pal_cor = colorRamp2(range_cor,colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
  
  cor_h = Heatmap(data,
                  col = col.pal_cor,
                  row_title = "repName",
                  row_title_side = "left",
                  show_row_names = F,
                  column_title = "Cor: Expression vs. Methylation\n(500bp ± 5' TEs)",
                  column_title_gp = gpar(fontsize = 18),
                  column_title_side = "top",
                  cluster_rows =  F,
                  cluster_columns = F,
                  clustering_method_columns = "complete",
                  row_names_gp = gpar(fontsize=18),
                  column_names_gp = gpar(fontsize=18),
                  cell_fun = function(j,i,x,y,width,height,fill){
                    grid.text(star[i,j],x,y,gp=gpar(fontsize=15))
                  },
                  rect_gp = gpar(col = "white", lwd = 2),
                  heatmap_legend_param = list(title = "Cor",legend_direction="horizontal"))
  
  cor_h
}