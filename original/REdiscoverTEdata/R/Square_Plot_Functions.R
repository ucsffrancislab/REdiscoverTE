#' Create subfigures for Figure 3
#'
#' @param data 
#' @param element a string denoting a TE subfamily name 
#'
#' @return a ggplot2 object
#' @import ggplot2
#' @export
#'
#' @examples
exp_plot = function(data,element){
  col = c("tumor" = "#B2182B","normal" = "#2166AC")
  g = ggplot(data,aes(x=as.factor(status),y=logCPM,colour=status)) + geom_boxplot(outlier.alpha = 0, fill="white") 
  g = g + geom_point(data=data,aes(x=as.factor(status),y=logCPM,colour=status,shape=paired),
                     position=position_jitter(width=0.2),alpha=0.75,size=3) 
  # g = g + geom_line(aes(group=patient_ID),colour="gray",linetype=2)
  g = g + scale_colour_manual(values=col) + scale_shape_manual(values=c("TRUE"=19,"FALSE"=1))
  g = g + theme_classic() + labs(x="Status",y="logCPM",
                                 title=paste0(element," Expression"))
  g = g + theme(legend.position = "left",
                #legend.justification = c(0,1),
                #legend.title = element_blank(),
                legend.text = element_text(size=28),
                legend.title = element_text(size=28),
                aspect.ratio = 1,
                plot.margin = margin(t=10,r=0,b=30,l=0,unit="pt"),
                axis.title.y = element_text(size=32,colour="black"),
                axis.title.x = element_text(size=32,colour="black"),
                #axis.title.x = element_blank(),
                axis.text = element_text(size=26,colour="black",angle=0),
                axis.text.x = element_text(size=26,colour="black",angle=90,hjust=1,vjust=0.5),
                #axis.text.x = element_blank(),
                plot.title = element_text(size=32,face="bold"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
  g 
}

#' Create subfigures for Figure 3
#'
#' @param data 
#' @param element a string denoting a TE subfamily name 
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
methy_plot = function(data,element){
  col = c("tumor" = "#B2182B","normal" = "#2166AC")
  g = ggplot(data,aes(x=status,y=beta,fill=status)) + geom_violin(scale="width",trim=FALSE) + geom_boxplot(width=0.05,fill="white",outlier.alpha = 0)
  g = g + scale_fill_manual(values=col) + ylim(0,1)
  g = g + theme_classic() + labs(x="Status",y="Methylation Beta",
                                 title=paste0("Methylation +/- 500bp ",element))
  g = g + theme(legend.position = "none",
                #legend.justification = c(0,1),
                legend.title = element_blank(),
                #legend.key.size = unit(0.0085,"npc"),
                legend.text = element_text(size=22),
                aspect.ratio = 1,
                plot.margin = margin(t=10,r=0,b=30,l=0,unit="pt"),
                axis.title.y = element_text(size=32,colour="black"),
                axis.title.x = element_text(size=32,colour="black"),
                #axis.title.x = element_blank(),
                axis.text = element_text(size=26,colour="black",angle=0),
                axis.text.x = element_text(size=26,colour="black",angle=90,hjust=1,vjust=0.5),
                #axis.text.x = element_blank(),
                plot.title = element_text(size=32,face="bold"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
  g 
}

#' Create subfigures for Figure 3
#'
#' @param data 
#' @param element a string denoting a TE subfamily name 
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
cor_plot = function(data,element){
  ordered_data = data[with(data,order(data$rescale_dis)),]
  bin = cut(ordered_data$rescale_dis,breaks=seq(-500,600,20))
  window = as.character(bin)
  window = sapply(window,function(x)strsplit(x,split=",",fixed=T)[[1]][1])
  window = substr(window, start=2,stop=nchar(window))
  ordered_data$window = as.numeric(window)
  tmp = ordered_data[1,]
  tmp$delta = mean(ordered_data$cor)
  
  g = ggplot(data=ordered_data,aes(x=window,y=cor)) + geom_smooth(se =T,method="loess",colour="#4292C6",
                                                                  size=0.75,span=0.2,alpha=0.25,fullrange=F) 
  g = g + geom_rect(data=tmp,aes(xmin=0,xmax=100,ymin=-Inf,ymax=Inf),fill="#DEEBF7",alpha=0.3,colour=NA)
  g = g + labs(x="Distance to TE",y="Cor (Exp vs. M Value)",
               title=paste0("Spatial Cor:",element,"\n Exp vs. Methy"))+ theme_classic()
  g = g + scale_x_continuous(breaks=c(-500,-200,-100,0,100,200,300,600),
                             labels=c("up_5kb","up_2kb","up_1kb","TE_start","TE_end","downn_1kb","down_2kb","down_5kb"))
  g = g + geom_hline(yintercept=0,colour="gray50",linetype=2,size=1)
  g = g + theme(aspect.ratio = 1,
                axis.title.y = element_text(size=32,colour="black"),
                axis.title.x = element_text(size=32,colour="black"),
                #axis.title.x = element_blank(),
                axis.text = element_text(size=26,colour="black",angle=0),
                axis.text.x = element_text(size=22,colour="black",angle=90,hjust=1,vjust=0.5),
                #axis.text.x = element_blank(),
                plot.title = element_text(size=32,face="bold"),
                panel.grid.major.x = element_line(linetype="dashed",colour="gray50",size=0.2),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank())
  g
}

#' Create subfigures for Figure 3
#'
#' @param data 
#' @param element a string denoting a TE subfamily name 
#'
#' @return ggplot object
#' @export
#'
#' @examples
cor_density_plot = function(data,element){
  g_hist = ggplot(data,aes(distance,colour=type,linetype=type,size=type)) 
  g_hist = g_hist + geom_rect(data=data[1:2,],aes(xmin=0,xmax=100,ymin=-Inf,ymax=Inf),fill="#DEEBF7",alpha=0.2,colour=NA) + geom_freqpoly(binwidth=30)
  g_hist = g_hist + scale_size_manual(values=c("negative_sig_cor"=1,"positive_sig_cor"=1,"non_significant"=0.5))
  g_hist = g_hist + scale_linetype_manual(values=c("negative_sig_cor"="solid","positive_sig_cor"="solid","non_significant"="dashed"))
  g_hist = g_hist + theme_classic() + labs(x="Distance to REs",y="Counts: Significant Correlated CpGs",title=paste0("Sig Correlated CpGs around ",element))
  g_hist = g_hist + scale_x_continuous(breaks=c(-500,-200,-100,0,100,200,300,600),
                                       labels=c("up_5kb","up_2kb","up_1kb","TE_start","TE_end","downn_1kb","down_2kb","down_5kb"))
  g_hist = g_hist + scale_colour_manual(values=c("negative_sig_cor"="#4393C3","positive_sig_cor"="#D6604D","non_significant"="gray"),
                                        breaks=c("negative_sig_cor","positive_sig_cor","non_significant"),
                                        labels=c("sig_negative_correlation","sig_positive_correlation","not_significant"))
  
  g_hist = g_hist + theme(legend.position = "right",
                          #legend.justification = c(1,1),
                          legend.title = element_blank(),
                          #legend.key.size = unit(0.0085,"npc"),
                          legend.text = element_text(size=18),
                          aspect.ratio = 1,
                          axis.title.y = element_text(size=18,colour="black",face="bold"),
                          axis.title.x = element_text(size=18,colour="black",face="bold"),
                          #axis.title.x = element_blank(),
                          axis.text = element_text(size=18,colour="black",angle=0),
                          axis.text.x = element_text(size=16,colour="black",angle=90,hjust=1,vjust=0.5),
                          #axis.text.x = element_blank(),
                          title = element_text(size=18,face="bold"),
                          panel.grid.major = element_line(linetype="dashed",colour="grey",size=0.1),
                          panel.grid.minor = element_blank()) + guides(size=FALSE,linetype=FALSE)
  g_hist
}

#' Create subfigures for Figure 3
#'
#' @param data 
#' @param element a string denoting a TE subfamily name 
#'
#' @return a ggplot object
#' @export
#'
#' @examples
delta_beta_plot = function(data,element){
  ordered_data = data[with(data,order(data$rescale_dis)),]
  bin = cut(ordered_data$rescale_dis,breaks=seq(-500,600,20))
  window = as.character(bin)
  window = sapply(window,function(x)strsplit(x,split=",",fixed=T)[[1]][1])
  window = substr(window, start=2,stop=nchar(window))
  ordered_data$window = as.numeric(window)
  tmp = ordered_data[1,]
  tmp$delta = mean(ordered_data$delta)
  
  g = ggplot(data=ordered_data,aes(x=window,y=delta))  + geom_smooth(se =T,method="loess",colour="#4292C6",size=0.75,span=0.2,alpha=0.25) + theme_classic()
  g = g + geom_rect(data=tmp,aes(xmin=0,xmax=100,ymin=-Inf,ymax=Inf),fill="#DEEBF7",alpha=0.3,colour=NA)
  g = g + labs(x="Distance to REs",y="Delta Beta Value (T-N)",title=paste0("Beta Value Change around ",element))
  g = g + scale_x_continuous(breaks=c(-500,-200,-100,0,100,200,300,600),
                             labels=c("up_5kb","up_2kb","up_1kb","TE_start","TE_end","downn_1kb","down_2kb","down_5kb"))
  g = g + theme(aspect.ratio = 1,
                axis.title.y = element_text(size=32,colour="black"),
                axis.title.x = element_text(size=32,colour="black"),
                #axis.title.x = element_blank(),
                axis.text = element_text(size=26,colour="black",angle=0),
                axis.text.x = element_text(size=22,colour="black",angle=90,hjust=1,vjust=0.5),
                #axis.text.x = element_blank(),
                plot.title = element_text(size=32,face="bold"),
                panel.grid.major = element_line(linetype="dashed",colour="grey",size=0.1),
                panel.grid.minor = element_blank())
  g
}

#' Create subfigures for Figure 3
#'
#' @param data 
#' @param element a string denoting a TE subfamily name 
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
plot_DMP_density = function(data,element){
  g_hist = ggplot(data,aes(rescale_dis,colour=CpG_Type,size=CpG_Type,linetype=CpG_Type)) 
  g_hist = g_hist + geom_rect(data=data[1:2,],aes(xmin=0,xmax=100,ymin=-Inf,ymax=Inf),fill="#DEEBF7",alpha=0.3,colour=NA) + geom_freqpoly(binwidth=30)
  g_hist = g_hist + scale_size_manual(values=c("hypo"=1,"hyper"=1,"null"=0.5))
  g_hist = g_hist + scale_linetype_manual(values=c("hypo"="solid","hyper"="solid","null"="dashed"))
  g_hist = g_hist + theme_classic() + labs(x="Distance to TE",y="Counts: DMCs",
                                           title=paste0("DMCs Distribution \n5kb around ",element))
  g_hist = g_hist + scale_x_continuous(breaks=c(-500,-200,-100,0,100,200,300,600),
                                       labels=c("up_5kb","up_2kb","up_1kb","TE_start","TE_end","downn_1kb","down_2kb","down_5kb"))
  g_hist = g_hist + scale_colour_manual(values=c("hypo"="#4D9221","hyper"="#C51B7D","null"="gray40"),
                                        breaks=c("hypo","hyper","null"),
                                        labels=c("hypomethylation","hypermethyltaion","no change"))
  g_hist = g_hist + theme(legend.position = "right",
                          #legend.justification = c(1,1),
                          legend.title = element_blank(),
                          #legend.key.size = unit(0.0085,"npc"),
                          legend.text = element_text(size=28),
                          aspect.ratio = 1,
                          axis.title.y = element_text(size=32,colour="black"),
                          axis.title.x = element_text(size=32,colour="black"),
                          #axis.title.x = element_blank(),
                          axis.text = element_text(size=26,colour="black",angle=0),
                          axis.text.x = element_text(size=22,colour="black",angle=90,hjust=1,vjust=0.5),
                          #axis.text.x = element_blank(),
                          plot.title = element_text(size=32,face="bold"),
                          #title = element_blank(),
                          panel.grid.major.x = element_line(linetype="dashed",colour="gray50",size=0.2),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank()) + guides(size=FALSE,linetype=FALSE)
  g_hist
}

#' Create subfigures for Figure 3
#'
#' @param data 
#' @param element a string denoting a TE subfamily name 
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
scatter_cor_plot = function(data,element){
  col = c("tumor" = "#B2182B","normal" = "#2166AC")
  g = ggplot(data,aes(x=exp,y=mVals)) + geom_smooth(method=lm,se = F) + geom_point(data=data,aes(colour=status,shape=paired),size=4)
  g = g + scale_colour_manual(values=col) + scale_shape_manual(values=c("TRUE"=19,"FALSE"=1))
  g = g + theme_classic() + labs(x="Expression (logCPM)",y="Methylation (M-Value)",
                                 title=paste0("Cor: Exp vs. Methy \n(+/- 500bp ",element,")"))
  g = g + annotate("text",x=Inf,y=Inf,hjust=1,vjust=1,
                   label=paste0("r = ",signif(cor(data[,"exp"],data[,"mVals"]),2)),size=14)
  g = g + theme(legend.position = "none",
                #legend.justification = c(0,0),
                legend.title = element_blank(),
                #legend.key.size = unit(0.0085,"npc"),
                legend.text = element_text(size=24),
                aspect.ratio = 1,
                plot.margin = margin(t=8,r=0,b=50,l=0,unit="pt"),
                axis.title.y = element_text(size=32,colour="black"),
                axis.title.x = element_text(size=32,colour="black"),
                #axis.title.x = element_blank(),
                axis.text = element_text(size=26,colour="black",angle=0),
                #axis.text.x = element_blank(),
                plot.title = element_text(size=32,face="bold"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
  
  g
}

#' Create subfigures for Figure 3
#'
#' @param data 
#' @param element a string denoting a TE subfamily name 
#'
#' @return a ggplot2 object
#' @import ggplot2
#' @export
#'
#' @examples
plot_normal_tumor_beta = function(data,element){
  ordered_data = data[with(data,order(data$rescale_dis)),]
  bin = cut(ordered_data$rescale_dis,breaks=seq(-500,600,20))
  window = as.character(bin)
  window = sapply(window,function(x)strsplit(x,split=",",fixed=T)[[1]][1])
  window = substr(window, start=2,stop=nchar(window))
  ordered_data$window = as.numeric(window)
  ordered_melt = melt(ordered_data[,c("normal_beta","tumor_beta","window")],id.vars = c("window"))
  ordered_melt$variable = gsub("_beta","",ordered_melt$variable)
  tmp = ordered_melt[1,]
  tmp$value = mean(ordered_melt$value)
  
  col = c("tumor" = "#B2182B","normal" = "#2166AC")
  g = ggplot(data=ordered_melt,aes(x=window,y=value,colour=variable,fill=variable)) + geom_smooth(se =T,method="loess",size=0.75,alpha=0.15,fullrange=T) + theme_classic()
  g = g + geom_rect(data=tmp,aes(xmin=0,xmax=100,ymin=-Inf,ymax=Inf),fill="#DEEBF7",alpha=0.3,colour=NA)
  g = g + labs(x="Distance to TE",y="Beta Value",title=paste0("Methylation around\n",element," (Beta)"))
  g = g + ylim(0,1)
  g = g + scale_colour_manual(values = col) + scale_fill_manual(values = col)
  g = g + scale_x_continuous(breaks=c(-500,-200,-100,0,100,200,300,600),
                             labels=c("up_5kb","up_2kb","up_1kb","TE_start","TE_end","downn_1kb","down_2kb","down_5kb"))
  g = g + theme(legend.position = "none",
                #legend.justification = c(0,1),
                #legend.title = element_blank(),
                #legend.key.size = unit(0.0085,"npc"),
                #legend.text = element_text(size=18),
                aspect.ratio = 1,
                axis.title.y = element_text(size=32,colour="black"),
                axis.title.x = element_text(size=32,colour="black"),
                #axis.title.x = element_blank(),
                axis.text = element_text(size=26,colour="black",angle=0),
                axis.text.x = element_text(size=22,colour="black",angle=90,hjust=1,vjust=0.5),
                #axis.text.x = element_blank(),
                plot.title = element_text(size=32,face="bold"),
                panel.grid.major.x = element_line(linetype="dashed",colour="gray50",size=0.2),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank())
  g
}
