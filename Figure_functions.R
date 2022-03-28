library(ggplot2)
library(viridis)
library(wesanderson)
library(ComplexHeatmap)
library(tidyr)
library(plotly)
library(gghighlight)
library(stringr)

cols = c(wes_palette("Darjeeling2")[2], "lightblue3", wes_palette("Darjeeling2")[3])

######################################################
#                                                    #
#  between treatment fold change line plots          #
#       for genes and exons                          #
######################################################

gene_exon_fc_plot = function(genes, title, highlight, contrast, use_direct_label){
  ids = read.delim("data/all_ids.txt", header = T, as.is = T)
  if(contrast == "control"){
    load("Cluster_analysis/2022_02_input_data/gene_and_exon_fc_between_CV_1.3fc.RData")
    dat = allfc[rownames(allfc) %in% genes,]
    d = as.data.frame(cbind(rep(0, nrow(dat)), dat))
    colnames(d) = c("t0", colnames(dat))
    dat = d
    ylab = "log2(SP+/V)"
    dat = as.data.frame(cbind(rownames(dat), dat))
    datLong = gather(dat, time, L2FC, t0:t24, factor_key = T)
    datLong = as.data.frame(cbind(datLong,
                                  rep(c(0, 0.5, 1, 2, 3, 4, 5, 6, 8, 12, 24), each = length(genes[genes %in% rownames(dat)])),
                                  rep(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12), each = length(genes[genes %in% rownames(dat)]))))
  }
  if(contrast == "mutant"){
      load("Cluster_analysis/2022_02_input_data/gene_and_exon_fc_between_SV.RData")
      dat = allfc[rownames(allfc) %in% genes,]
      d = as.data.frame(cbind(rep(0, nrow(dat)), dat))
      colnames(d) = c("t0", colnames(dat))
      dat = d
      ylab = "log2fc S- vs V"
      dat = as.data.frame(cbind(rownames(dat), dat))
      datLong = gather(dat, time, L2FC, t0:t24, factor_key = T)
      datLong = as.data.frame(cbind(datLong,
                                    rep(c(0, 0.5, 1, 2, 3, 4, 5, 6, 8, 12, 24), each = length(genes[genes %in% rownames(dat)])),
                                    rep(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12), each = length(genes[genes %in% rownames(dat)]))))
    }
  
  colnames(datLong) = c("gene", "timepoint", "log2FC",  "time", "breaks") 
  
  getSymbols = unlist(sapply(datLong[,1], function(x) strsplit(x, split = ":")[[1]][1]))
  exonInfo = unlist(sapply(datLong[,1], function(x) strsplit(x, split = ":")[[1]][2]))
  idx = match(getSymbols, ids[,2])
  datLong$symbol = paste(ids[idx, 1], exonInfo, sep = ":")
  h = datLong[,6] %in% highlight
  datLong = as.data.frame(cbind(datLong, h))
  
  g = 
    ggplot(datLong, aes(x = as.numeric(as.character(breaks)), 
                          y = as.numeric(as.character(log2FC)), group = symbol, color = symbol)) + 
    geom_line(aes(group = symbol), size=1.5, alpha = 0.5) + 
    geom_hline(yintercept = 0)+
    gghighlight(h, use_direct_label = use_direct_label,label_params = list(size = 10), unhighlighted_params = list(size = 0.75)) + # label_params = list(size = 8),
    geom_point(aes(group = symbol), size=0.75, alpha = 0.5) + 
    theme_classic(base_line_size = 1.2,
                  base_rect_size = 1.2)+
    theme(plot.title=element_text(size=24),
          axis.text=element_text(size=24, colour = "black"),
          legend.text=element_text(size=24),
          legend.title=element_text(size=24),
          axis.title.x=element_text(size=24, margin=margin(15,0,0,0)),
          axis.title.y=element_text(size=24, margin=margin(0,5,0,0))) +
    scale_x_continuous("hours", breaks=as.numeric((datLong$breaks)), 
                       labels = (datLong$time))+
    ylab(ylab)+
    ggtitle(title)
  plot(g)
}

