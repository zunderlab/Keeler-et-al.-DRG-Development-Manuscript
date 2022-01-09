#smgoggin
#07.23.20
#plotting (and prep for plotting) functions for protein/RNA comparison

library(ggrepel)
library(ggplot2)


##Line plot with dots, rna + protein separate plots, showing grayed out negative and colored for positive based on sort_labels
rna_SortedLabel_LineAndDot <- function(rna.in,rna.metadata.in,sub.list.common.names.for.sort, WRITE_PATH) {
  for (marker in sub.list.common.names.for.sort) {
    #Make DF with labeled colnames
    rna.in.df <- data.frame("SortLabel" = rna.metadata.in[,marker],
                            "Timepoint" = rna.metadata.in$timepoint_num,
                            "Expression" = rna.in[,marker])
    rna.in.df$SortLabel <- as.factor(rna.in.df$SortLabel)
    
    # Add in missing timepoints with NA and/or remove excess timepoints: #########################################################################
    total.time.points <- c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
    for (time.point in total.time.points) {
      if (time.point %in% unique(rna.in.df$Timepoint) == FALSE) {
        rna.in.df[nrow(rna.in.df)+1,] <- NA
        rna.in.df[nrow(rna.in.df),"Timepoint"] <- time.point
      }
    }
    ##Remove extra rna timepoint
    rna.in.df <- rna.in.df[-which(rna.in.df$Timepoint == 40),]
    
    # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
    # Generate the summary stats for the POSITIVE POPULATIONS ONLY!!!
    sumstat.rna.df <- summarySE(rna.in.df[which(rna.in.df$SortLabel == 1),], measurevar="Expression", groupvars=c("Timepoint","SortLabel"))
    
    #generate plot
    p_points <- ggplot(data=rna.in.df, aes(x=Timepoint, y=Expression, color=SortLabel)) + # data=rna.in.df, aes(x=Timepoint, y=Expression, color=Type) factor()  , group=Type
      geom_jitter(data=rna.in.df, aes(x=Timepoint, y=Expression, color=SortLabel), size=0.1) + #, alpha=0.05
      ###POINTS WITH LINE AND ERROR BARS FOR SUMMARY STATS:
      geom_line(data=sumstat.rna.df, aes(x=Timepoint, y=mean), position=position_dodge(0.2), color="black", size=2) +
      geom_point(data=sumstat.rna.df, aes(x=Timepoint, y=mean), position=position_dodge(0.2), color="black", size=4) +
      scale_x_continuous("Timepoint", breaks = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25),
                         labels=c("", "E11","E12","E13","E14","E15","E16","E17","E18","P0","P1","P2","P3","P4","P5", "")) +
      scale_color_manual(values = c("0"="grey", "1"="blue")) +
      coord_cartesian(xlim = c(10,25), #10.55,24.00
                      ylim = c(-0.2, max(sumstat.rna.df$mean + (sumstat.rna.df$sd))+2), #2*
                      expand = FALSE,default = FALSE,clip = "on") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.ticks = element_line(size = 1.5),
        legend.text = element_text(size = 7, face = "bold"),
        legend.title = element_text(size = 8, face = "bold"),
        axis.ticks.length=unit(.25, "cm"),
        axis.text.x = element_text(size = 7, face = "bold", angle = 45, vjust=-0.5),
        axis.text.y = element_text(size = 6, face = "bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.y = element_text(size = 9, face = "bold"),
        axis.title.x = element_text(size = 9, face = "bold"),
        plot.title = element_text(size=13, face = "bold", hjust = 0.5)) +
      ggtitle(marker) 
    
    ##SAVE PLOT!
    ggsave(paste0(WRITE_PATH, "figures/sorted_pop_timecourse_plots_pdfs/", "rna_sortLabels_", marker, "_timecourse_dot_plot.pdf"),  
           device = 'pdf',
           plot = p_points,
           width = 10,
           height = 5)
  }
}


##############################################################################################################################################################################################
##Line plot with dots, rna + protein separate plots, showing grayed out negative and colored for positive based on sort_labels
#' @import ggplot2

protein_SortedLabel_LineAndDot <- function(protein.in,protein.metadata.in,sub.list.common.names.for.sort, WRITE_PATH) {
  for (marker in sub.list.common.names.for.sort) {
    #Make DF with labeled colnames
    protein.in.df <- data.frame("SortLabel" = protein.metadata.in[,marker],
                                "Timepoint" = protein.metadata.in$timepoint_num,
                                "Expression" = protein.in[,marker])
    protein.in.df$SortLabel <- as.factor(protein.in.df$SortLabel)
    
    # Add in missing timepoints with NA and/or remove excess timepoints: #########################################################################
    total.time.points <- c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
    for (time.point in total.time.points) {
      if (time.point %in% unique(protein.in.df$Timepoint) == FALSE) {
        protein.in.df[nrow(protein.in.df)+1,] <- NA
        protein.in.df[nrow(protein.in.df),"Timepoint"] <- time.point
      }
    }
    
    # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
    # Generate the summary stats for the POSITIVE POPULATIONS ONLY!!!
    sumstat.protein.df <- summarySE(protein.in.df[which(protein.in.df$SortLabel == 1),], measurevar="Expression", groupvars=c("Timepoint","SortLabel")) 
    sumstat.protein.df$SortLabel <- as.factor(sumstat.protein.df$SortLabel)
    #generate plot
    p_points <- ggplot(data=protein.in.df, aes(x=Timepoint, y=Expression, color=SortLabel)) + 
      geom_jitter(data=protein.in.df, aes(x=Timepoint, y=Expression, color=SortLabel), size=0.1) + 
      ###POINTS WITH LINE AND ERROR BARS FOR SUMMARY STATS:
      geom_line(data=sumstat.protein.df, aes(x=Timepoint, y=mean), position=position_dodge(0.2), color="black", size = 2) + 
      geom_point(data=sumstat.protein.df, aes(x=Timepoint, y=mean), position=position_dodge(0.2), color="black", size = 4) + 
      scale_x_continuous("Timepoint", breaks = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25),
                         labels=c("", "E11","E12","E13","E14","E15","E16","E17","E18","P0","P1","P2","P3","P4","P5", "")) +
      scale_color_manual(values = c("0"="grey", "1"="orange")) +
      coord_cartesian(xlim = c(10,25), #10.55,24.00
                      ylim = c(-0.02, max(sumstat.protein.df$mean + (sumstat.protein.df$sd))+2), #2*
                      expand = FALSE,default = FALSE,clip = "on") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.ticks = element_line(size = 1.5),
        legend.text = element_text(size = 7, face = "bold"),
        legend.title = element_text(size = 8, face = "bold"),
        axis.ticks.length=unit(.25, "cm"),
        axis.text.x = element_text(size = 7, face = "bold", angle = 45, vjust=-0.5), 
        axis.text.y = element_text(size = 6, face = "bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.y = element_text(size = 9, face = "bold"),
        axis.title.x = element_text(size = 9, face = "bold"),
        plot.title = element_text(size=13, face = "bold", hjust = 0.5)) +
      ggtitle(marker) 
    
    ##SAVE PLOT!
    ggsave(paste0(WRITE_PATH, "figures/sorted_pop_timecourse_plots_pdfs/", "protein_sortLabels_", marker, "_timecourse_dot_plot.pdf"), 
           device = 'pdf',
           plot = p_points,
           width = 10,
           height = 5)
  }
}

##############################################################################################################################################################################################
##Line plot summarizinf over time expression of sorted populations, with violins showing distributions at each timepoint
rna_SortedLabel_Violin <- function(rna.in,rna.metadata.in, WRITE_PATH, plot_label, axis.side, violin.color, line.color) {
  #Make DF with labeled colnames #############################################################################################################
  rna.in.df <- data.frame(#"Type" = rep("RNA", length(rna.in)),
    "Timepoint" = rna.metadata.in$timepoint_num,
    "Expression" = rna.in)
  
  ##Remove extra rna timepoint: ###############################################################################################################
  tmpnt_40 <- which(rna.in.df$Timepoint == 40)
  if (length(tmpnt_40) >= 1) {
    rna.in.df <- rna.in.df[-tmpnt_40,]
  }
  
  # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
  # Generate the summary stats
  # Do this before adding missing timepoints, so that the line connects between non-adjacent timepoints
  sumstat.rna.df <- summarySE(rna.in.df, measurevar="Expression", groupvars=c("Timepoint"))
  
  # Add in missing timepoints with NA and/or remove excess timepoints: #########################################################################
  total.time.points <- c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
  for (time.point in total.time.points) {
    if (time.point %in% unique(rna.in.df$Timepoint) == FALSE) {
      rna.in.df[nrow(rna.in.df)+1,] <- NA
      rna.in.df[nrow(rna.in.df),"Timepoint"] <- time.point
    }
  }
  
  #generate plot for violins with discrete x-axis scaling: #########################################################################
  p_violin <- ggplot(data=rna.in.df, aes(x=as.factor(Timepoint), y=Expression)) +
    geom_violin(data=rna.in.df, aes(x=as.factor(Timepoint), y=Expression), scale = "width", fill=violin.color, color=violin.color, size = 0.25, alpha = 0.25) + 
    scale_x_discrete("Timepoint", breaks = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25),
                     labels=c("", "E11","E12","E13","E14","E15","E16","E17","E18","P0","P1","P2","P3","P4","P5", "")) +
    scale_y_continuous(position = axis.side) +
    coord_cartesian(ylim = c(0, max(sumstat.rna.df$mean + (sumstat.rna.df$sd))+2), expand = FALSE,default = FALSE,clip = "on") + 
    theme_minimal_hgrid(11, rel_small = 1) +
    theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      legend.position = "none") +
    ggtitle(plot_label) 
  
  #generate plot for line with continuous x-axis scaling: #########################################################################
  p_line <- ggplot(data=rna.in.df, aes(x=Timepoint, y=Expression)) + 
    geom_line(data=sumstat.rna.df, aes(x=Timepoint, y=mean), position=position_dodge(0.2), color=line.color, size = 1.5) +
    geom_point(data=sumstat.rna.df, aes(x=Timepoint, y=mean), position=position_dodge(0.2), color=line.color, size = 3.5) +
    scale_x_continuous("Timepoint", breaks = c(11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24),
                       labels=c("E11","E12","E13","E14","E15","E16","E17","E18","P0","P1","P2","P3","P4","P5")) +
    scale_y_continuous(position = axis.side) +
    coord_cartesian(xlim = c(10,25), 
                    ylim = c(0, max(sumstat.rna.df$mean + (sumstat.rna.df$sd))+2), #2*
                    expand = FALSE,default = FALSE,clip = "on") + 
    theme_half_open(11, rel_small = 1) +
    theme(legend.position = "none") +
    ggtitle(plot_label) 

  return(list(p_violin, p_line))
  
}

##############################################################################################################################################################################################
##Line plot with dots, rna + protein separate plots, showing grayed out negative and colored for positive based on sort_labels
protein_SortedLabel_Violin <- function(protein.in,protein.metadata.in, WRITE_PATH, plot_label, axis.side, violin.color, line.color) {
  #Make DF with labeled colnames #############################################################################################################
  protein.in.df <- data.frame("Timepoint" = protein.metadata.in$timepoint_num,"Expression" = protein.in)
  
  # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
  # Generate the summary stats
  # Do this before adding missing timepoints, so that the line connects between non-adjacent timepoints
  sumstat.protein.df <- summarySE(protein.in.df, measurevar="Expression", groupvars=c("Timepoint"))
  
  # Add in missing timepoints with NA and/or remove excess timepoints: #########################################################################
  total.time.points <- c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
  for (time.point in total.time.points) {
    if (time.point %in% unique(protein.in.df$Timepoint) == FALSE) {
      protein.in.df[nrow(protein.in.df)+1,] <- NA
      protein.in.df[nrow(protein.in.df),"Timepoint"] <- time.point
    }
  }
  
  #generate plot for violins with discrete x-axis scaling: #########################################################################
  p_violin <- ggplot(data=protein.in.df, aes(x=as.factor(Timepoint), y=Expression)) + 
    geom_violin(data=protein.in.df, aes(x=as.factor(Timepoint), y=Expression), scale = "width", fill=violin.color, color=violin.color, size = 0.25, alpha = 0.25) + 
    scale_x_discrete("Timepoint", breaks = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25),
                     labels=c("", "E11","E12","E13","E14","E15","E16","E17","E18","P0","P1","P2","P3","P4","P5", "")) +
    scale_y_continuous(position = axis.side) +
    coord_cartesian(ylim = c(0, max(sumstat.protein.df$mean + (sumstat.protein.df$sd))+2),expand = FALSE,default = FALSE,clip = "on") + 
    theme_minimal_hgrid(11, rel_small = 1) +
    theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      legend.position = "none") +
    ggtitle(plot_label) 
  
  #generate plot for line with continuous x-axis scaling: #########################################################################
  p_line <- ggplot(data=protein.in.df, aes(x=Timepoint, y=Expression)) + 
    geom_line(data=sumstat.protein.df, aes(x=Timepoint, y=mean), position=position_dodge(0.2), color=line.color, size = 1.5) +
    geom_point(data=sumstat.protein.df, aes(x=Timepoint, y=mean), position=position_dodge(0.2), color=line.color, size = 3.5) +
    scale_x_continuous("Timepoint", breaks = c(11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24),
                       labels=c("E11","E12","E13","E14","E15","E16","E17","E18","P0","P1","P2","P3","P4","P5")) +
    scale_y_continuous(position = axis.side) +
    coord_cartesian(xlim = c(10,25), #10.55,24.00
                    ylim = c(0, max(sumstat.protein.df$mean + (sumstat.protein.df$sd))+2), 
                    expand = FALSE,default = FALSE,clip = "on") + 
    theme_half_open(11, rel_small = 1) +
    theme(legend.position = "none") +
   ggtitle(plot_label) 

  return(list(p_violin, p_line))
  
}


##############################################################################################################################################################################################
#' @param sum.mat expression matrix summarized per feature (gene/protein) and consisting of gene levels, timepoint, data type
# generates circle plots in which expression level is color coded by dot size and percentage of expressing cells by dot color

plot.mean.exp.combined <- function(sum.mat, name){
  
  # plot using ggplot
  
  gg <- ggplot(data = as.data.frame(sum.mat), aes(x = factor(age), y = datasetID)) +
    geom_point(aes(size = percent, color = measurement, alpha = dataBinary )) + #, alpha = percent , shape = dataBinary datasetID fill = measurement, 
    scale_size_continuous(range = c(3,12)) +
    scale_alpha_manual(values = c("1" = 1, "0" = 0), guide = "none") +
    scale_color_viridis_c() +
    scale_y_discrete(position = "right") +
    theme_classic() +
    ylab("Dataset") +
    xlab("embryonic days") +
    facet_wrap(~genename, ncol = 1, strip.position = "left") +
    ggtitle(name) +
    labs(size="ratio expressing cells", color="normalized mean expression") +
    theme(panel.background = element_rect(fill = "grey92", colour = NA),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_blank(), #element_text(size = 11),
          panel.spacing.y = unit(0.5, "lines"),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_text(size = 11),
          legend.title = element_text(size = 8),
          legend.title.align=0.5,
          legend.text = element_text(size = 8),
          strip.text = element_text(size = 8, face = "bold"),
          strip.placement = "outside",
          strip.background = element_rect(fill = NULL,size = NULL,linetype = NULL,color = "white"), #"white"
          plot.title = element_text(face = "bold", size = 16))
  
  return(gg)
}

