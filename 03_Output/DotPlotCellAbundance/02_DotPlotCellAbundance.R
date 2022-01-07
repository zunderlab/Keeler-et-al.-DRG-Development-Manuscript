#Kristen Fread, University of Virginia
#20 Oct 2021
#code to make plots for figure 1 of paper with replicate data

rm(list = ls(all = TRUE))
#.libPaths("/project/zunderlab/R/3.6.3")

library(ZunderPipelineFunctions)
library(ggfortify)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(viridis)
library(cowplot)

## Input parameters
INPUT.FOLDER <- getwd()
FILES.FILENAME <- "filteredmatrix_file_numbers.csv" #to get this file from pipeline 1, run Save_file_nums_pipeline1.R to separate the filesnumbs from concattransformed
METADATA.FILENAME <- "metadata.csv"
CLUSTERS.FILENAME <- "/clusters.csv"
OUTPUT.FILENAME.Totalcells <- "total_cells_per_timepoint.png"
CLUSTER.NAMES.ORDER.COLOR <- "cluster_names_order_color.csv" #has to be individually crafted with labels and desired colors
OUTPUT.FILENAME.DOTPLOT.1 <- "dotplot_cell_abundance_per_timepoint1.png"
OUTPUT.FILENAME.DOTPLOT.2 <- "dotplot_cell_abundance_per_timepoint2.png"
OUTPUT.FILENAME.DOTPLOT.3 <- "dotplot_cell_abundance_per_timepoint3.png"


## Read needed files
fil_nums <- fread(paste0(INPUT.FOLDER, "/", FILES.FILENAME), stringsAsFactors = F)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)
# read metadata
md <- fread(paste0(INPUT.FOLDER, "/", METADATA.FILENAME), stringsAsFactors = F)
# assign file numbering
md <- cbind(md, "File" = 1:nrow(md)) 
#read this in for dot plot ordering
ord_name_color <- fread(paste0(INPUT.FOLDER, "/", CLUSTER.NAMES.ORDER.COLOR), stringsAsFactors = F)


## Prep dataframe for plotting
plotting_df <- as.data.frame(cbind(fil_nums,Cluster = clusters.in))

#join metadata to the dataframe with file and cluster number
cluster_file_metadata <- plotting_df %>% 
  left_join(md, by="File")

#calculate number of cells per age
totalcellstimepoint <- cluster_file_metadata %>% 
  group_by(age) %>% 
  summarise(num_cells=n())

total_cells <- sum(totalcellstimepoint$num_cells)

####### Plot the total cells per timepoint #######
p <- ggplot(totalcellstimepoint, aes(age, num_cells, fill = age)) + geom_col() + 
  scale_fill_viridis(discrete=TRUE, name = "Timepoint", guide=FALSE) +
  theme_bw() +
  ggtitle("Total cells per timepoint (532,926 cells total)") +
  labs(x = "Timepoint", y = "Number of cells") +
  theme(plot.title = element_text(size=30), axis.title = element_text(size = 30), 
        axis.text = element_text(size = 22), axis.text.x = element_text(angle = 70, hjust=1)) 

p

# save plot
ggsave(OUTPUT.FILENAME.Totalcells, 
       plot = p, width = 7, height = 9, units = c("in"), dpi = 600,
       device = "png"
)


#### Trying to figure out dot plot code ####

#calculate values to plot:
#group matrix by cluster and file then calculate how much of each file in the cluster
clusterz <- plotting_df %>%
  group_by(File, Cluster) %>% 
  summarize(num_cells=n())

#calculate the percent breakdown of each cluster contributing to a file/timepoint (on a per file basis, find the makeup of each file)
percent_similar <- clusterz %>%
  group_by(File) %>%
  summarize(total_cells_per_file = sum(num_cells)) %>%
  inner_join(clusterz) %>%
  mutate(fraction_cells_in_file = num_cells/total_cells_per_file) %>% 
  mutate(percent_cells_in_file = fraction_cells_in_file*100)

#calculate the normalized fraction of each timepoint within each cluster
cluster_plot_values <- percent_similar %>% 
  group_by(Cluster) %>% 
  summarize(total_percents_in_cluster = sum(percent_cells_in_file)) %>% 
  inner_join(percent_similar) %>% 
  mutate(normed_fraction_file_in_cluster = percent_cells_in_file/total_percents_in_cluster) %>% 
  mutate(normed_percent_file_in_cluster = normed_fraction_file_in_cluster*100)

#pull out matrix of just required info for plotting
normed_fractions <- cluster_plot_values %>% 
  select(Cluster, File, normed_fraction_file_in_cluster)

#add the day labels
days_fractions <- normed_fractions %>% 
  left_join(md[,c("File", "age")], by="File")

#sum each day together, this is specific for my plots since they have multiple replicates
sumeed_fractions <- ddply(days_fractions,.(Cluster,age),summarize,sum=sum(normed_fraction_file_in_cluster),number=length(age))

ordered_plotting <- sumeed_fractions %>% 
  left_join(ord_name_color, by="Cluster")

#make the dotplot
#inputs
dot.scale = 8
scale.min = 0
scale.max = 1
scale.by = 'size'

scale.func <- switch(
  EXPR = scale.by,
  'size' = scale_size,
  'radius' = scale_radius,
  stop("'scale.by' must be either 'size' or 'radius'")
)

#plot just the data without any modifications
plot <- ggplot(data = normed_fractions, mapping = aes_string(x = 'File', y = 'Cluster')) + #put in matrix containing x (day) and y (cluster)
  geom_point(mapping = aes_string(size = 'normed_fraction_file_in_cluster', color = 'Cluster')) +
  scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + #optional scaling method
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +  #sets x and y axis title to be blank
  guides(size = guide_legend(title = 'Percent Expressed')) +  #sets the legend title
  labs(x = 'File', y = 'Cluster') +
  theme_cowplot()

plot

# save plot
ggsave(OUTPUT.FILENAME.DOTPLOT.1, 
       plot = plot, width = 12, height = 12, units = c("in"), dpi = 300,
       device = "png")

#plotting with summed fractions of data per timepoint instead of individual file
plot2 <- ggplot(data = sumeed_fractions, mapping = aes_string(x = 'age', y = 'Cluster')) + #put in matrix containing x (day) and y (cluster)
  geom_point(mapping = aes_string(size = 'sum', color = 'Cluster')) +
  scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + #optional scaling method
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +  #sets x and y axis title to be blank
  guides(size = guide_legend(title = 'Percent Expressed')) +  #sets the legend title
  labs(x = 'Timepoint', y = 'Cluster') +
  theme_cowplot()
plot2

# save plot
ggsave(OUTPUT.FILENAME.DOTPLOT.2, 
       plot = plot2, width = 7, height = 9, units = c("in"), dpi = 600,
       device = "png")

#plotting trying to order by manually put in numbers
order_for_plotting <- ordered_plotting %>% 
  arrange(Plotting_Order)
order_for_plotting$Dot_Plot_Name <- factor(order_for_plotting$Dot_Plot_Name, levels=rev
                                           (unique(order_for_plotting$Dot_Plot_Name)[unique(order_for_plotting$Plotting_Order)]), ordered=TRUE)

coloring <- rev(unique(order_for_plotting$Color))

plot3 <- ggplot(data = order_for_plotting, mapping = aes_string(x = 'age', y = 'Dot_Plot_Name')) + #put in matrix containing x (day) and y (cluster)
  geom_point(mapping = aes_string(size = 'sum', color = 'Dot_Plot_Name')) +
  scale_color_manual(values = coloring) +
  scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + #optional scaling method
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +  #sets x and y axis title to be blank
  guides(size = guide_legend(title = 'Percent Expressed')) +  #sets the legend title
  labs(x = 'Timepoint', y = 'Cluster') +
  theme_cowplot()+
  theme(plot.title = element_text(size=30), axis.text.x = element_text(angle = 70, hjust=1))
plot3

# save plot
ggsave(OUTPUT.FILENAME.DOTPLOT.3, 
       plot = plot3, width = 7, height = 10, units = c("in"), dpi = 600,
       device = "png")