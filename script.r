# Set the working directory to save the final results from this script. It is different where structure results are located 
setwd("C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Fruit Morphology Projects/GWAS/structure/STR/RUN4")

# Install necessary packages for analysis and visualization
install.packages(c("devtools","ggplot2","gridExtra","gtable","label.switching","tidyr"), dependencies = TRUE)

# Load required libraries for the analysis
library(devtools)  # For installing and managing development tools
library(gridExtra)  # For arranging multiple grid-based plots
library(ggplot2)    # For creating visualizations

install.packages(c("rlang", "vctrs", "purrr"))

# install pophelper package from GitHub
remotes::install_github('royfrancis/pophelper')
library(pophelper)  # For working with STRUCTURE output files


# Define custom color palettes for plotting STRUCTURE results
clist <- list(
  "shiny" = c("#1D72F5", "#DF0101", "#77CE61", "#FF9326", "#A945FF", "#0089B2", "#FDF060", "#FFA6B2", "#BFF217", "#60D5FD", "#CC1577", "#F2B950", "#7FB21D", "#EC496F", "#326397", "#B26314", "#027368", "#A4A4A4", "#610B5E"),
  "strong" = c("#11A4C8", "#63C2C5", "#1D4F9F", "#0C516D", "#2A2771", "#396D35", "#80C342", "#725DA8", "#B62025", "#ED2224", "#ED1943", "#ED3995", "#7E277C", "#F7EC16", "#F8941E", "#8C2A1C", "#808080"),
  "oceanfive" = c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
  "keeled" = c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28"),
  "vintage" = c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
  "muted" = c("#46BDDD", "#82DDCE", "#F5F06A", "#F5CC6A", "#F57E6A"),
  "teal" = c("#CFF09E", "#A8DBA8", "#79BD9A", "#3B8686", "#0B486B"),
  "merry" = c("#5BC0EB", "#FDE74C", "#9BC53D", "#E55934", "#FA7921"),
  "funky" = c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762", "#ED8F47", "#9471B4"),
  "retro" = c("#01948E", "#A9C4E2", "#E23560", "#01A7B3", "#FDA963", "#323665", "#EC687D"),
  "cb_paired" = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
  "cb_set3" = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"),
  "morris" = c("#4D94CC", "#34648A", "#8B658A", "#9ACD32", "#CC95CC", "#9ACD32", "#8B3A39", "#CD6601", "#CC5C5B", "#8A4500"),
  "wong" = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#006699", "#D55E00", "#CC79A7"),
  "krzywinski" = c("#006E82", "#8214A0", "#005AC8", "#00A0FA", "#FA78FA", "#14D2DC", "#AA0A3C", "#FA7850", "#0AB45A", "#F0F032", "#A0FA82", "#FAE6BE")
)

# Calculate the length of each palette and update their names
lengths <- sapply(clist, length)
names(clist) <- paste0(names(clist), "_", lengths)

# Set graphical parameters for displaying palettes
par(mar = c(0.2, 6, 0.2, 0))  # Margins
par(mfrow = c(length(clist), 1))  # Arrange plots in a single column

# Loop through the palettes and create bar plots for visualization
for (i in 1:length(clist)) {
  barplot(
    rep(1, max(lengths)),
    col = c(clist[[i]], rep("white", max(lengths) - length(clist[[i]]))),
    axes = FALSE,
    border = FALSE
  )
  text(x = -0.1, y = 0.5, adj = 1, label = names(clist)[i], xpd = TRUE, cex = 1.2)
}

# Load STRUCTURE output files from the specified directory
sfiles <- list.files(
  path = "C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Fruit Morphology Projects/GWAS/structure/STR/RUN4/RUN4/Results", full.names = TRUE)

# Read STRUCTURE output files into a list
slist <- readQ(files = sfiles)
readQ(files = sfiles, filetype = "structure")  # Specify file type explicitly

# Summarize and tabulate the STRUCTURE results
sr1 <- summariseQ(tabulateQ(slist))

# Export Evanno method summary results to a CSV
write.csv(evannoMethodStructure(sr1), "evannoMethodStructure.csv", na= "")

# Generate and display the Evanno method plot for STRUCTURE analysis
p <- evannoMethodStructure(
  data = sr1,
  exportplot = FALSE,
  returnplot = TRUE,
  returndata = FALSE,
  basesize = 12,
  linesize = 0.7
)
grid.arrange(p)  # Arrange and display the plot

clusters <- evannoMethodStructure(data = sr1, exportplot = F, returnplot = T, returndata = T)
str(clusters)

# Extract the gtable and convert it to a grob to save the plot 
plot_grob <- grid::grobTree(clusters$plot)

# Save the grob as an image
ggsave(
  filename = "Cluster_K_fixed.png",
  plot = plot_grob,
  width = 10,
  height = 6,
  dpi = 300
)

#  Check the number of clusters (columns) in each entry of 'slist'
# 'slist' stores the STRUCTURE outputs, and each entry corresponds to a different K.
sapply(slist, ncol)  
# Load required packages
library(pophelper)  # For STRUCTURE plot visualization

# NOTE: Before running this script, ensure that the optimal number of clusters (K)
# is determined using the Evanno method. The chosen K should correspond to the highest Delta K value.


# Compute TRUE K for each run
run_K_table <- data.frame(Run = names(slist),K = sapply(slist, ncol),
  stringsAsFactors = FALSE)

# View mapping
print(run_K_table)

# Decide which K to plot - based on evanno method 
K_selected <- 3   

# Find runs that correspond to selected K
selected_runs <- run_K_table$Run[run_K_table$K == K_selected]
print(selected_runs)

if (length(selected_runs) == 0) {
  stop(paste("No STRUCTURE runs found for K =", K_selected))
}

# Use first run for that K (or average later if needed)
selected_run_name <- selected_runs[1]
selected_index <- which(names(slist) == selected_run_name)

cat("Selected run:", selected_run_name, " | K =", K_selected, "\n")

# extract Q-matrix
qmatrix <- slist[[selected_index]]

stopifnot(ncol(qmatrix) == K_selected)

# Assign individuals to dominant cluster
# For each individual, it finds the cluster with the highest ancestry proportion
max_cluster <- apply(qmatrix, 1, which.max)

# Sort individuals by cluster
sorted_indices <- order(max_cluster)
qmatrix_sorted <- qmatrix[sorted_indices, ]

# Convert to qlist 
qlist_sorted <- as.qlist(
  setNames(list(qmatrix_sorted), paste0("K", K_selected))
)

# Plot structure 
p1_sorted <- plotQ(
  qlist_sorted,
  returnplot = TRUE,
  exportplot = FALSE,
  basesize = 11,
  clustercol = clist$shiny,
  showlegend = TRUE,
  indlabsize = 14,        # genotype label size
  legendkeysize = 20,
  legendtextsize = 20,
  showindlab = TRUE,      # <<< SHOW NAMES
  splabsize = 19,
  splab = paste0("K = ", K_selected)
)


# Display plot
p1_sorted$plot[[1]]

# save the plot 
ggsave(
  filename = paste0("STRUCTURE_K", K_selected, "_Sorted.png"),
  plot = p1_sorted$plot[[1]],
  width = 18,
  height = 7,
  dpi = 300
)


library(writexl)   # if not installed: install.packages("writexl")
# Genotype names (must exist as rownames)
geno_names <- rownames(qmatrix)

# Dominant cluster assignment
dominant_cluster <- apply(qmatrix, 1, which.max)

# Build dataframe
cluster_df <- data.frame(
  Genotype = geno_names,
  Cluster  = dominant_cluster,
  stringsAsFactors = FALSE
)

# Split by cluster
cluster_list <- split(cluster_df, cluster_df$Cluster)

# Export to Excel (one sheet per cluster)
write_xlsx(cluster_list, path = "Genotypes_By_Cluster.xlsx")

# NOTE:
# STRUCTURE only remembers the order of genotypes, not their names.
# So row 1 in the results = the first genotype in your input file (e.g., GBS001),
# row 2 = the second genotype (GBS002), and so on
# To make it clear which genotype belongs to which cluster, always attach the
# real names from your original file. This keeps everything correct and easy to read



