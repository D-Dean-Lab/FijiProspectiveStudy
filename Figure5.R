library(readxl)
library(tidyverse)
library(dendextend)
library(gplots)
library(vegan)
library(gridExtra)
library(grid)
library(gridGraphics)
library(dendextend)

setwd("~/path/") #point to the correct path

heatmap_grob <- function(df, title, left_side = TRUE) {
  print(title)
  gv_abundance <- df[, "Gardnerella_vaginalis"]
  li_abundance <- df[, "Lactobacillus_iners"]
  lc_abundance <- df[, "Lactobacillus_crispatus"]
  lactobacillus_abundance <- li_abundance + lc_abundance
  df_t <- as.data.frame(t(df))
  total_abundance <- apply(df_t, 1, sum)
  top_species <- order(total_abundance, decreasing = TRUE)[1:25]
  df_t <- df[, top_species]
  names(df_t) <- gsub("_", " ", names(df_t))
  italic_labels <- as.expression(lapply(colnames(df_t), function(x) bquote(italic(.(x)))))
  
  dist_mat_samples <- vegdist(df, method = 'bray')
  hclust_samples <- hclust(dist_mat_samples, method = 'complete')
  dendro_samples <- as.dendrogram(hclust_samples)
  
  dendro_samples_colored <- as.dendrogram(hclust_samples)
  
  
  if (title != "Baseline_pos_to_pos Heatmap" && title != "Followup_neg_to_neg Heatmap") {
    dendro_samples_colored <- reorder(dendro_samples_colored, wts = gv_abundance[order.dendrogram(dendro_samples_colored)], agglo.FUN = mean)
    dendro_samples_colored <- reorder(dendro_samples_colored, wts = lactobacillus_abundance[order.dendrogram(dendro_samples_colored)], agglo.FUN = mean)
    dendro_samples_colored <- rev(dendro_samples_colored)
  }
  if (title == "Baseline_pos_to_pos Heatmap") {
    dendro_samples_colored <- reorder(dendro_samples_colored, wts = gv_abundance[order.dendrogram(dendro_samples_colored)], agglo.FUN = mean)
    dendro_samples_colored <- rev(dendro_samples_colored)
  }
  if (title == "Followup_neg_to_neg Heatmap") {
    dendro_samples_colored <- reorder(dendro_samples_colored, wts = lactobacillus_abundance[order.dendrogram(dendro_samples_colored)], agglo.FUN = mean)
    dendro_samples_colored <- rev(dendro_samples_colored)
  }
  
  # Dynamically determine margins
  longest_col_label <- max(nchar(gsub("_", " ", names(df_t))))
  bottom_margin <- min(30, max(12, round(longest_col_label * 0.6)))
  row_label_lengths <- nchar(rownames(df))
  left_margin <- if (left_side) min(25, max(10, round(max(row_label_lengths) * 0.5))) else 5
  
  # Assign sample label colors based on suffix
  sample_names <- rownames(df)
  sample_colors <- sapply(sample_names, function(name) {
    if (grepl("C$", name)) {
      "red"
    } else if (grepl("V$", name)) {
      "red"
    } else if (grepl("R$", name)) {
      "blue"
    } else {
      "black"
    }
  })
  
  par(font = 2)
  
  heatmap <- heatmap.2(as.matrix(df_t),
                       dendrogram = "row",
                       Rowv = dendro_samples_colored,
                       Colv = "NA",
                       col = colorRampPalette(c("darkblue", "cyan", "green", "yellow", "orange", "red"))(100),
                       trace = "none",
                       margins = c(bottom_margin, left_margin),
                       cexRow = 1.4,
                       cexCol = 1.4,
                       srtCol = 90, adjCol = c(1, 1),
                       colRow = sample_colors,  # <-- apply color
                       keysize = 1,
                       density.info = "none",
                       key = FALSE,
                       labCol = italic_labels,
                       lhei = c(1, 30)  
  )
  
  
  grid.echo()
  grob <- grid.grab()
  return(grob)
}

# Load datasets and save each one individually
datasets <- list(
  "Baseline_pos_to_pos",
  "Followup_pos_to_pos",
  "Baseline_pos_to_neg",
  "Followup_pos_to_neg",
  "Baseline_neg_to_neg",
  "Followup_neg_to_neg"
)

for (i in seq_along(datasets)) {
  sheet_name <- datasets[[i]]
  df <- as.data.frame(read_xlsx("spreadsheet.xlsx", sheet = sheet_name)) #replace with the correct excel spreadsheet
  row.names(df) <- df$Sample_ID
  df <- df %>% select(-Dataset, -TimePoint, -Sample_ID, -`sum`)
  
  is_left <- i %% 2 == 1
  title_label <- gsub("_", " → ", sheet_name)
  heatmap <- heatmap_grob(df, paste(sheet_name, "Heatmap"), left_side = is_left)
  
 
  # Save to PNG
  png_filename <- paste0(sheet_name, "_heatmap.png")
  png(png_filename, width = 4000, height = 7500, res = 400) #adjust as needed
  grid.newpage()
  grid.draw(arrangeGrob(
    textGrob(title_label, gp = gpar(fontsize = 16, fontface = "bold")),
    heatmap,
    ncol = 1,
    heights = unit.c(unit(1.5, "lines"), unit(1, "null"))
  ))
  dev.off()
  
  print(paste("Saved:", png_filename))
}
