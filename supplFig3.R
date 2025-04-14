# Install/load circlize and RColorBrewer
if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}
library(circlize)

# Set working directory
setwd("PATH")

# Read data
data <- read.csv("Ct_persistence_rectal.csv") #change to the correct spreadsheet

# Define permanent, consistent CST color palette
CST_COLORS <- c(
  "ET-B"   = "#E41A1C",  # red
  "ET-B/P" = "#377EB8",  # blue
  "ET-F"   = "#4DAF4A",  # green
  "ET-P"   = "#6e515c",  # purple
  "I-A"    = "#FF7F00",  # orange
  "I-B"    = "#00FFFF",  # 
  "III-A"  = "#d703fc",  # 
  "III-B"  = "#CC0066",  # 
  "IV-A"   = "#66C2A5",  # teal
  "IV-B"   = "#808000",  # 
  "IV-C"   = "#8DA0CB",  # steel blue
  "IV-D0"  = "#8A2BE2",  # 
  "IV-D1"  = "#800000",  # maroon
  "IV-D2"  = "#12063d",  # yellow
  "IV-E"   = "#0000FF",  # blue
  "N"      = "#fcba03",  # 
  "V"      = "#D95F02"   # burnt orange
)

# Explicitly separate sectors for Baseline and Follow-up
data$Baseline_Sector <- paste0("B_", data$Baseline_CST)
data$Followup_Sector <- paste0("F_", data$Followup_CST)

# Sector order: baseline (left), followup (right)
sector_order <- c(
  paste0("B_", sort(unique(data$Baseline_CST))), 
  paste0("F_", sort(unique(data$Followup_CST)))
)

# Set gaps: big gap (15 degrees) between baseline and followup halves
gap_vector <- c(rep(2, length(unique(data$Baseline_CST)) - 1), 70,
                rep(2, length(unique(data$Followup_CST)) - 1), 70)

# Get CSTs in current dataset
cst_levels <- unique(c(data$Baseline_CST, data$Followup_CST))

# Subset permanent palette
cst_palette <- CST_COLORS[cst_levels]

# Match colors to sectors with prefixes
sector_colors <- setNames(
  cst_palette[gsub("^(B_|F_)", "", sector_order)],
  sector_order
)

# Set output file
png("Ct_persistence_rectal_plot.png", width = 2000, height = 2000, res = 300)

# Clear previous plot settings
circos.clear()

# Set plot parameters with correct gaps
circos.par(start.degree = -125, gap.after = gap_vector)

# Draw chord diagram
chordDiagram(
  data[, c("Baseline_Sector", "Followup_Sector")],
  grid.col = sector_colors,
  order = sector_order,
  annotationTrack = "grid",
  preAllocateTracks = 1,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  diffHeight = -0.001,
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.largest.ontop = TRUE
)

# Add sector labels neatly (without B_/F_ prefix for readability)
circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_label <- gsub("^(B_|F_)", "", CELL_META$sector.index)
    circos.text(
      CELL_META$xcenter, 
      CELL_META$ylim[1], 
      sector_label, 
      facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
      cex = 0.8, font = 2
    )
  },
  bg.border = NA
)

title("CST at Baseline/Followup: Ct persistence rectal samples")

# Save and close the PNG
dev.off()
