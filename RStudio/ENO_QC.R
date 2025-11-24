library(data.table)
library(Gviz)
library(RColorBrewer)
library(GenomicRanges)
library(grid)

# Set working directory
setwd("/project/sysviro/users/Max/analyses/QC")

# Genome/region to plot
myChr <- "ENO2"
myStart <- 1
myEnd <- 1314

# Bedgraph files
files <- list(
  "HSV2_333_ARPE19_WDX-3.sup-allMods.2h.trimAdapters.dorado.1.1.1.filtered.aligned.ENO2.primary.bedgraph",
  "HSV2_333_ARPE19_WDX-4.sup-allMods.4h_DMSO.trimAdapters.dorado.1.1.1.filtered.aligned.ENO2.primary.bedgraph",
  "HSV2_333_ARPE19_WDX-4.sup-allMods.4h_STM2457.trimAdapters.dorado.1.1.1.filtered.aligned.ENO2.primary.bedgraph",
  "HSV2_333_ARPE19_WDX-4.sup-allMods.6h.trimAdapters.dorado.1.1.1.filtered.aligned.ENO2.primary.bedgraph",
  "HSV2_333_ARPE19_WDX-4.sup-allMods.8h.trimAdapters.dorado.1.1.1.filtered.aligned.ENO2.primary.bedgraph",
  "HSV2_333_ARPE19_WDX-3.sup-allMods.CHX.trimAdapters.dorado.1.1.1.filtered.aligned.ENO2.primary.bedgraph"
)

# Labels for legend
labels <- c("2h", "4h DMSO", "4h STM2457", "6h", "8h", "CHX")

# Read, filter, and normalize data
tracks_data <- lapply(files, function(f) {
  dt <- fread(f, col.names = c("chromosome", "start", "end", "value"))
  dt <- dt[start > myStart & end < myEnd]
  dt$value <- dt$value / max(dt$value, na.rm = TRUE)  # Normalize to [0,1]
  return(dt)
})

# Colors for tracks
colors <- brewer.pal(n = length(tracks_data), "Blues")

# Ensure UCSC naming off
options(ucscChromosomeNames = FALSE)

# Create normalized DataTracks with shared ylim
tracks <- mapply(function(dt, col) {
  DataTrack(
    range = GRanges(seqnames = dt$chromosome,
                    ranges = IRanges(start = dt$start, end = dt$end),
                    score = dt$value),
    type = "l",
    chromosome = myChr,
    genome = "ENO2",
    col = col,
    ylim = c(0, 1),  # Shared normalized range
    col.axis = "black",
    background.title = "transparent",
    name = "Normalized Coverage"
  )
}, tracks_data, colors, SIMPLIFY = FALSE)

# Combine into overlay
overlay <- OverlayTrack(trackList = tracks, col.axis = "black", background.title = "transparent")

# Genome axis track
gtrack <- GenomeAxisTrack(col = "black")

# Output PDF (COVERAGE ONLY)
pdf("/project/sysviro/users/Max/analyses/Conditions/Results/QC/ENO2_normalized_overlay.pdf",
    width = 8, height = 5)

plotTracks(
  c(overlay, gtrack),
  from = myStart, to = myEnd,
  type = "l",
  col.histogram = NA,
  cex.title = 0,
  cex.axis = 1,
  title.width = 1,
  sizes = c(0.2, 0.05)
)

dev.off()


# Standalone legend as PNG
png("/project/sysviro/users/Max/analyses/Conditions/Results/QC/ENO2_normalized_overlay_legend.png",
    width = 800, height = 600, res = 300)

grid.newpage()

legend_x <- 0.1        # a bit in from left
legend_y <- 0.85       # near top
line_spacing <- 0.1   # vertical spacing between lines

for (i in seq_along(labels)) {
  y_pos <- legend_y - (i - 1) * line_spacing
  
  # line swatch
  grid.lines(
    x = unit(c(legend_x, legend_x + 0.1), "npc"),
    y = unit(rep(y_pos, 2), "npc"),
    gp = gpar(col = colors[i], lwd = 3)
  )
  
  # text label
  grid.text(
    labels[i],
    x = legend_x + 0.12,
    y = y_pos,
    just = "left",
    gp = gpar(fontsize = 11)
  )
}

dev.off()