library(shiny)
library(ggplot2)
library(dplyr)
library(plotgardener)
library(grid)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(tidyverse)

manhattan_plot <- function(input_file, plotting_file, enhancer_file) {
  # Read input data
  input_df <- read.csv(input_file)
  input_ordered <- input_df[, c("chrom", "pos", "p", "snp", "Sig.region")]

  # Read plotting data
  Sig_plotting_df <- read.csv(plotting_file, header = TRUE)
  Sig.regions <- unique(Sig_plotting_df$Sig.region)

  Sig.region.toplot <- Sig_plotting_df[which(Sig_plotting_df$Sig.region == Sig.regions[1]),]
  input <- input_ordered[which(input_ordered$Sig.region == Sig.regions[1]),]

  # Promoter information
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  promoters_txdb <- promoters(genes(txdb, single.strand.genes.only=FALSE))
  promoter_df <- unlist(promoters_txdb)

  # Enhancer information
  enhancer_df <- read.csv(enhancer_file)
  remove(promoters_txdb)

  # Create plotgardener page
  pageCreate(width = 20, height = 8, default.units = "inches", showGuides = FALSE)

  # Plot GWAS data
  manhattanPlot <- plotManhattan(
    data = input_ordered,
    chrom =  input_ordered$chrom[1],
    chromstart = Sig.region.toplot$CHR_llim[1],
    chromend = Sig.region.toplot$CHR_ulim[1],
    assembly = "hg19",
    ymax = 1.1,
    cex = 0.5,
    fill = "black",
    sigLine = TRUE,
    lty = 2,
    range = c(0, max((-log10(input$p) + 1))),
    x = 4,
    y = 0.5,
    width = 12,
    height = 3,
    just = c("left", "top"),
    default.units = "inches"
  )

  # Annotate the plot
  annoGenomeLabel(plot = manhattanPlot, x = 4, y = 3.5, fontsize = 15, just = c("left", "top"), default.units = "inches")
  plotText(label = "Chromosome", fontsize = 15, x = 10, y = 3.85, just = "center", default.units = "inches")
  annoYaxis(plot = manhattanPlot, at = seq(0, max(-log10(input$p), na.rm = TRUE), 2), axisLine = TRUE, fontsize = 15)
  plotText(label = "-log10(p-value)", x = 0.1, y = 1.75, rot = 90, fontsize = 15, fontface = "bold", just = "center", default.units = "inches")

  # Plot gene track
  Genes2plot <- as.data.frame(Sig_plotting_df[, c(1, 10)])
  colnames(Genes2plot) <- c("gene", "color")
  plotGenes(
    chrom = input_ordered$chrom[1],
    chromstart = Sig.region.toplot$CHR_llim[1],
    chromend = Sig.region.toplot$CHR_ulim[1],
    assembly = "hg19",
    x = 4,
    y = 4,
    width = 12,
    height = 1,
    just = c("left", "top"),
    default.units = "inches",
    geneHighlights = Genes2plot,
    geneBackground = '#666666'
  )

  # Plot promoters
  plotRanges(
    promoter_df,
    chrom = "chr22",
    assembly = "hg19",
    x = -2,
    y = 4,
    height = 2,
    width = 18,
    collapse = FALSE,
    boxHeight = unit(0.25, "inches"),
    fill = "#41ab5d",
    spaceWidth = 0,
    spaceHeight = 0,
    default.units = "inches"
  )

  # Plot enhancers
  plotRanges(
    enhancer_df,
    chrom = "chr22",
    assembly = "hg19",
    x = -2,
    y = 4.75,
    height = 2,
    width = 18,
    collapse = FALSE,
    boxHeight = unit(0.25, "inches"),
    fill = "#00adff",
    spaceWidth = 0,
    spaceHeight = 0,
    default.units = "inches"
  )
}
