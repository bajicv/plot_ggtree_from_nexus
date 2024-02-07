################################################################################
# Script to plot ggtree from nexus file with mutations on branches
#
# Author:
#   Vladimir Bajić
#
# Date:
#   January 2024
#
# Usage:
#   Rscript --vanilla plot_ggtree.R -i input.nexus
#
#   TODO:
#   - add possibility to choose a root
#   - adjust colors
#
################################################################################

# Packages ---------------------------------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(ape))
suppressMessages(library(treeio))
library(optparse)

# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", metavar = "character",
        help = "Input file in nexus format with mutations [REQUIRED]\n"
    ),
    make_option(c("-o", "--output"),
        type = "character", metavar = "character",
        help = "Output path [OPTIONAL]\n"
    ),
    make_option(c("-r", "--root"),
        type = "character", metavar = "character",
        help = "ID that should be used to root the tree [OPTIONAL]\n"
    )
)

# Parsing options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check the provided option and execute the corresponding code -----------------

if (is.null(opt$i)) {
    print_help(opt_parser)
    stop("Input file must be provided.")
}

if (is.null(opt$o)) {
    opt$o <- tools::file_path_sans_ext(opt$i)
}

if (is.null(opt$r)) {
    cat("Using default root.\n")
}


# Load nexus -------------------------------------------------------------------
iqtree <- read.mega(opt$i)

# make mutation labels to be used in plot
iqtree@data$mutations_labels <- unlist(lapply(iqtree@data$mutations, paste, collapse = ", "))


# Plot tree --------------------------------------------------------------------

#  Adjustment for height and width of the plot to ensure that all the lables fit
plot_height <- 5 + (length(iqtree@phylo$tip.label) * 0.2)
plot_width <- 30 + (max(nchar(iqtree@phylo$tip.label)) * 0.3)

# Making tree
gg_tree <-
    iqtree %>%
    ggtree(layout = "rectangular") +
    geom_tiplab(align = TRUE, color = "darkblue") +
    coord_cartesian(clip = "off") +
    theme(plot.margin = unit(c(10, (plot_width * 3.5), 10, 5), "mm")) +
    geom_tippoint(color = "darkblue", size = 4) +
    geom_label(aes(x = branch, label = mutations_labels), fill = "lightgreen") +
    geom_rootedge()


# Save plot --------------------------------------------------------------------
ggsave(paste(opt$o, ".png", sep = ""), gg_tree, height = plot_height, width = plot_width, limitsize = FALSE)
ggsave(paste(opt$o, ".pdf", sep = ""), gg_tree, height = plot_height, width = plot_width, limitsize = FALSE)
