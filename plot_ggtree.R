################################################################################
# Script to plot ggtree from nexus file with mutations on branches
#
# Author:
#   Vladimir Bajić
#
# Last update:
#   16 February 2024
#
# Usage:
#   - Minimal example run requires a nexus file containing mutations
#   Rscript --vanilla plot_ggtree.R -i input.nexus
#
#   - Plot and output in specific location. Do not provide extension.
#     It automatically saves .jpeg .pdf and .png files.
#   Rscript --vanilla plot_ggtree.R -i input.nexus -o out_path_without_extension
#
#   - Plot with specific outgroup at the bottom of the tree.
#     sampleID's to be used for rooting should be within quotation marks("sampleID").
#   Rscript --vanilla plot_ggtree.R -i input.nexus -r "sampleID"
#
#   TODO:
#   - adjust colors
#
################################################################################


# Packages ---------------------------------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(ape))
suppressMessages(library(treeio))
library(optparse)


# Define functions -------------------------------------------------------------

## Function to plot ggtree -----------------------------------------------------
plot_ggtree <- function(tree) {
    gg_tree <-
        tree %>%
        ggtree(layout = "rectangular", ladderize = TRUE) +
        geom_tiplab(align = TRUE, color = "#081d58") +
        coord_cartesian(clip = "off") +
        theme(plot.margin = unit(c(10, (plot_width * 3.5), 10, 5), "mm")) +
        geom_tippoint(color = "#081d58", size = 4) +
        geom_label2(aes(x = branch, label = mutations_labels, subset = mutations_labels != ""), fill = "#a1d99b") +
        geom_rootedge()
    return(gg_tree)
}

## Function to plot rooted tree ------------------------------------------------
plot_ggtree_with_outgroup <- function(tree) {
    ### Choose outgroup
    tree@phylo <-
        tree@phylo %>%
        root(outgroup = opt$r, edgelabel = TRUE)

    ### Making temporary tree
    gg_tree_tmp <- ggtree(tree, layout = "rectangular", ladderize = FALSE)

    ### Reordering branches so that outgroup is at the very bottom
    init_taxa_order <- get_taxa_name(tree_view = gg_tree_tmp)
    outgroup_index <- grep(opt$r, init_taxa_order, fixed = TRUE)
    final_taxa_order <- append(init_taxa_order[-outgroup_index], init_taxa_order[outgroup_index])
    rev_final_taxa_order <- rev(final_taxa_order)

    ### rotate nodes so that outgroup is at the bottom
    tree@phylo <- ape::rotateConstr(tree@phylo, rev_final_taxa_order)

    ### Plot tree with the outgroup at the bottom
    gg_tree <-
        tree %>%
        ggtree(layout = "rectangular", ladderize = TRUE) +
        geom_tiplab(align = TRUE, color = "#081d58") +
        geom_tiplab(aes(subset = (label == opt$r)), align = TRUE, color = "#b30000") +
        coord_cartesian(clip = "off") +
        theme(plot.margin = unit(c(10, (plot_width * 3.5), 10, 5), "mm")) +
        geom_tippoint(color = "#081d58", size = 4) +
        geom_tippoint(aes(subset = (label == opt$r)), color = "#b30000", size = 4) +
        geom_label2(aes(x = branch, label = mutations_labels, subset = mutations_labels != ""), fill = "#a1d99b") +
        geom_point2(aes(subset = (label == opt$r)), shape = 23, size = 5, fill = "#b30000") +
        geom_rootedge()

    return(gg_tree)
}

## Function to adjust plot margins to ensure that all the labels fit ----------
adjust_plot_margin <- function(iqtree) {
    ### assign plot_height and plot_width to global environment
    plot_height <<- 5 + (length(iqtree@phylo$tip.label) * 0.2)
    plot_width <<- 30 + (max(nchar(iqtree@phylo$tip.label)) * 0.3)
}

## Function to save plots in different formats ---------------------------------
save_ggtree <- function(gg_tree) {
    ggsave(paste(opt$o, ".png", sep = ""), gg_tree, height = plot_height, width = plot_width, limitsize = FALSE)
    ggsave(paste(opt$o, ".pdf", sep = ""), gg_tree, height = plot_height, width = plot_width, limitsize = FALSE)
    ggsave(paste(opt$o, ".jpeg", sep = ""), gg_tree, height = plot_height, width = plot_width, limitsize = FALSE)
}


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

# Parsing options --------------------------------------------------------------
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check the provided option and execute the corresponding code -----------------

if (is.null(opt$i)) {
    print_help(opt_parser)
    stop("Input file must be provided.\n")
}

# Load nexus -------------------------------------------------------------------
cat("Reading in the tree.\n")
iqtree <- read.mega(opt$i)

# Make mutation labels to be used in plot --------------------------------------
cat("Making mutation labels.\n")
iqtree@data$mutations_labels <- unlist(lapply(iqtree@data$mutations, paste, collapse = ", "))

# Adjust plot margin to ensure that labels fit ----------------------------------
cat("Adjusting plot margins.\n")
adjust_plot_margin(iqtree)

if (is.null(opt$o)) {
    opt$o <- tools::file_path_sans_ext(opt$i)
    cat("Output not specified.\nOutput will be saved as:", opt$o, ".{jpeg,pdf,png}\n", sep = "")
}

if (is.null(opt$r)) {
    cat("Plotting unrooted tree.\n")
    ## Plot and save ggtree without specifying outgroup
    suppressWarnings({
        plot_ggtree(iqtree) %>% save_ggtree()
    })
} else {
    cat("Plotting tree with an outgroup:", opt$r, "\n")
    ## Plot and save ggtree without specifying outgroup
    suppressWarnings({
        plot_ggtree_with_outgroup(iqtree) %>% save_ggtree()
    })
}
