#!/mnt/sequence/mhu2/miniconda2/envs/r-4.0/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  stop("Please provide: newick_file, cdd_query, cdd_domain, output_pdf, sequence_lengths", call.=FALSE)
} 

library(tidyverse)
library(ggplot2)
library(ggtree)
library(treeio)
library(aplot)
library(gridExtra)

tr <- read.newick(args[1])


fasta_len <- read.table(args[5], header = F, sep = "\t", stringsAsFactors = F)
colnames(fasta_len) <- c("entryName", "end")

gene_len <- fasta_len
gene_len$begin <- 1
gene_len$length <- gene_len$end - gene_len$begin
gene_len$type <- "CHAIN"
gene_len$description <- gene_len$entryName

gene_len <- gene_len[, c("description", "begin", "end", "entryName", "type", "length")]


query_data <- read.table(args[2], header = F, sep = "\t", fill = TRUE)
domain_data <- tryCatch({
  read.table(args[3], header = F, sep = "\t", fill = TRUE)
}, error = function(e) data.frame())

if(nrow(domain_data) > 0) {

  domain_merged <- left_join(domain_data[, c(2,5,6,10)], query_data[, c(2,5)], by = c("V2" = "V2"))
  

  data_domain <- domain_merged[, c(4, 2, 3, 5)]
  colnames(data_domain) <- c("description", "begin", "end", "entryName")
  data_domain$type <- "DOMAIN"
  data_domain$length <- data_domain$end - data_domain$begin
  

  data_domain <- data_domain[!is.na(data_domain$entryName) & !is.na(data_domain$begin), ]
  
  
  data <- rbind(data_domain, gene_len)
} else {
  
  data <- gene_len
}

data$species   <- sub("^(.*?)_(.*)$", "\\1", data$entryName)
data$accession <- sub("^(.*?)_(.*)$", "\\2", data$entryName)

tree_order <- fortify(tr) %>% subset(isTip) %>% arrange(desc(y))

# reorder gene according to the order of phylogenetic tree
data$entryName <- factor(data$entryName, levels = rev(tree_order$label))
data$order <- as.integer(data$entryName)


data <- data[!is.na(data$order), ]


p <- ggplot() +
  
  geom_segment(data = subset(data, type == "CHAIN"),
               aes(x = begin, xend = end, y = order, yend = order),
               color = "darkgrey", linewidth = 1)


if("DOMAIN" %in% data$type) {
  p <- p + geom_rect(data = subset(data, type == "DOMAIN"),
                     aes(xmin = begin, xmax = end, ymin = order - 0.15, ymax = order + 0.15, fill = description))
}

p <- p + 
  scale_y_continuous(limits = c(0.5, max(data$order, na.rm=T) + 0.5), expand = c(0, 0)) +
  labs(fill = "Motif / Domain") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "right")


p_tree <- ggtree(tr, branch.length = "none", linewidth = 1)


p_tree$data$label <- sapply(p_tree$data$label, function(x) {
  if (is.na(x)) return(x) 
  if (grepl("^Uniprot\\|", x)) {
    return(x) 
  } else if(grepl("^Ref\\|", x)){
    return(x) 
  } else if (grepl("\\|", x)) {
    
    parts <- strsplit(x, "\\|")[[1]]
    id_part <- parts[1]
    gene_name <- parts[2]
    species_code <- sub("_.*", "", id_part)
    return(paste0(species_code, "|", gene_name))
  }
  return(x)
})
# ===============================================

p_tree <- p_tree + 
  geom_tiplab(align = TRUE, size = 5, offset = 0.6) +
  scale_x_continuous(expand = expansion(mult = c(0, 3))) 

pdf(args[4], height = 10, width = 15)
grid.arrange(p_tree, p, ncol = 2, widths = c(1, 2.25))
dev.off()