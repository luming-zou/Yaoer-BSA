#!/usr/bin/env Rscript

# Prepare data

# Read in chromosome information from the reference genome
fai <- read.table("REFERENCE.fa.fai", header = FALSE, stringsAsFactors = FALSE)[1:31,]
starts <- c(1, cumsum(fai[,2])[1:30] + 1)
ends <- cumsum(fai[,2])
x <- (starts + ends) / 2
chr <- fai[,1]
chr_alias <- paste("chr", 1:31, sep = "")
fai <- data.frame(chr = chr, x = x, chr_alias = chr_alias, starts = starts, chr_nu = 1:31, stringsAsFactors = FALSE)
fai$color[fai$chr_nu %% 2 == 1] <- "color1"
fai$color[fai$chr_nu %% 2 == 0] <- "color2"

# Read in delta SNP index data
index <- read.table("delta_snp_index.csv", header = TRUE, sep = "\t")
colnames(index) <- c("chr", "start", "end", "value")
index$x <- (index$start + index$end) / 2

# Map chromosome information to delta SNP index data
for (i in 1:31){
  index$color[index$chr == fai$chr[i]] <- fai$color[i]
  index$chr_alias[index$chr == fai$chr[i]] <- fai$chr_alias[i]
  index$x[index$chr == fai$chr[i]] <- index$x[index$chr == fai$chr[i]] + fai$starts[i] - 1
  index$chr_nu[index$chr == fai$chr[i]] <- i
}
index <- index[!is.na(index$color),]

# Load ggplot2 library
library(ggplot2)

# Define threshold for highlighting
thread <- quantile(index$value, .99)

# Create a ggplot object
p <- ggplot(fai, aes(x = x, y = 0)) +
  geom_point(data = fai, aes(x = x, y = 0, colour = color), size = 0, show.legend = FALSE, stroke = 0) +
  geom_point(data = index, aes(x = x, y = value, colour = color), size = .6, show.legend = FALSE, stroke = 0) +
  geom_hline(yintercept = thread, color = 'red4', linetype = 2, linewidth = .4) +
  scale_x_continuous(breaks = fai$x, labels = fai$chr_nu, expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  scale_colour_manual(values = c(color1 = "#04544c", color2 = "#74bc54")) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 6, margin = margin(r = 0, l = 0)),
    axis.text.y = element_text(size = 10)
  ) +
  theme_bw(base_line_size = 0) +
  labs(x = "Chromosome", y = "Average delta index")

# Add smoothed lines for each chromosome
for (chr in fai$chr){
  gs <- index[index$chr == chr,]
  p <- p + geom_smooth(data = gs, aes(x = x, y = value), method = "gam", se = FALSE, linewidth = .5, colour = "red4")
}

# Save the plot as a JPEG file
ggsave("delta_index.jpg", plot = p, width = 22, height = 8, units = "cm", dpi = 800)
