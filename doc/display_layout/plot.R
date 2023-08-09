# ---------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
# ---------------------------------------------------------------------------------------------------

library(ggplot2)
library(here)
library(patchwork)
library(RColorBrewer)
library(scales)

common_theme<-  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.8/.pt),
        axis.ticks = element_line(colour = "black", linewidth = 1/.pt),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black", face="bold"),
        plot.title = element_text(colour = "black", face="bold"),
        text=element_text(colour="black"),
        plot.margin = unit(c(0, 0, -6, 2), "pt"), # -6pt crop for bottom
        panel.spacing = unit(c(0, 0, 0, 0), "null"))
my.colors <- c("All k-mers.merged" = "#998ec3", "Shared k-mers.merged" = "#f1a340", "All k-mers.split" = "grey", "Shared k-mers.split" = "grey")

filenames <- list("mantis_40k.1000.tmax64.stats", "mantis_40k.1000.tmax64.norearrange.stats", "mantis_40k.1000.tmax192.stats", "RefSeq.k20.stats", "RefSeq.k20.prepare.stats")

for (filename in filenames)
{
  layout <- read.table(here(filename), header = T, row.names = 1, comment = '#')
  layout$index <- seq(1,nrow(layout))
  layout.t <- cbind(layout$index, layout$size, rep("All k-mers", nrow(layout)), layout$kind)
  layout.t <- rbind(layout.t, cbind(layout$index, layout$shared_size, rep("Shared k-mers", nrow(layout)), layout$kind))
  layout.t <- data.frame(layout.t)
  colnames(layout.t) <- c("index", "size", "type", "kind")
  layout.t$index <- as.numeric(layout.t$index)
  layout.t$size <- as.numeric(layout.t$size)
  layout.t$type <- as.factor(layout.t$type)
  layout.t$kind <- as.factor(layout.t$kind)

  g1 <- ggplot(data=layout.t, aes(x=index, y=size, fill=interaction(type,kind))) +
    geom_bar(stat="identity", position="identity") +
    common_theme +
    scale_fill_manual(values = my.colors, labels=c("All k-mers", "Shared k-mers", "Split bin"), breaks = c("All k-mers.merged", "Shared k-mers.merged", "All k-mers.split")) +
    ylab("Size of TBs\n") +
    labs (x = NULL) +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))

  g2 <- ggplot(data=layout, aes(x=index, y=ub_count)) +
    geom_bar(stat="identity", show.legend = FALSE) +
    common_theme +
    ylab("Amount of UBs\n") +
    labs (x = NULL)

  playout <- "
  AAA
  BBB
  "

  gg <- g1 + g2 + plot_layout(design = playout, guides = "collect") +
    plot_annotation(filename) &
    common_theme + theme(plot.title=element_text(hjust=0.5),
                         legend.position = "bottom",
                         legend.direction = "horizontal",
                         legend.title = element_blank(),
                         legend.text=element_text(colour="black"))

  pdf(here(paste(filename, ".pdf", sep="")), width = 7, height = 3.5)
  print(gg)
  invisible(dev.off())
}
