#R scripts for generating figures and tables in the manuscript.
##Fig 1A. Species composition
library(ggplot2)
library(ggalluvial)
library(ggsci)
dat <- read.delim("protists.rb.txt", header = T)
dat$taxon <- factor(dat$taxon, levels = rev(unique(dat$taxon)))
dat$Group <- factor(dat$Group, levels = unique(dat$Group))
## plot figure
p <- ggplot(data = dat, aes(x = Group, y = value, alluvium = taxon, stratum = taxon))
p1 <- p + geom_alluvium(aes(fill = taxon), alpha = .5, width = 0.6) +
  geom_stratum(aes(fill = taxon), width = 0.6, alpha = .8)
p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p3 <- p2 + scale_fill_manual(values = rev(c(
  "#046586", "#28A9A1", "#C9A77C",
  "#F4A016", "#F6BBC6", "#E71F19", "gray57"
)))
p4 <- p3 + theme_classic()
p5 <- p4 + theme(axis.text.x = element_text(size = 16, family = "serif", angle = 45, hjust = 1, colour = "black")) +
  theme(axis.text.y = element_text(size = 16, family = "serif", colour = "black")) +
  theme(axis.title.x = element_text(size = 16, family = "serif", colour = "black")) +
  theme(axis.title.y = element_text(size = 16, family = "serif", colour = "black")) +
  theme(legend.text = element_text(size = 15, family = "serif", colour = "black")) +
  theme(legend.title = element_text(size = 12, hjust = 0, family = "serif", colour = "black")) +
  theme(strip.text = element_text(size = 14, family = "serif"))
p6 <- p5 + scale_y_continuous(expand = c(0, 0))
p7 <- p6 + theme(legend.text = element_text(colour = "black", size = 15)) +
  theme(legend.title = element_text(size = 16, colour = "black"))
p8 <- p7 + theme(text = element_text(family = "Times"))
p8
##Fig 1D. Principle coordinate analysis (PCoA)
library(ape)
library(ggplot2)
library(grid)
library(ggsci)
library(ggrepel)
otu <- read.table(file = "otu.table.protists.txt", header = T, row.names = 1, sep = "\t")
otu <- log(otu + 1)
groups <- read.delim("group.txt", header = T, row.names = 1)
samp.fg <- colnames(otu)
samp.env <- rownames(groups)
my.env <- match(samp.fg, samp.env)
group <- na.omit(groups[my.env, ]) # omit the NA rows if without fg data
samp.env <- rownames(group)
my.fg <- match(samp.env, samp.fg)
otu2 <- otu[, my.fg]
group$Depth <- factor(group$Depth, levels = c("Photic", "Aphotic"))
library(vegan)
library(maptools)
otu.dist <- vegdist(t(otu2), scale = T, method = "bray", na.rm = T)
# Pcoa
PCOA <- pcoa(otu.dist, correction = "none", rn = NULL) # ????PCOA()<U+05B8>????pcoa????
result <- PCOA$values[, "Relative_eig"]
pro1 <- as.numeric(sprintf("%.3f", result[1])) * 100
pro2 <- as.numeric(sprintf("%.3f", result[2])) * 100
x <- PCOA$vectors
sample_names <- rownames(x)
pc <- as.data.frame(PCOA$vectors)
pc$names <- sample_names
legend_title <- ""
pc$Axis.1 <- pc$Axis.1
pc$Axis.2 <- pc$Axis.2
pc$Depth <- group$Depth
xlab <- paste("PCOA1(", pro1, "%)", sep = "")
ylab <- paste("PCOA2(", pro2, "%)", sep = "")
# ggplot
p <- ggplot(pc, aes(Axis.1, Axis.2)) +
  geom_point(aes(color = Depth, shape = Depth), size = 4) +
  labs(x = xlab, y = ylab, title = "PCOA", color = legend_title) +
  geom_hline(yintercept = 0, linetype = 4, color = "grey") +
  geom_vline(xintercept = 0, linetype = 4, color = "grey") +
  theme_bw()
p <- p + theme(axis.text.x = element_text(size = 16, family = "serif", color = "black")) +
  theme(axis.text.y = element_text(size = 16, family = "serif", color = "black")) +
  theme(axis.title.x = element_text(size = 16, family = "serif", color = "black")) +
  theme(axis.title.y = element_text(size = 16, family = "serif", color = "black")) +
  theme(legend.text = element_text(size = 14, family = "serif", color = "black")) +
  theme(legend.title = element_text(size = 14, hjust = 0, family = "serif", color = "black")) +
  theme(strip.text = element_text(size = 16, family = "serif", color = "black"))
p <- p + scale_color_manual(values = c("#CD5555", "#0000AA"))
p

## Fig 2 variation partitioning analysis
R codes for variation partitioning analysis modified according to Wu et al. (2018)
##Fig 3A. Fitting Sloan’s neutral model
library(Hmisc)
library(minpack.lm)
library(stats4)
# spp: 	A community data matrix with each row for each species, and each column for each sample.
spp <- read.csv("otu.table.txt", head = T, stringsAsFactors = F, row.names = 1, sep = "\t")
spp <- t(spp)
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m / N
spp.bi <- 1 * (spp > 0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by = 0)
C <- C[order(C[, 2]), ]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]
p <- C.0[, 2]
freq <- C.0[, 3]
names(p) <- C.0[, 1]
names(freq) <- C.0[, 1]
d <- 1 / N
m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))
m.fit
m.ci <- confint(m.fit, "m", level = 0.95)
freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
Rsqr
pdf("NCM.pdf", width = 7, height = 6)
bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[, 2:3])
inter.col <- rep("black", nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- "#A52A2A"
inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- "#29A6A6"
library(grid)
grid.newpage()
pushViewport(viewport(h = 0.6, w = 0.6))
pushViewport(dataViewport(xData = range(log10(bacnlsALL$p)), yData = c(0, 1.02), extension = c(0.02, 0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, gp = gpar(col = inter.col, cex = 0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp = gpar(col = "blue", lwd = 2), default = "native")

grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")
grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")
grid.text(y = unit(0, "npc") - unit(2.5, "lines"), label = "Mean Relative Abundance (log10)", gp = gpar(fontface = 2))
grid.text(x = unit(0, "npc") - unit(3, "lines"), label = "Frequency", gp = gpar(fontface = 2), rot = 90)
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=", round(Rsqr, 3), "\n", "m=", round(coef(m.fit), 3)), x = x[j], y = y[i], just = just)
}
x <- unit(1:4 / 5, "npc")
y <- unit(1:4 / 5, "npc")
draw.text(c("centre", "bottom"), 4, 1)
dev.off()
##Fig 3G. The C-score metric
library(EcoSimR) # load EcoSimR library
library(devEMF)
set.seed(1234) # for reproducible results

T_list <- c("otu.table.protists.txt", "otu.table.FL.txt", "otu.table.FL.txt")

for (n in T_list) {
  OTU_Abu <- read.table(n, header = T)
  OTU_Abu[OTU_Abu > 0] <- 1
  
  # Filter out empty rows
  OTU_Abu.nonzerorow <- OTU_Abu[which(rowSums(OTU_Abu) > 0), ]
  OTU_Abu <- OTU_Abu.nonzerorow
  
  csModel <- cooc_null_model(OTU_Abu,
                             algo = "sim9", metric = "c_score",
                             nReps = 30000, saveSeed = FALSE, burn_in = 500, algoOpts = list(),
                             metricOpts = list(), suppressProg = FALSE
  )
  n
  summary(csModel)
  
  write.table(n, "c-score.N1w1.xls", append = TRUE)
  sink("c-score.N1w1.xls", append = TRUE)
  summary(csModel)
  sink(NULL)
  
  emf(
    file = sprintf("%s.c-score.hist.30000.emf", n), width = 7, height = 7,
    bg = "transparent", fg = "black", pointsize = 12,
    family = "Helvetica", custom.lty = FALSE
  )
  plot(csModel, type = "hist")
  dev.off()
  
  emf(
    file = sprintf("%s.c-score.burnin.30000.emf", n), width = 2.1, height = 2.1,
    bg = "transparent", fg = "black", pointsize = 12,
    family = "Helvetica", custom.lty = FALSE
  )
  plot(csModel, type = "burn_in")
  dev.off()
}
##Fig 4. (A-B) The Co-occurrence networks were constructed by sparcc algorithm.
##Fig 4C. 
library(ggsci)
library(ggplot2)
library(ggthemes)
library(ggpmisc)
dis <- read.csv("subnet-protists.csv", header = T)
dis$Group <- factor(dis$Group, levels = c("Total", "Photic", "Aphotic"))
names(dis)
p <- ggplot(data = dis, aes(x = Richness.log, y = Node.log, color = Group, shape = Group, group = Group))
p <- p + geom_point(size = 4, alpha = 0.5) + theme_few()
p <- p + labs(x = "OTU richness (log-transformed)", y = "Node numbers (log-transformed)", title = "")
p <- p + theme(axis.text.x = element_text(size = 16, family = "serif", color = "black")) +
  theme(axis.text.y = element_text(size = 16, family = "serif", color = "black")) +
  theme(axis.title.x = element_text(size = 16, family = "serif", color = "black")) +
  theme(axis.title.y = element_text(size = 16, family = "serif", color = "black")) +
  theme(legend.text = element_text(size = 14, family = "serif", color = "black")) +
  theme(legend.title = element_text(size = 14, hjust = 0, family = "serif", color = "black")) +
  theme(strip.text = element_text(size = 16, family = "serif", color = "black"))
p <- p + geom_smooth(aes(group = Group), method = "lm", formula = y ~ x, size = 2, se = F, alpha = 0.5)
p <- p + stat_poly_eq(
  formula = y ~ x,
  aes(label = paste(..adj.rr.label.., stat(p.value.label), sep = "*\" \"*")),
  parse = TRUE, label.x.npc = "right", label.y.npc = "bottom", size = 5, family = "serif", lineheight = 1.5
)
p <- p + scale_color_manual(values = c("gray44", "#CD5555", "#0000AA"))
p
##Fig 5A. Zi-Pi plot
#Zi and Pi were calculated using the Cytoscape plugin GIANT (Cumbo et al. 2014).
library(ggsci)
library(ggplot2)
library(ggthemes)
dis <- read.csv("network.phoitc.zi-pi.csv", header = T)
dis$Group <- factor(dis$Group, levels = unique(dis$Group))
p <- ggplot(data = dis, aes(x = Participation.Coefficient, y = Z.Score, color = Group)) +
  geom_point(size = 4) +
  theme_few()
p <- p + labs(x = "Among-module connectivity", y = "Within-module connectivity", title = "")
p <- p + theme(axis.text.x = element_text(size = 16, family = "serif", color = "black")) +
  theme(axis.text.y = element_text(size = 16, family = "serif", color = "black")) +
  theme(axis.title.x = element_text(size = 16, family = "serif", color = "black")) +
  theme(axis.title.y = element_text(size = 16, family = "serif", color = "black")) +
  theme(legend.text = element_text(size = 14, family = "serif", color = "black")) +
  theme(legend.title = element_text(size = 14, hjust = 0, family = "serif", color = "black")) +
  theme(strip.text = element_text(size = 16, family = "serif", color = "black"))
p <- p + geom_hline(yintercept = 2.5, linetype = 2) + geom_vline(xintercept = 0.62, linetype = 2)
p <- p + scale_color_manual(values = c("gray50", "black"))
p
