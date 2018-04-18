# test

motifdb <- convert_motifs(MotifDb::MotifDb)
motifdbraw <- lapply(motifdb, function(x) motif_slots(x, "motif"))

jaspar.scores <- MotIV::readDBScores(file.path(find.package("MotIV"),
                                               "extdata", "jaspar2010_PCC_SWU.scores"))

d <- MotIV::motifDistances(motifdbraw, DBscores = jaspar.scores)
hc <- MotIV::motifHclust(d, method = "average")

mot_org <- vapply(motifdb, function(x) motif_slots(x, "extrachar")["extrachar.organism"],
                  character(1))

library(RColorBrewer)
n <- length(mot_org) 
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

mot_cols <- vector(mode = "character", length = length(mot_org))
for (i in seq_along(unique(mot_org))) {
  mot_cols[mot_org == unique(mot_org)[i]] <- col_vector[i]
}

library(ape)
plot(as.phylo(hc), type = "unrooted", cex = 0.1, no.margin = TRUE,
     show.tip.label = FALSE, 
     edge.color = sample(mot_cols, length(as.phylo(hc)$edge) / 2, replace = TRUE))

library(dendextend)
library(circlize)

motnames <- names(motifdbraw)
mot_org2 <- mot_org[match(labels(dend), motnames)]

# num_clades <- length(mot_org)

hc2 <- hc
hc2$labels <- mot_org

dend <- as.dendrogram(hc2)
dend <- dend %>% 
  # color_branches(k = num_clades, col = rainbow) #%>%
  color_branches(col = 1:length(unique(mot_org)), 
                 clusters = as.integer(as.factor(mot_org)))
  # color_labels(col = 1:length(unique(mot_org)), labels = labels(dend))
  # labels_cex(value = 0.1)
  # color_labels(k = num_clades, col = rainbow) %>%
  # set("leaves_pch", 19) %>%
  # set("leaves_col", c(1:length(mot_org)))

# par(mar = rep(0, 4))
circlize_dendrogram(dend, dend_track_height = 0.8, labels = FALSE)
