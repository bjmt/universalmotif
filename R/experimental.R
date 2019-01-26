# UNIVERSALMOTIF CLUSTERS

# Idea:
# Have an object with multiple motifs, and some additional meta info. This meta
# info includes total max/min width, max/min inter-motif width, motif order.
# The idea is to allow for scanning and enriching for clusters of motifs.

universalcluster <- setClass("universalcluster",
                             slots = list(motifs = "list",
                                          meta = "character"))

setValidity("universalcluster",
            function(object) {

              msg <- vector()
              valid <- TRUE

              motifs <- object@motifs
              motif.types <- vapply(motifs, class, character(1))
              if (!all(motif.types == "universalcluster")) {
                msg <- c(msg, "All motifs must be 'universalmotif' objects")
                valid <- FALSE
              }

            })

setMethod("initialize", signature = "universalcluster",
          definition = function(.Object, motifs, meta) {

            .Object@motifs <- motifs
            .Object@meta <- meta

            .Object

          })

setMethod("show", signature = "universalcluster",
          definition = function(object) {
          
            cat("universalcluster object with", length(object@motifs),
                "motifs\n")
            cat("Meta info:\n")
            print(object@meta)
          
          })
