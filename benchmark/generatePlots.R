#!/usr/bin/env Rscript
require("digest")
require("labeling")
require("withr")
require("crayon")
require("ggplot2")

pdf(NULL)

#http://mkweb.bcgsc.ca/biovis2012/
palette <- c("#004949", "#009292", "#ff6db6", "#ffb6db", "#490092",
             "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", "#920000",
             "#924900", "#db6d00", "#24ff24", "#ffff6d", "#000000")

#Paul Tol
paul_tol <- c("#332288", "#117733", "#44aa99", "#88ccee",
              "#ddcc77", "#cc6677", "#aa4499", "#882255")

ibmpalette <- c("#648fff", "#785ef0",
                "#dc267f", "#fe6100",
                "#ffb000"           )


data <- read.csv("results/mkStats.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

precmk <- data.frame(algorithms = "MK",
            value = data$V1 / (data$V1 + data$V2))
recmk <- data.frame(algorithms = "MK",
            value = data$V1 / (data$V1 + data$V3))
precmk <- replace(precmk, is.na(precmk), 0)

data <- read.csv("results/mpStats.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

precmp <- data.frame(algorithms = "MP",
            value = data$V1 / (data$V1 + data$V2))
recmp <- data.frame(algorithms = "MP",
           value = data$V1 / (data$V1 + data$V3))
precmp <- replace(precmp, is.na(precmp), 0)

data <- read.csv("results/emmaStats.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

precemma <- data.frame(algorithms = "EMMA",
              value = data$V1 / (data$V1 + data$V2))
recemma <- data.frame(algorithms = "EMMA",
              value = data$V1 / (data$V1 + data$V3))
precemma <- replace(precemma, is.na(precemma), 0)

data <- read.csv("results/emma2rStats.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

precemma2r <- data.frame(algorithms = "EMMA_2r",
                value = data$V1 / (data$V1 + data$V2))
recemma2r <- data.frame(algorithms = "EMMA_2r",
               value = data$V1 / (data$V1 + data$V3))
precemma2r <- replace(precemma2r, is.na(precemma2r), 0)

data <- read.csv("results/gvStats.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

precgv <- data.frame(algorithms = "GV",
            value = data$V1 / (data$V1 + data$V2))
recgv <- data.frame(algorithms = "GV",
           value = data$V1 / (data$V1 + data$V3))
precgv <- replace(precgv, is.na(precgv), 0)

data <- read.csv("results/lmStats.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

preclm <- data.frame(algorithms = "LM",
            value = data$V1 / (data$V1 + data$V2))
reclm <- data.frame(algorithms = "LM",
           value = data$V1 / (data$V1 + data$V3))
preclm <- replace(preclm, is.na(preclm), 0)

data <- read.csv("results/scanmkStats.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

precscanmk <- data.frame(algorithms = "ScanMK",
                value = data$V1 / (data$V1 + data$V2))
recscanmk <- data.frame(algorithms = "ScanMK",
               value = data$V1 / (data$V1 + data$V3))
precscanmk <- replace(precscanmk, is.na(precscanmk), 0)

data <- read.csv("results/clustermkStats.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

precclustermk <- data.frame(algorithms = "ClusterMK",
                   value = data$V1 / (data$V1 + data$V2))
recclustermk <- data.frame(algorithms = "ClusterMK",
                  value = data$V1 / (data$V1 + data$V3))
precclustermk <- replace(precclustermk, is.na(precclustermk), 0)

data <- read.csv("results/setfinderStats.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

precsetfinder <- data.frame(algorithms = "SetFinder",
                   value = data$V1 / (data$V1 + data$V2))
recsetfinder <- data.frame(algorithms = "SetFinder",
                  value = data$V1 / (data$V1 + data$V3))
precsetfinder <- replace(precsetfinder, is.na(precsetfinder), 0)

plot.data <- rbind(precmk, precmp, precemma, precemma2r, precgv, preclm,
               precscanmk, precclustermk, precsetfinder)

p <- ggplot(plot.data, aes(x = algorithms, y = value, fill = algorithms)) +
       stat_boxplot(geom = "errorbar") +
       labs(x = "Algorithm", y = "Precision", fill = "Algorithms") +
       stat_boxplot(geom = "errorbar") +
       geom_boxplot() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("precPlot.pdf", p, width = 10, height = 4.5, units = "cm", scale = 3)

plot.data <- rbind(recmk, recmp, recemma, recemma2r, recgv, reclm, recscanmk,
               recclustermk, recsetfinder)

p <- ggplot(plot.data, aes(x = algorithms, y = value, fill = algorithms)) +
        stat_boxplot(geom = "errorbar") +
        labs(x = "Algorithm", y = "Recall", fill = "Algorithms") +
        geom_boxplot() +
        scale_colour_manual(values = palette) +
        scale_fill_manual(values = palette) +
        theme(text = element_text(size = 18))
    
ggsave("recPlot.pdf", p, width = 10, height = 4.5, units = "cm", scale = 3)

fscoremk <- data.frame(algorithms = "MK",
              value = 2 * (precmk$value * recmk$value) /
                (precmk$value + recmk$value))
fscoremk <- replace(fscoremk, is.na(fscoremk), 0)

fscoremp <- data.frame(algorithms = "MP",
              value = 2 * (precmp$value * recmp$value) /
                (precmp$value + recmp$value))
fscoremp <- replace(fscoremp, is.na(fscoremp), 0)

fscoreemma <- data.frame(algorithms = "EMMA",
                value = 2 * (precemma$value * recemma$value) /
                         (precemma$value + recemma$value))
fscoreemma <- replace(fscoreemma, is.na(fscoreemma), 0)

fscoreemma2r <- data.frame(
                  algorithms = "EMMA_2r",
                  value = 2 * (precemma2r$value * recemma2r$value) /
                           (precemma2r$value + recemma2r$value))
fscoreemma2r <- replace(fscoreemma2r, is.na(fscoreemma2r), 0)

fscoregv <- data.frame(algorithms = "GV",
              value = 2 * (precgv$value * recgv$value) /
                (precgv$value + recgv$value))
fscoregv <- replace(fscoregv, is.na(fscoregv), 0)

fscorelm <- data.frame(algorithms = "LM",
              value = 2 * (preclm$value * reclm$value) /
                (preclm$value + reclm$value))
fscorelm <- replace(fscorelm, is.na(fscorelm), 0)

fscorescanmk <- data.frame(algorithms = "ScanMK",
                  value = 2 * (precscanmk$value * recscanmk$value) /
                           (precscanmk$value + recscanmk$value))
fscorescanmk <- replace(fscorescanmk, is.na(fscorescanmk), 0)

fscoreclustermk <- data.frame(
                     algorithms = "ClusterMK",
                     value = 2 * (precclustermk$value * recclustermk$value) /
                              (precclustermk$value + recclustermk$value))
fscoreclustermk <- replace(fscoreclustermk, is.na(fscoreclustermk), 0)

fscoresetfinder <- data.frame(algorithms = "SetFinder",
                     value = 2 * (precsetfinder$value * recsetfinder$value) /
                              (precsetfinder$value + recsetfinder$value))
fscoresetfinder <- replace(fscoresetfinder, is.na(fscoresetfinder), 0)

plot.data <- rbind(fscoremk, fscoremp, fscoreemma, fscoreemma2r, fscoregv,
               fscorelm, fscorescanmk, fscoreclustermk, fscoresetfinder)
plot.data$value[is.nan(plot.data$value)] <- 0

p <- ggplot(plot.data, aes(x = algorithms, y = value, fill = algorithms)) +
        stat_boxplot(geom = "errorbar") +
        labs(x = "Algorithm", y = "F1", fill = "Algorithms") +
        geom_boxplot() +
        scale_colour_manual(values = palette) +
        scale_fill_manual(values = palette) +
        theme(text = element_text(size = 18))
    
ggsave("fScorePlot.pdf", p, width = 10, height = 4.5, units = "cm", scale = 3)

precmk1 <- precmk[1:128, ]
recmk1 <- recmk[1:128, ]
fscoremk1 <- fscoremk[1:128, ]

precmp1 <- precmp[1:128, ]
recmp1 <- recmp[1:128, ]
fscoremp1 <- fscoremp[1:128, ]

precemma1 <- precemma[1:128, ]
recemma1 <- recemma[1:128, ]
fscoreemma1 <- fscoreemma[1:128, ]

precgv1 <- precgv[1:128, ]
recgv1 <- recgv[1:128, ]
fscoregv1 <- fscoregv[1:128, ]

preclm1 <- preclm[1:128, ]
reclm1 <- reclm[1:128, ]
fscorelm1 <- fscorelm[1:128, ]

precscanmk1 <- precscanmk[1:128, ]
recscanmk1 <- recscanmk[1:128, ]
fscorescanmk1 <- fscorescanmk[1:128, ]

precclustermk1 <- precclustermk[1:128, ]
recclustermk1 <- recclustermk[1:128, ]
fscoreclustermk1 <- fscoreclustermk[1:128, ]

precsetfinder1 <- precsetfinder[1:128, ]
recsetfinder1 <- recsetfinder[1:128, ]
fscoresetfinder1 <- fscoresetfinder[1:128, ]

precmkpertype <- data.frame(algorithms = "MK",
                    aggregate(list(mean = precmk1[, 2]),
                      list(typeprec = rep(cbind( "box",
                                                 "triangle",
                                                 "semicircle",
                                                 "trapezoid",
                                                 "positiveflank",
                                                 "negativeflank",
                                                 "sine",
                                                 "cosine"
                                               ), rep(16, 8))), mean))
recmkpertype <- data.frame(algorithms = "MK",
                    aggregate(list(mean = recmk1[, 2]),
                      list(typeprec = rep(cbind( "box",
                                                 "triangle",
                                                 "semicircle",
                                                 "trapezoid",
                                                 "positiveflank",
                                                 "negativeflank",
                                                 "sine",
                                                 "cosine"
                                               ), rep(16, 8))), mean))
fscoremkpertype <- data.frame(algorithms = "MK",
                      aggregate(list(mean = fscoremk1[, 2]),
                        list(typeprec = rep(cbind( "box",
                                                   "triangle",
                                                   "semicircle",
                                                   "trapezoid",
                                                   "positiveflank",
                                                   "negativeflank",
                                                   "sine",
                                                   "cosine"
                                                 ), rep(16, 8))), mean))

precmppertype <- data.frame(algorithms = "MP",
                    aggregate(list(mean = precmp1[, 2]),
                      list(typeprec = rep(cbind( "box",
                                                 "triangle",
                                                 "semicircle",
                                                 "trapezoid",
                                                 "positiveflank",
                                                 "negativeflank",
                                                 "sine",
                                                 "cosine"
                                               ), rep(16, 8))), mean))
recmppertype <- data.frame(algorithms = "MP",
                   aggregate(list(mean = recmp1[, 2]),
                     list(typeprec = rep(cbind( "box",
                                                "triangle",
                                                "semicircle",
                                                "trapezoid",
                                                "positiveflank",
                                                "negativeflank",
                                                "sine",
                                                "cosine"
                                              ), rep(16, 8))), mean))
fscoremppertype <- data.frame(algorithms = "MP",
                      aggregate(list(mean = fscoremp1[, 2]),
                        list(typeprec = rep(cbind( "box",
                                                   "triangle",
                                                   "semicircle",
                                                   "trapezoid",
                                                   "positiveflank",
                                                   "negativeflank",
                                                   "sine",
                                                   "cosine"
                                                 ), rep(16, 8))), mean))

precemmapertype <- data.frame(algorithms = "EMMA",
                      aggregate(list(mean = precemma1[, 2]),
                        list(typeprec = rep(cbind( "box",
                                                   "triangle",
                                                   "semicircle",
                                                   "trapezoid",
                                                   "positiveflank",
                                                   "negativeflank",
                                                   "sine",
                                                   "cosine"
                                                 ), rep(16, 8))), mean))
recemmapertype <- data.frame(algorithms = "EMMA",
                     aggregate(list(mean = recemma1[, 2]),
                       list(typeprec = rep(cbind( "box",
                                                  "triangle",
                                                  "semicircle",
                                                  "trapezoid",
                                                  "positiveflank",
                                                  "negativeflank",
                                                  "sine",
                                                  "cosine"
                                                ), rep(16, 8))), mean))
fscoreemmapertype <- data.frame(algorithms = "EMMA",
                        aggregate(list(mean = fscoreemma1[, 2]),
                          list(typeprec = rep(cbind( "box",
                                                     "triangle",
                                                     "semicircle",
                                                     "trapezoid",
                                                     "positiveflank",
                                                     "negativeflank",
                                                     "sine",
                                                     "cosine"
                                                   ), rep(16, 8))), mean))

precgvpertype <- data.frame(algorithms = "GV",
                    aggregate(list(mean = precgv1[, 2]),
                      list(typeprec = rep(cbind( "box",
                                                 "triangle",
                                                 "semicircle",
                                                 "trapezoid",
                                                 "positiveflank",
                                                 "negativeflank",
                                                 "sine",
                                                 "cosine"
                                               ), rep(16, 8))), mean))
recgvpertype <- data.frame(algorithms = "GV",
                   aggregate(list(mean = recgv1[, 2]),
                     list(typeprec = rep(cbind( "box",
                                                "triangle",
                                                "semicircle",
                                                "trapezoid",
                                                "positiveflank",
                                                "negativeflank",
                                                "sine",
                                                "cosine"
                                              ), rep(16, 8))), mean))
fscoregvpertype <- data.frame(algorithms = "GV",
                      aggregate(list(mean = fscoregv1[, 2]),
                        list(typeprec = rep(cbind( "box",
                                                   "triangle",
                                                   "semicircle",
                                                   "trapezoid",
                                                   "positiveflank",
                                                   "negativeflank",
                                                   "sine",
                                                   "cosine"
                                                 ), rep(16, 8))), mean))

preclmpertype <- data.frame(algorithms = "LM",
                    aggregate(list(mean = preclm1[, 2]),
                      list(typeprec = rep(cbind( "box",
                                                 "triangle",
                                                 "semicircle",
                                                 "trapezoid",
                                                 "positiveflank",
                                                 "negativeflank",
                                                 "sine",
                                                 "cosine"
                                               ), rep(16, 8))), mean))
reclmpertype <- data.frame(algorithms = "LM",
                   aggregate(list(mean = reclm1[, 2]),
                     list(typeprec = rep(cbind( "box",
                                                "triangle",
                                                "semicircle",
                                                "trapezoid",
                                                "positiveflank",
                                                "negativeflank",
                                                "sine",
                                                "cosine"
                                              ), rep(16, 8))), mean))
fscorelmpertype <- data.frame(algorithms = "LM",
                      aggregate(list(mean = fscorelm1[, 2]),
                        list(typeprec = rep(cbind( "box",
                                                   "triangle",
                                                   "semicircle",
                                                   "trapezoid",
                                                   "positiveflank",
                                                   "negativeflank",
                                                   "sine",
                                                   "cosine"
                                                 ), rep(16, 8))), mean))

precscanmkpertype <- data.frame(algorithms = "ScanMK",
                        aggregate(list(mean = precscanmk1[, 2]),
                          list(typeprec = rep(cbind( "box",
                                                     "triangle",
                                                     "semicircle",
                                                     "trapezoid",
                                                     "positiveflank",
                                                     "negativeflank",
                                                     "sine",
                                                     "cosine"
                                                   ), rep(16, 8))), mean))
recscanmkpertype <- data.frame(algorithms = "ScanMK",
                       aggregate(list(mean = recscanmk1[, 2]),
                         list(typeprec = rep(cbind( "box",
                                                    "triangle",
                                                    "semicircle",
                                                    "trapezoid",
                                                    "positiveflank",
                                                    "negativeflank",
                                                    "sine",
                                                    "cosine"
                                                  ), rep(16, 8))), mean))
fscorescanmkpertype <- data.frame(algorithms = "ScanMK",
                          aggregate(list(mean = fscorescanmk1[, 2]),
                            list(typeprec = rep(cbind( "box",
                                                       "triangle",
                                                       "semicircle",
                                                       "trapezoid",
                                                       "positiveflank",
                                                       "negativeflank",
                                                       "sine",
                                                       "cosine"
                                                     ), rep(16, 8))), mean))

precclustermkpertype <- data.frame(algorithms = "ClusterMK",
                           aggregate(list(mean = precclustermk1[, 2]),
                             list(typeprec = rep(cbind( "box",
                                                        "triangle",
                                                        "semicircle",
                                                        "trapezoid",
                                                        "positiveflank",
                                                        "negativeflank",
                                                        "sine",
                                                        "cosine"
                                                      ), rep(16, 8))), mean))
recclustermkpertype <- data.frame(algorithms = "ClusterMK",
                          aggregate(list(mean = recclustermk1[, 2]),
                            list(typeprec = rep(cbind( "box",
                                                       "triangle",
                                                       "semicircle",
                                                       "trapezoid",
                                                       "positiveflank",
                                                       "negativeflank",
                                                       "sine",
                                                       "cosine"
                                                     ), rep(16, 8))), mean))
fscoreclustermkpertype <- data.frame(algorithms = "ClusterMK",
                             aggregate(list(mean = fscoreclustermk1[, 2]),
                               list(typeprec = rep(cbind( "box",
                                                          "triangle",
                                                          "semicircle",
                                                          "trapezoid",
                                                          "positiveflank",
                                                          "negativeflank",
                                                          "sine",
                                                          "cosine"
                                                        ), rep(16, 8))), mean))

precsetfinderpertype <- data.frame(algorithms = "SetFinder",
                           aggregate(list(mean = precsetfinder1[, 2]),
                             list(typeprec = rep(cbind( "box",
                                                        "triangle",
                                                        "semicircle",
                                                        "trapezoid",
                                                        "positiveflank",
                                                        "negativeflank",
                                                        "sine",
                                                        "cosine"
                                                       ), rep(16, 8))), mean))
recsetfinderpertype <- data.frame(algorithms = "SetFinder",
                          aggregate(list(mean = recsetfinder1[, 2]),
                            list(typeprec = rep(cbind( "box",
                                                       "triangle",
                                                       "semicircle",
                                                       "trapezoid",
                                                       "positiveflank",
                                                       "negativeflank",
                                                       "sine",
                                                       "cosine"
                                                      ), rep(16, 8))), mean))
fscoresetfinderpertype <- data.frame(algorithms = "SetFinder",
                             aggregate(list(mean = fscoresetfinder1[, 2]),
                               list(typeprec = rep(cbind( "box",
                                                          "triangle",
                                                          "semicircle",
                                                          "trapezoid",
                                                          "positiveflank",
                                                          "negativeflank",
                                                          "sine",
                                                          "cosine"
                                                        ), rep(16, 8))), mean))

plot.data <- rbind(precmkpertype, precmppertype, precemmapertype,
               precgvpertype, preclmpertype, precscanmkpertype,
               precclustermkpertype, precsetfinderpertype)

p <- ggplot(plot.data, aes(x = typeprec, y = mean, colour = algorithms,
       linetype = algorithms, group = algorithms)) +
       labs(x = "Type", y = "Precision",
         colour = "Algorithms", linetype = "Algorithms") +
       geom_line() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("precTypePlot.pdf", p, width = 10, height = 4.5, units = "cm",
  scale = 3)

plot.data <- rbind(recmkpertype, recmppertype, recemmapertype, recgvpertype,
               reclmpertype, recscanmkpertype, recclustermkpertype,
               recsetfinderpertype)

p <- ggplot(plot.data, aes(x = typeprec, y = mean, colour = algorithms,
       linetype = algorithms, group = algorithms)) +
       labs(x = "Type", y = "Recall",
         colour = "Algorithms", linetype = "Algorithms") +
       geom_line() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("recTypePlot.pdf", p, width = 10, height = 4.5, units = "cm",
  scale = 3)

plot.data <- rbind(fscoremkpertype, fscoremppertype, fscoreemmapertype,
               fscoregvpertype, fscorelmpertype, fscorescanmkpertype,
               fscoreclustermkpertype, fscoresetfinderpertype)

p <- ggplot(plot.data, aes(x = typeprec, y = mean, colour = algorithms,
       linetype = algorithms, group = algorithms)) +
       labs(x = "Type", y = "F1",
         colour = "Algorithms", linetype = "Algorithms") +
       geom_line() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("fScoreTypePlot.pdf", p, width = 10, height = 4.5, units = "cm",
  scale = 3)

precmk2 <- precmk[129:256, ]
recmk2 <- recmk[129:256, ]
fscoremk2 <- fscoremk[129:256, ]

precmp2 <- precmp[129:256, ]
recmp2 <- recmp[129:256, ]
fscoremp2 <- fscoremp[129:256, ]

precemma2 <- precemma[129:256, ]
recemma2 <- recemma[129:256, ]
fscoreemma2 <- fscoreemma[129:256, ]

precgv2 <- precgv[129:256, ]
recgv2 <- recgv[129:256, ]
fscoregv2 <- fscoregv[129:256, ]

preclm2 <- preclm[129:256, ]
reclm2 <- reclm[129:256, ]
fscorelm2 <- fscorelm[129:256, ]

precscanmk2 <- precscanmk[129:256, ]
recscanmk2 <- recscanmk[129:256, ]
fscorescanmk2 <- fscorescanmk[129:256, ]

precclustermk2 <- precclustermk[129:256, ]
recclustermk2 <- recclustermk[129:256, ]
fscoreclustermk2 <- fscoreclustermk[129:256, ]

precsetfinder2 <- precsetfinder[129:256, ]
recsetfinder2 <- recsetfinder[129:256, ]
fscoresetfinder2 <- fscoresetfinder[129:256, ]

precmkpernoise <- data.frame(algorithms = "MK",
                    aggregate(list(mean = precmk2[, 2]),
                      list(noiseprec = rep(1:8, rep(16, 8))), mean))
recmkpernoise <- data.frame(algorithms = "MK",
                    aggregate(list(mean = recmk2[, 2]),
                      list(noiseprec = rep(1:8, rep(16, 8))), mean))
fscoremkpernoise <- data.frame(algorithms = "MK",
                      aggregate(list(mean = fscoremk2[, 2]),
                        list(noiseprec = rep(1:8, rep(16, 8))), mean))

precmppernoise <- data.frame(algorithms = "MP",
                    aggregate(list(mean = precmp2[, 2]),
                      list(noiseprec = rep(1:8, rep(16, 8))), mean))
recmppernoise <- data.frame(algorithms = "MP",
                   aggregate(list(mean = recmp2[, 2]),
                     list(noiseprec = rep(1:8, rep(16, 8))), mean))
fscoremppernoise <- data.frame(algorithms = "MP",
                      aggregate(list(mean = fscoremp2[, 2]),
                        list(noiseprec = rep(1:8, rep(16, 8))), mean))

precemmapernoise <- data.frame(algorithms = "EMMA",
                      aggregate(list(mean = precemma2[, 2]),
                        list(noiseprec = rep(1:8, rep(16, 8))), mean))
recemmapernoise <- data.frame(algorithms = "EMMA",
                     aggregate(list(mean = recemma2[, 2]),
                       list(noiseprec = rep(1:8, rep(16, 8))), mean))
fscoreemmapernoise <- data.frame(algorithms = "EMMA",
                        aggregate(list(mean = fscoreemma2[, 2]),
                          list(noiseprec = rep(1:8, rep(16, 8))), mean))

precgvpernoise <- data.frame(algorithms = "GV",
                    aggregate(list(mean = precgv2[, 2]),
                      list(noiseprec = rep(1:8, rep(16, 8))), mean))
recgvpernoise <- data.frame(algorithms = "GV",
                   aggregate(list(mean = recgv2[, 2]),
                     list(noiseprec = rep(1:8, rep(16, 8))), mean))
fscoregvpernoise <- data.frame(algorithms = "GV",
                      aggregate(list(mean = fscoregv2[, 2]),
                        list(noiseprec = rep(1:8, rep(16, 8))), mean))

preclmpernoise <- data.frame(algorithms = "LM",
                    aggregate(list(mean = preclm2[, 2]),
                      list(noiseprec = rep(1:8, rep(16, 8))), mean))
reclmpernoise <- data.frame(algorithms = "LM",
                   aggregate(list(mean = reclm2[, 2]),
                     list(noiseprec = rep(1:8, rep(16, 8))), mean))
fscorelmpernoise <- data.frame(algorithms = "LM",
                      aggregate(list(mean = fscorelm2[, 2]),
                        list(noiseprec = rep(1:8, rep(16, 8))), mean))

precscanmkpernoise <- data.frame(algorithms = "ScanMK",
                        aggregate(list(mean = precscanmk2[, 2]),
                          list(noiseprec = rep(1:8, rep(16, 8))), mean))
recscanmkpernoise <- data.frame(algorithms = "ScanMK",
                       aggregate(list(mean = recscanmk2[, 2]),
                         list(noiseprec = rep(1:8, rep(16, 8))), mean))
fscorescanmkpernoise <- data.frame(algorithms = "ScanMK",
                          aggregate(list(mean = fscorescanmk2[, 2]),
                            list(noiseprec = rep(1:8, rep(16, 8))), mean))

precclustermkpernoise <- data.frame(algorithms = "ClusterMK",
                           aggregate(list(mean = precclustermk2[, 2]),
                             list(noiseprec = rep(1:8, rep(16, 8))), mean))
recclustermkpernoise <- data.frame(algorithms = "ClusterMK",
                          aggregate(list(mean = recclustermk2[, 2]),
                            list(noiseprec = rep(1:8, rep(16, 8))), mean))
fscoreclustermkpernoise <- data.frame(algorithms = "ClusterMK",
                             aggregate(list(mean = fscoreclustermk2[, 2]),
                               list(noiseprec = rep(1:8, rep(16, 8))), mean))

precsetfinderpernoise <- data.frame(algorithms = "SetFinder",
                           aggregate(list(mean = precsetfinder2[, 2]),
                             list(noiseprec = rep(1:8, rep(16, 8))), mean))
recsetfinderpernoise <- data.frame(algorithms = "SetFinder",
                          aggregate(list(mean = recsetfinder2[, 2]),
                            list(noiseprec = rep(1:8, rep(16, 8))), mean))
fscoresetfinderpernoise <- data.frame(algorithms = "SetFinder",
                             aggregate(list(mean = fscoresetfinder2[, 2]),
                               list(noiseprec = rep(1:8, rep(16, 8))), mean))

plot.data <- rbind(precmkpernoise, precmppernoise, precemmapernoise,
               precgvpernoise, preclmpernoise, precscanmkpernoise,
               precclustermkpernoise, precsetfinderpernoise)

p <- ggplot(plot.data, aes(x = noiseprec, y = mean, colour = algorithms,
       linetype = algorithms)) +
       labs(x = "Noise in %", y = "Precision",
         colour = "Algorithms", linetype = "Algorithms") +
       geom_line() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("precNoisePlot.pdf", p, width = 10, height = 4.5, units = "cm",
  scale = 3)

plot.data <- rbind(recmkpernoise, recmppernoise, recemmapernoise, recgvpernoise,
               reclmpernoise, recscanmkpernoise, recclustermkpernoise,
               recsetfinderpernoise)

p <- ggplot(plot.data, aes(x = noiseprec, y = mean, colour = algorithms,
       linetype = algorithms)) +
       labs(x = "Noise in %", y = "Recall",
         colour = "Algorithms", linetype = "Algorithms") +
       geom_line() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("recNoisePlot.pdf", p, width = 10, height = 4.5, units = "cm",
  scale = 3)

plot.data <- rbind(fscoremkpernoise, fscoremppernoise, fscoreemmapernoise,
               fscoregvpernoise, fscorelmpernoise, fscorescanmkpernoise,
               fscoreclustermkpernoise, fscoresetfinderpernoise)

p <- ggplot(plot.data, aes(x = noiseprec, y = mean, colour = algorithms,
       linetype = algorithms)) +
       labs(x = "Noise in %", y = "F1",
         colour = "Algorithms", linetype = "Algorithms") +
       geom_line() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("fScoreNoisePlot.pdf", p, width = 10, height = 4.5, units = "cm",
  scale = 3)

precmk3 <- precmk[257:384, ]
recmk3 <- recmk[257:384, ]
fscoremk3 <- fscoremk[257:384, ]

precmp3 <- precmp[257:384, ]
recmp3 <- recmp[257:384, ]
fscoremp3 <- fscoremp[257:384, ]

precemma3 <- precemma[257:384, ]
recemma3 <- recemma[257:384, ]
fscoreemma3 <- fscoreemma[257:384, ]

precgv3 <- precgv[257:384, ]
recgv3 <- recgv[257:384, ]
fscoregv3 <- fscoregv[257:384, ]

preclm3 <- preclm[257:384, ]
reclm3 <- reclm[257:384, ]
fscorelm3 <- fscorelm[257:384, ]

precscanmk3 <- precscanmk[257:384, ]
recscanmk3 <- recscanmk[257:384, ]
fscorescanmk3 <- fscorescanmk[257:384, ]

precclustermk3 <- precclustermk[257:384, ]
recclustermk3 <- recclustermk[257:384, ]
fscoreclustermk3 <- fscoreclustermk[257:384, ]

precsetfinder3 <- precsetfinder[257:384, ]
recsetfinder3 <- recsetfinder[257:384, ]
fscoresetfinder3 <- fscoresetfinder[257:384, ]

precmkpersize <- data.frame(algorithms = "MK",
                    aggregate(list(mean = precmk3[, 2]),
                      list(sizeprec = rep(3:10, rep(16, 8))), mean))
recmkpersize <- data.frame(algorithms = "MK",
                    aggregate(list(mean = recmk3[, 2]),
                      list(sizeprec = rep(3:10, rep(16, 8))), mean))
fscoremkpersize <- data.frame(algorithms = "MK",
                      aggregate(list(mean = fscoremk3[, 2]),
                        list(sizeprec = rep(3:10, rep(16, 8))), mean))

precmppersize <- data.frame(algorithms = "MP",
                    aggregate(list(mean = precmp3[, 2]),
                      list(sizeprec = rep(3:10, rep(16, 8))), mean))
recmppersize <- data.frame(algorithms = "MP",
                   aggregate(list(mean = recmp3[, 2]),
                     list(sizeprec = rep(3:10, rep(16, 8))), mean))
fscoremppersize <- data.frame(algorithms = "MP",
                      aggregate(list(mean = fscoremp3[, 2]),
                        list(sizeprec = rep(3:10, rep(16, 8))), mean))

precemmapersize <- data.frame(algorithms = "EMMA",
                      aggregate(list(mean = precemma3[, 2]),
                        list(sizeprec = rep(3:10, rep(16, 8))), mean))
recemmapersize <- data.frame(algorithms = "EMMA",
                     aggregate(list(mean = recemma3[, 2]),
                       list(sizeprec = rep(3:10, rep(16, 8))), mean))
fscoreemmapersize <- data.frame(algorithms = "EMMA",
                        aggregate(list(mean = fscoreemma3[, 2]),
                          list(sizeprec = rep(3:10, rep(16, 8))), mean))

precgvpersize <- data.frame(algorithms = "GV",
                    aggregate(list(mean = precgv3[, 2]),
                      list(sizeprec = rep(3:10, rep(16, 8))), mean))
recgvpersize <- data.frame(algorithms = "GV",
                   aggregate(list(mean = recgv3[, 2]),
                     list(sizeprec = rep(3:10, rep(16, 8))), mean))
fscoregvpersize <- data.frame(algorithms = "GV",
                      aggregate(list(mean = fscoregv3[, 2]),
                        list(sizeprec = rep(3:10, rep(16, 8))), mean))

preclmpersize <- data.frame(algorithms = "LM",
                    aggregate(list(mean = preclm3[, 2]),
                      list(sizeprec = rep(3:10, rep(16, 8))), mean))
reclmpersize <- data.frame(algorithms = "LM",
                   aggregate(list(mean = reclm3[, 2]),
                     list(sizeprec = rep(3:10, rep(16, 8))), mean))
fscorelmpersize <- data.frame(algorithms = "LM",
                      aggregate(list(mean = fscorelm3[, 2]),
                        list(sizeprec = rep(3:10, rep(16, 8))), mean))

precscanmkpersize <- data.frame(algorithms = "ScanMK",
                        aggregate(list(mean = precscanmk3[, 2]),
                          list(sizeprec = rep(3:10, rep(16, 8))), mean))
recscanmkpersize <- data.frame(algorithms = "ScanMK",
                       aggregate(list(mean = recscanmk3[, 2]),
                         list(sizeprec = rep(3:10, rep(16, 8))), mean))
fscorescanmkpersize <- data.frame(algorithms = "ScanMK",
                          aggregate(list(mean = fscorescanmk3[, 2]),
                            list(sizeprec = rep(3:10, rep(16, 8))), mean))

precclustermkpersize <- data.frame(algorithms = "ClusterMK",
                           aggregate(list(mean = precclustermk3[, 2]),
                             list(sizeprec = rep(3:10, rep(16, 8))), mean))
recclustermkpersize <- data.frame(algorithms = "ClusterMK",
                          aggregate(list(mean = recclustermk3[, 2]),
                            list(sizeprec = rep(3:10, rep(16, 8))), mean))
fscoreclustermkpersize <- data.frame(algorithms = "ClusterMK",
                             aggregate(list(mean = fscoreclustermk3[, 2]),
                               list(sizeprec = rep(3:10, rep(16, 8))), mean))

precsetfinderpersize <- data.frame(algorithms = "SetFinder",
                           aggregate(list(mean = precsetfinder3[, 2]),
                             list(sizeprec = rep(3:10, rep(16, 8))), mean))
recsetfinderpersize <- data.frame(algorithms = "SetFinder",
                          aggregate(list(mean = recsetfinder3[, 2]),
                            list(sizeprec = rep(3:10, rep(16, 8))), mean))
fscoresetfinderpersize <- data.frame(algorithms = "SetFinder",
                             aggregate(list(mean = fscoresetfinder3[, 2]),
                               list(sizeprec = rep(3:10, rep(16, 8))), mean))

plot.data <- rbind(precmkpersize, precmppersize, precemmapersize,
               precgvpersize, preclmpersize, precscanmkpersize,
               precclustermkpersize, precsetfinderpersize)

p <- ggplot(plot.data, aes(x = sizeprec, y = mean, colour = algorithms,
       linetype = algorithms)) +
       labs(x = "Number of Injected Subsequences", y = "Precision",
         colour = "Algorithms", linetype = "Algorithms") +
       geom_line() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("precSizePlot.pdf", p, width = 10, height = 4.5, units = "cm",
  scale = 3)

plot.data <- rbind(recmkpersize, recmppersize, recemmapersize, recgvpersize,
               reclmpersize, recscanmkpersize, recclustermkpersize,
               recsetfinderpersize)

p <- ggplot(plot.data, aes(x = sizeprec, y = mean, colour = algorithms,
       linetype = algorithms)) +
       labs(x = "Number of Injected Subsequences", y = "Recall",
         colour = "Algorithms", linetype = "Algorithms") +
       geom_line() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("recSizePlot.pdf", p, width = 10, height = 4.5, units = "cm",
  scale = 3)

plot.data <- rbind(fscoremkpersize, fscoremppersize, fscoreemmapersize,
               fscoregvpersize, fscorelmpersize, fscorescanmkpersize,
               fscoreclustermkpersize, fscoresetfinderpersize)

p <- ggplot(plot.data, aes(x = sizeprec, y = mean, colour = algorithms,
       linetype = algorithms)) +
       labs(x = "Number of Injected Subsequences", y = "F1",
         colour = "Algorithms", linetype = "Algorithms") +
       geom_line() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))

ggsave("fScoreSizePlot.pdf", p, width = 10, height = 4.5, units = "cm",
  scale = 3)

data <- read.csv("results/mkRuntimes.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

timemk <- data.frame(algorithms = "MK", value = data$V1)

data <- read.csv("results/mpRuntimes.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

timemp <- data.frame(algorithms = "MP", value = data$V1)

data <- read.csv("results/emmaRuntimes.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

timeemma <- data.frame(algorithms = "EMMA", value = data$V1)

data <- read.csv("results/gvRuntimes.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

timegv <- data.frame(algorithms = "GV", value = data$V1)

data <- read.csv("results/lmRuntimes.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

timelm <- data.frame(algorithms = "LM", value = data$V1)

data <- read.csv("results/scanmkRuntimes.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

timescanmk <- data.frame(algorithms = "ScanMK", value = data$V1)

data <- read.csv("results/clustermkRuntimes.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

timeclustermk <- data.frame(algorithms = "ClusterMK", value = data$V1)

data <- read.csv("results/setfinderRuntimes.csv", header = FALSE, sep = ",",
          quote = "\"", dec = ".", fill = TRUE, comment.char = "")

timesetfinder <- data.frame(algorithms = "SetFinder", value = data$V1)

plot.data <- rbind(timemk, timemp, timeemma, timegv, timelm, timescanmk,
               timeclustermk, timesetfinder)

p <- ggplot(plot.data, aes(x = algorithms, y = value, fill = algorithms)) +
       stat_boxplot(geom = "errorbar") +
       labs(x = "Algorithm", y = "Runtime in sec", fill = "Algorithms") +
       scale_y_continuous(trans = "log2") +
       geom_boxplot() +
       scale_colour_manual(values = palette) +
       scale_fill_manual(values = palette) +
       theme(text = element_text(size = 18))
    
ggsave("timePlot.pdf", p, width = 10, height = 4.5, units = "cm", scale = 3)
