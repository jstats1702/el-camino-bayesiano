# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(suppressWarnings(library(igraph)))
suppressMessages(suppressWarnings(library(igraphdata)))
suppressMessages(suppressWarnings(library(corrplot)))

# Datos ------------------------------------------------------------------------

data("karate", package = "igraphdata")

karate <- upgrade_graph(karate)

# Versión binaria de la red
g <- delete_edge_attr(karate, "weight")

# Matriz de adyacencia
Y <- as.matrix(
     as_adjacency_matrix(
          g,
          sparse = FALSE
     )
)

dimnames(Y) <- list(
     seq_len(vcount(g)),
     seq_len(vcount(g))
)

# Representación circular ------------------------------------------------------

pdf(
     file      = "karate_grafo_circular.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(0, 0, 0, 0)
)

plot(
     g,
     layout             = layout_in_circle(g),
     vertex.label       = seq_len(vcount(g)),
     vertex.label.color = "black",
     vertex.label.cex   = 0.9,
     vertex.size        = 12,
     vertex.color       = "white",
     vertex.frame.color = "black",
     edge.color         = adjustcolor("blue4", 0.5),
     edge.width         = 0.5,
     edge.curved        = 0,
     main               = ""
)

dev.off()

# Matriz de adyacencia ---------------------------------------------------------

pdf(
     file      = "karate_matriz_adyacencia.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(0, 0, 0, 0)
)

corrplot(
     corr        = Y,
     is.corr     = FALSE,
     method      = "color",
     type        = "full",
     col         = c(
          "white",
          "blue4"
     ),
     col.lim     = c(0, 1),
     tl.pos      = "lt",
     tl.col      = "black",
     tl.cex      = 0.55,
     tl.srt      = 90,
     addgrid.col = NA,
     cl.pos      = "n",
     diag        = TRUE,
     mar         = c(0, 0, 1.5, 0),
     title       = ""
)

rect(
     xleft   = 0.5,
     ybottom = 0.5,
     xright  = ncol(Y) + 0.5,
     ytop    = nrow(Y) + 0.5,
     border  = "black",
     lwd     = 1,
     xpd     = NA
)

dev.off()

# Comparación de diseños -------------------------------------------------------

layout_circle <- layout_in_circle(g)

set.seed(1234)
layout_fr <- layout_with_fr(g)

set.seed(1234)
layout_kk <- layout_with_kk(g)

set.seed(1234)
layout_dh <- layout_with_dh(g)

layouts <- list(
     circular              = layout_circle,
     fruchterman_reingold  = layout_fr,
     kamada_kawai          = layout_kk,
     davidson_harel        = layout_dh
)

# Visualizaciones
for (layout_name in names(layouts)) {
     pdf(
          file      = paste0("karate_diseno_", layout_name, ".pdf"),
          width     = 5,
          height    = 5,
          pointsize = 15
     )
     
     par(
          mar = c(0, 0, 0, 0)
     )
     
     plot(
          g,
          layout             = layouts[[layout_name]],
          vertex.label       = seq_len(vcount(g)),
          vertex.label.color = "black",
          vertex.label.cex   = 0.9,
          vertex.size        = 12,
          vertex.color       = "white",
          vertex.frame.color = "black",
          edge.color         = adjustcolor("blue4", 0.5),
          edge.width         = 0.8,
          edge.curved        = 0
     )
     
     dev.off()
}

# Visualización decorada -------------------------------------------------------

set.seed(123)

l <- layout_with_dh(karate)

# Etiquetas
V(karate)$label <- sub(
     pattern     = "Actor ",
     replacement = "",
     x           = V(karate)$name
)

# Formas
V(karate)$shape <- "circle"
V(karate)[c("Mr Hi", "John A")]$shape <- "rectangle"

# Colores de los vértices
V(karate)[Faction == 1]$color <- "red"
V(karate)[Faction == 2]$color <- "dodgerblue"

# Colores de las aristas
F1 <- V(karate)[Faction == 1]
F2 <- V(karate)[Faction == 2]

E(karate)[F1 %--% F1]$color <- "pink"
E(karate)[F2 %--% F2]$color <- "lightblue"
E(karate)[F1 %--% F2]$color <- "yellow"

# Tamaños de los vértices y las aristas
V(karate)$size <- 7 * sqrt(degree(karate))
E(karate)$width <- E(karate)$weight

# Visualización
pdf(
     file      = paste0("karate_decoracion.pdf"),
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(0, 0, 0, 0)
)

plot(
     karate,
     layout             = l,
     vertex.frame.color = "black",
     vertex.label.color = "black"
)

dev.off()

# Fin --------------------------------------------------------------------------