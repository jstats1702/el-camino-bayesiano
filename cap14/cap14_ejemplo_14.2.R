# Settings ---------------------------------------------------------------------

rm(list = ls())

# Directorio de trabajo
setwd("~/Dropbox/UN/bayes_book")

# Librerías
suppressMessages(suppressWarnings(library(igraph)))
suppressMessages(suppressWarnings(library(igraphdata)))

# Datos ------------------------------------------------------------------------

data("karate", package = "igraphdata")

karate <- upgrade_graph(karate)

# Versión binaria de la red
g <- delete_edge_attr(karate, "weight")

# Medidas de centralidad -------------------------------------------------------

centralities <- list(
     grado = degree(
          graph      = g,
          normalized = TRUE
     ),
     cercania = closeness(
          graph      = g,
          normalized = TRUE
     ),
     intermediacion = betweenness(
          graph      = g,
          normalized = TRUE
     ),
     propia = eigen_centrality(
          graph    = g,
          directed = FALSE
     )$vector
)

# Actores con los cinco valores más altos
sapply(
     centralities,
     function(x) order(x, decreasing = TRUE)[1:5]
)

# Visualizaciones --------------------------------------------------------------

# Diseño
set.seed(123)

layout_centrality <- layout_with_dh(g)

for (centrality_name in names(centralities)) {
     centrality_values <- centralities[[centrality_name]]
     
     vertex_size <- 0.5 +
          18 * sqrt(centrality_values / max(centrality_values))
     
     pdf(
          file = paste0(
               "karate_centralidad_",
               centrality_name,
               ".pdf"
          ),
          width     = 5,
          height    = 5,
          pointsize = 15
     )
     
     par(
          mar = c(0, 0, 0, 0)
     )
     
     plot(
          g,
          layout             = layout_centrality,
          vertex.label       = NA,
          vertex.size        = vertex_size,
          vertex.color       = adjustcolor(
               "royalblue",
               alpha.f = 0.05
          ),
          vertex.frame.color = "royalblue",
          edge.color         = adjustcolor(
               "black",
               alpha.f = 0.5
          ),
          edge.width         = 0.5,
          edge.curved        = 0
     )
     
     dev.off()
}

# Fin --------------------------------------------------------------------------