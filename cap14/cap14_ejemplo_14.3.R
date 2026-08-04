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

# Medidas estructurales --------------------------------------------------------

density_value <- edge_density(
     graph = g,
     loops = FALSE
)

round(density_value, 3)

transitivity_value <- transitivity(
     graph = g,
     type  = "global"
)

round(transitivity_value, 3)

assortativity_degree_value <- assortativity_degree(
     graph    = g,
     directed = FALSE
)

round(assortativity_degree_value, 3)

# Detección de comunidades -----------------------------------------------------

fast_greedy <- cluster_fast_greedy(graph = g)

community_membership <- igraph::membership(fast_greedy)
community_groups     <- igraph::groups(fast_greedy)
community_sizes      <- igraph::sizes(fast_greedy)

modularity_value <- modularity(fast_greedy)

round(modularity_value, 3)

# Resultados -------------------------------------------------------------------

network_summary <- data.frame(
     medida = c(
          "Densidad",
          "Transitividad global",
          "Asortatividad por grado",
          "Número de comunidades",
          "Modularidad"
     ),
     valor = c(
          density_value,
          transitivity_value,
          assortativity_degree_value,
          length(community_groups),
          modularity_value
     )
)

print(network_summary)
print(community_sizes)

# Elementos visuales -----------------------------------------------------------

# Diseño
set.seed(123)

layout_communities <- layout_with_fr(graph = g)

n_communities <- length(community_groups)

community_palette <- hcl.colors(
     n       = n_communities,
     palette = "Dark 3"
)

vertex_colors <- adjustcolor(
     community_palette[community_membership],
     alpha.f = 0.35
)

vertex_frames <- community_palette[community_membership]

vertex_sizes <- 6 + 3 * sqrt(
     degree(g)
)

edge_ends <- ends(
     graph = g,
     es    = E(g),
     names = FALSE
)

edge_community_1 <- community_membership[edge_ends[, 1]]
edge_community_2 <- community_membership[edge_ends[, 2]]

edge_colors <- ifelse(
     edge_community_1 == edge_community_2,
     adjustcolor(
          community_palette[edge_community_1],
          alpha.f = 0.45
     ),
     adjustcolor(
          "gray40",
          alpha.f = 0.35
     )
)

group_colors <- adjustcolor(
     community_palette,
     alpha.f = 0.10
)

# Visualización ----------------------------------------------------------------

pdf(
     file      = "karate_comunidades_fast_greedy.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(0, 0, 0, 0)
)

plot(
     g,
     layout             = layout_communities,
     vertex.label       = seq_len(vcount(g)),
     vertex.label.color = "black",
     vertex.label.cex   = 0.65,
     vertex.size        = vertex_sizes,
     vertex.color       = vertex_colors,
     vertex.frame.color = vertex_frames,
     edge.color         = edge_colors,
     edge.width         = 1,
     edge.curved        = 0
)

dev.off()

# Fin --------------------------------------------------------------------------