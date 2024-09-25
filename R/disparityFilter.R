install.packages("disparityfilter")
library(igraph)
library(disparityfilter)
g <- sample_pa(n = 250, m = 5, directed = FALSE)
E(g)$weight <- sample(1:25, ecount(g), replace = TRUE)
backbone(g)




devtools::load_all("~/syncomBuildR")
taxa <- c(
  "Bacteria", "Proteobacteria", "Betaproteobacteria", "Burkholderiales",
  "Burkholderiaceae", "Paraburkholderia", NA
)
taxa <- matrix(rep(taxa, 10), nrow = 10, byrow = TRUE)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxa) <- paste0("ASV", 1:10)
# taxonomy data if used should have ASV names explicitly as a column
taxa_df <- as.data.frame(taxa)
taxa_df$asv <- rownames(taxa_df)
sp_dist <- asvDist(asv, method = "spearman", clr_transform = TRUE, edgeFilter = 0.5)
net_data <- asvNet(sp_dist, taxa_df, edge = "spearman")
net_data <- netClust(net = net_data, "components")
asv$biomass_z <- rnorm(nrow(asv))

# syncombuildr disparity filter option

disp_filt <- function() {
  # params
  net_data <- net_data
  weights <- NULL
  alpha = 0.05
  node_metric <- "degree"
  
  # calculated values
  weights <- .check_weights(net_data, weights)
  use_names <- igraph::is_named(G)
  net_met <- .check_node_metric(net_data, node_metric)
  
  d = degree(net_data$graph, mode = mode)
  e = cbind(as_data_frame(G)[, 1:2 ], weight = weights, alpha = NA)
  
  
}

.check_node_metric <- function() {
  
}

.check_weights <- function(net_data, weights) {
  if (is.null(weights)) { # put this in a helper function
    weights <- igraph::E(net_data$graph)$weight
  } else if (length(weights) == 1 && is.character(weights)) {
    weights <- net_data$edges[[weights]]
  } else if (!is.numeric(weights)) {
    stop("weights must be a numeric vector, NULL, or a single column name from the edges data")
  }
  return(weights)
}


G <- net_data$graph
#G <- g
weights = E(G)$weight
directed = igraph::is_directed(G)
alpha = 0.05
mode = "all"
disparity_filter(G, weights, mode, alpha = 0.1)

disparity_filter <- function(G, weights, mode = "all", alpha = 0.05) {
  
  d = degree(G, mode = mode)
  e = cbind(as_data_frame(G)[, 1:2 ], weight = weights, alpha = NA)
  if (mode == "all") {
    e = rbind(e, data.frame(from = e[["to"]], to = e[["from"]], e[, c("weight", "alpha")]))
  }
  
  for (u in which(d > 1)) {
    
    w = switch(substr(mode, 1, 1),
               a = which(e[, 1] == u | e[, 2] == u),
               i = which(e[, 2] == u),
               o = which(e[, 1] == u)
    )
    w = sum(e$weight[ w ]) / (1 + (mode == "all"))
    
    k = d[u]
    
    for (v in igraph::ego(G, 1, u, mode)[[1]][-1]) {
      
      ij = switch(substr(mode, 1, 1),
                  a = which(e[, 1] == u & e[, 2] == v),
                  i = which(e[, 1] == v & e[, 2] == u),
                  o = which(e[, 1] == u & e[, 2] == v)
      )
      e$alpha[ ij ] = (1 - e$weight[ ij ] / w) ^ (k - 1)
      
    }
    
  }
  
  return(e[ !is.na(e$alpha) & e$alpha < alpha, 1:4 ])
  
}


