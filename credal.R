# credal.R
# credal_ranking_helpers.R — no axes, no titles, offset arrows

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggrepel)
  library(grid)
})

build_credal_graph <- function(stan_fit_object, names_vector, confidence_threshold = 0.7){
  alpha_draws <- as.matrix(stan_fit_object, pars = "alpha")
  if (!is.null(colnames(alpha_draws)) && all(names_vector %in% colnames(alpha_draws))) {
    alpha_draws <- alpha_draws[, names_vector, drop = FALSE]
  }
  m <- ncol(alpha_draws)
  prob_matrix <- matrix(NA_real_, m, m, dimnames = list(names_vector, names_vector))
  for (i in seq_len(m)) for (k in seq_len(m)) if (i != k) prob_matrix[i, k] <- mean(alpha_draws[, i] > alpha_draws[, k])
  edge_idx <- which(prob_matrix > confidence_threshold, arr.ind = TRUE)
  edges <- tibble(
    from       = rownames(prob_matrix)[edge_idx[, 1]],
    to         = colnames(prob_matrix)[edge_idx[, 2]],
    confidence = round(prob_matrix[edge_idx], 2)
  )
  nodes <- tibble(name = names_vector)
  tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
}

credal_palette <- c(c1 = "#43A047", c2 = "#304FFE", c3 = "#FFC107", c4 = "#E53935", c5 = "#9E9E9E")
pretty_label <- function(x){
  map <- c(
    c1 = "Current AI use", c2 = "AI as tool", c3 = "AI as assistant", c4 = "AI as threat", c5 = "AI impact",
    c1_v1 = "Awareness of AI Tools", c1_v2 = "AI Tool Competency", c1_v4 = "AI in Library Operations", c1_v7 = "AI Training for Staff",
    c2_v14 = "AI's Impact on Information Literacy", c2_v11 = "AI in Literature Review", c2_v13 = "AI as a Subscription Service", c2_v9 = "Human-in-the-Loop Indexing",
    c2_v8 = "Personalized AI Services", c2_v15 = "Transformative AI Impact", c2_v12 = "Prompt Engineering Collections", c2_v10 = "Automated Cataloguing",
    c3_v17 = "AI in Scientific Writing", c3_v20 = "AI in Library Management", c3_v19 = "AI Agent as Staff Assistant", c3_v18 = "AI Agent in Resource Management",
    c4_v23 = "Resistance to AI Adoption", c4_v25 = "Library Marginalization", c4_v24 = "Staffing Impact",
    c5_v28 = "Importance of Soft Skills", c5_v29 = "Librarian as Critical AI User", c5_v31 = "AI Tool Assessment Service", c5_v30 = "Focus on Human Interaction",
    c5_v26 = "Librarian as AI Consultant", c5_v27 = "AI Tool Co-development"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}
node_domain <- function(x){ sub("^((c[1-5])).*$", "\\1", x) }

plot_credal_ranking <- function(stan_fit_object, names_vector,
                                confidence_threshold = 0.70,
                                node_size = 10,
                                gap_mm = 2.5,
                                arrow_length = unit(4, "mm"),
                                repel_force = 1.2,
                                layout_algo = "fr") {
  g_tbl <- build_credal_graph(stan_fit_object, names_vector, confidence_threshold)
  if (igraph::gsize(g_tbl) == 0) return(invisible(NULL))
  g_tbl <- g_tbl %>%
    activate(nodes) %>% mutate(label = pretty_label(name), dom = node_domain(name)) %>%
    activate(edges) %>% mutate(conf_lab = round(confidence, 2))
  
  node_radius_mm <- node_size / 2
  endcap <- circle(node_radius_mm + gap_mm, "mm")
  
  p <- ggraph(g_tbl, layout = layout_algo) +
    geom_edge_link(aes(label = conf_lab), angle_calc = "along",
                   label_dodge = unit(3, "mm"),
                   arrow = arrow(length = arrow_length, type = "closed"),
                   start_cap = endcap, end_cap = endcap, width = 0.5) +
    geom_node_point(aes(color = dom), size = node_size, show.legend = FALSE) +
    geom_text_repel(aes(x = x, y = y, label = label, color = dom),
                    size = 4.5, box.padding = 0.6, point.padding = 0.5,
                    max.overlaps = Inf, force = repel_force, seed = 42, show.legend = FALSE) +
    scale_color_manual(values = credal_palette, guide = "none") +
    theme_void(base_size = 12) +
    theme(plot.margin = margin(5, 5, 5, 5))
  print(p)
  invisible(p)
}

plot_main_credal <- function(bbwm_rds = "bbwm_stage.rds", confidence_threshold = 0.7,
                             file = "main_credal_ranking.png", width = 7, height = 6.5, res = 300,
                             ...) {
  B <- readRDS(bbwm_rds)
  png(filename = file, width = width, height = height, units = "in", res = res)
  plot_credal_ranking(B$fit_global, B$domains,
                      confidence_threshold = confidence_threshold, ...)
  dev.off()
}

plot_domain_credal <- function(domain, bbwm_rds = "bbwm_stage.rds", confidence_threshold = 0.9,
                               file = NULL, width = 9, height = 9, res = 300, ...) {
  B <- readRDS(bbwm_rds)
  stopifnot(domain %in% names(B$fits_local))
  if (is.null(file)) file <- paste0("domain_credal_ranking_", domain, ".png")
  png(filename = file, width = width, height = height, units = "in", res = res)
  plot_credal_ranking(B$fits_local[[domain]], B$indicators_local[[domain]],
                      confidence_threshold = confidence_threshold, ...)
  dev.off()
}




# Example calls (uncomment to run)
 plot_main_credal(node_size = 12, gap_mm = 4)
# plot_domain_credal("c1")
# plot_domain_credal("c2")
# plot_domain_credal("c3")
# plot_domain_credal("c4")
# plot_domain_credal("c5")
