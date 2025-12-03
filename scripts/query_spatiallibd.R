#!/usr/bin/env Rscript
# Query spatialLIBD for Aging Biology Questions
#
# This script provides functions to query the spatialLIBD brain spatial
# transcriptomics dataset for aging-related genes WITHOUT downloading
# the full 2GB dataset.
#
# Usage:
#   Rscript scripts/query_spatiallibd.R
#   or source in R: source("scripts/query_spatiallibd.R")

# Check and install required packages
check_and_install <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("Installing %s...\n", pkg))
      if (pkg %in% c("spatialLIBD", "SpatialExperiment")) {
        BiocManager::install(pkg, update = FALSE)
      } else {
        install.packages(pkg)
      }
    }
  }
}

# Required packages (minimal set)
required_pkgs <- c("BiocManager")
check_and_install(required_pkgs)

# Load minimal libraries
suppressPackageStartupMessages({
  library(BiocManager)
})

#' Query spatialLIBD for Aging-Related Questions
#'
#' This function provides a framework for asking questions about
#' aging genes in brain spatial transcriptomics data.
#'
#' @param question Character. The question to ask.
#' @param genes Character vector. Gene symbols to query.
#' @param use_web Logical. If TRUE, use web interface (no download).
#'                If FALSE, requires downloading ~2GB dataset.
#' @return Results based on question type
query_aging_genes <- function(question, genes = NULL, use_web = TRUE) {

  cat("\n=== spatialLIBD Aging Biology Query System ===\n\n")
  cat(sprintf("Question: %s\n", question))

  if (is.null(genes)) {
    cat("\nNo genes specified. Loading example aging genes from GenAge...\n")
    genes <- get_example_aging_genes()
  }

  cat(sprintf("\nGenes to query: %s\n", paste(genes, collapse = ", ")))

  if (use_web) {
    cat("\n⚠️  Web interface mode: Using web browser for queries\n")
    cat("    (No 2GB download required)\n\n")

    # Generate web URLs for each gene
    for (gene in genes) {
      url <- sprintf("http://spatial.libd.org/spatialLIBD/?gene=%s", gene)
      cat(sprintf("  • %s: %s\n", gene, url))
    }

    cat("\nRecommendation: Open URLs above to explore spatial expression\n")
    cat("You can:\n")
    cat("  - View spatial expression patterns\n")
    cat("  - Check layer enrichment\n")
    cat("  - Compare across samples\n")
    cat("  - Export figures\n\n")

    return(invisible(list(
      method = "web",
      genes = genes,
      urls = sapply(genes, function(g) sprintf("http://spatial.libd.org/spatialLIBD/?gene=%s", g))
    )))

  } else {
    cat("\n⚠️  DOWNLOAD MODE: This will download ~2GB of data!\n")
    cat("    Are you sure you want to proceed? (y/N): ")

    response <- readline()
    if (tolower(response) != "y") {
      cat("\nCancelled. Use use_web=TRUE for lightweight queries.\n")
      return(invisible(NULL))
    }

    cat("\nDownloading spatialLIBD data (~2GB)...\n")
    cat("This may take several minutes...\n\n")

    # Install spatialLIBD if needed
    if (!requireNamespace("spatialLIBD", quietly = TRUE)) {
      BiocManager::install("spatialLIBD", update = FALSE)
    }

    library(spatialLIBD)
    spe <- fetch_data(type = "spe")

    cat("\nData loaded. Analyzing genes...\n")
    analyze_genes_in_spe(spe, genes, question)
  }
}

#' Get Example Aging Genes from GenAge
get_example_aging_genes <- function() {
  # Top aging-related genes from GenAge
  return(c("TP53", "FOXO3", "IGF1R", "SIRT1", "APOE", "TERT"))
}

#' Analyze genes in SpatialExperiment object
#' (Only runs if full data is downloaded)
analyze_genes_in_spe <- function(spe, genes, question) {

  if (!requireNamespace("spatialLIBD", quietly = TRUE)) {
    stop("spatialLIBD package required for full analysis")
  }

  library(spatialLIBD)
  library(SpatialExperiment)

  cat("\n=== Analysis Results ===\n\n")

  # Check which genes are present
  gene_data <- rowData(spe)
  present_genes <- genes[genes %in% gene_data$gene_name]
  missing_genes <- genes[!genes %in% gene_data$gene_name]

  if (length(missing_genes) > 0) {
    cat(sprintf("⚠️  Genes not found in dataset: %s\n\n",
                paste(missing_genes, collapse = ", ")))
  }

  if (length(present_genes) == 0) {
    cat("No genes found in dataset.\n")
    return(invisible(NULL))
  }

  cat(sprintf("✓ Found %d genes in dataset\n\n", length(present_genes)))

  # Analyze each gene
  for (gene in present_genes) {
    cat(sprintf("--- %s ---\n", gene))

    # Get gene expression
    gene_idx <- which(gene_data$gene_name == gene)
    expr <- assay(spe, "logcounts")[gene_idx, ]

    # Basic statistics
    cat(sprintf("  Mean expression: %.2f\n", mean(expr, na.rm = TRUE)))
    cat(sprintf("  Spots detected: %d / %d (%.1f%%)\n",
                sum(expr > 0), length(expr),
                100 * sum(expr > 0) / length(expr)))

    # Layer enrichment (if available)
    if ("layer_guess_reordered" %in% colnames(colData(spe))) {
      layer_expr <- tapply(expr, colData(spe)$layer_guess_reordered, mean)
      cat("  Layer expression:\n")
      for (layer in names(layer_expr)) {
        cat(sprintf("    %s: %.2f\n", layer, layer_expr[layer]))
      }
    }
    cat("\n")
  }

  return(invisible(list(
    method = "full",
    genes = present_genes,
    spe = spe
  )))
}

#' Predefined Aging Biology Questions
#'
#' Common questions we might ask about aging genes in brain
aging_questions <- list(
  q1 = list(
    question = "Are aging genes layer-specific in the brain?",
    genes = c("TP53", "FOXO3", "IGF1R", "SIRT1"),
    description = "Check if GenAge aging-related genes show preferential expression in specific cortical layers"
  ),

  q2 = list(
    question = "Where are cellular senescence genes expressed in DLPFC?",
    genes = c("CDKN2A", "TP53", "RB1", "CDKN1A", "BCL6"),
    description = "Map CellAge senescence markers to spatial locations in prefrontal cortex"
  ),

  q3 = list(
    question = "Do longevity genes show distinct spatial patterns?",
    genes = c("FOXO3", "APOE", "SIRT1", "TERT"),
    description = "Examine spatial distribution of pro-longevity genes"
  ),

  q4 = list(
    question = "Are immune/inflammation genes enriched in specific layers?",
    genes = c("IL6", "TNF", "NFKB1", "CD68"),
    description = "Spatial pattern of aging-associated inflammation markers"
  ),

  q5 = list(
    question = "How is BCL6 (senescence inhibitor) distributed in brain?",
    genes = c("BCL6"),
    description = "Following up on Shvarts 2002 paper about BCL6 and senescence"
  )
)

#' Run Predefined Question
#'
#' @param question_id Character. Question ID (e.g., "q1", "q2")
#' @param use_web Logical. Use web interface (TRUE) or download data (FALSE)
run_aging_question <- function(question_id, use_web = TRUE) {

  if (!question_id %in% names(aging_questions)) {
    cat("Available questions:\n")
    for (qid in names(aging_questions)) {
      q <- aging_questions[[qid]]
      cat(sprintf("  %s: %s\n", qid, q$question))
    }
    stop(sprintf("Invalid question_id. Use one of: %s",
                 paste(names(aging_questions), collapse = ", ")))
  }

  q <- aging_questions[[question_id]]

  cat("\n" , strrep("=", 70), "\n")
  cat(sprintf("QUESTION %s\n", question_id))
  cat(strrep("=", 70), "\n")
  cat(sprintf("\n%s\n\n", q$question))
  cat(sprintf("Description: %s\n", q$description))
  cat(strrep("-", 70), "\n\n")

  query_aging_genes(
    question = q$question,
    genes = q$genes,
    use_web = use_web
  )
}

#' List All Available Questions
list_aging_questions <- function() {
  cat("\n=== Available Aging Biology Questions for spatialLIBD ===\n\n")

  for (qid in names(aging_questions)) {
    q <- aging_questions[[qid]]
    cat(sprintf("[%s] %s\n", qid, q$question))
    cat(sprintf("     Genes: %s\n", paste(q$genes, collapse = ", ")))
    cat(sprintf("     %s\n\n", q$description))
  }

  cat("Usage:\n")
  cat('  run_aging_question("q1", use_web = TRUE)   # Web interface\n')
  cat('  run_aging_question("q1", use_web = FALSE)  # Full analysis (2GB download)\n\n')
}

#' Query Custom Genes
#'
#' @param genes Character vector of gene symbols
#' @param use_web Logical. Use web interface or download
query_custom_genes <- function(genes, use_web = TRUE) {
  query_aging_genes(
    question = "Custom gene query",
    genes = genes,
    use_web = use_web
  )
}

# ============================================================================
# MAIN EXECUTION (if run as script)
# ============================================================================

if (!interactive()) {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║   spatialLIBD Aging Biology Query System                      ║\n")
  cat("║   Query brain spatial transcriptomics for aging genes         ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n")
  cat("\n")

  # Show available questions
  list_aging_questions()

  cat("Examples:\n")
  cat("  # In R:\n")
  cat('  source("scripts/query_spatiallibd.R")\n')
  cat('  run_aging_question("q1")  # Layer-specific aging genes\n')
  cat('  run_aging_question("q5")  # BCL6 senescence inhibitor\n')
  cat('  query_custom_genes(c("TP53", "CDKN2A"))  # Custom query\n')
  cat("\n")

  cat("💡 TIP: All queries use web interface by default (no download)\n")
  cat("   To do full analysis: run_aging_question('q1', use_web = FALSE)\n")
  cat("   (Warning: Downloads 2GB of data)\n\n")
}

# Print usage message when sourced
cat("\n✓ spatialLIBD query functions loaded.\n")
cat("  Type: list_aging_questions() to see available questions\n")
cat('  Run:  run_aging_question("q1") to execute a query\n\n')
