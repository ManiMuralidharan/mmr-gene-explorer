###############################################################################
# utils.R : shared helper functions
###############################################################################

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

save_session <- function(script_name) {
  f <- file.path(DIR_SESSION, paste0(script_name, "_sessionInfo.txt"))
  writeLines(c(paste("Script:", script_name), paste("Run:", Sys.time()), "",
               capture.output(sessionInfo())), f)
}

# Record each downloaded dataset in metadata/data_sources.csv
record_source <- function(dataset, release, url, access = "open", note = "") {
  f <- file.path(DIR_META, "data_sources.csv")
  row <- data.frame(dataset = dataset, release = release, url = url,
                    access_level = access, access_date = as.character(Sys.Date()),
                    note = note)
  if (file.exists(f)) {
    old <- read.csv(f, stringsAsFactors = FALSE)
    old <- old[old$dataset != dataset, , drop = FALSE]
    row <- rbind(old, row)
  }
  write.csv(row, f, row.names = FALSE)
}

# Find a file inside a folder by regex; error clearly if absent
find_file <- function(dir, pattern, required = TRUE) {
  hits <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!length(hits)) {
    if (required) stop("No file matching '", pattern, "' in ", dir,
                       "\nRun R/00_download_data.R or check the folder contents.")
    return(NA_character_)
  }
  hits[1]
}

# Read a cBioPortal clinical file (skips the 4 '#' header lines)
read_cbio_clinical <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  data.table::fread(text = paste(lines, collapse = "\n"), sep = "\t",
                    data.table = FALSE, check.names = FALSE,
                    na.strings = c("", "NA", "[Not Available]", "[Not Applicable]"))
}

# Read a cBioPortal gene x sample matrix; returns numeric matrix with symbols as rows
read_cbio_matrix <- function(path, genes = NULL) {
  dt <- data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE)
  sym_col <- intersect(c("Hugo_Symbol", "Composite.Element.REF", "Gene"), names(dt))[1]
  if (is.na(sym_col)) sym_col <- names(dt)[1]
  drop <- intersect(c("Hugo_Symbol", "Entrez_Gene_Id", "Composite.Element.REF",
                      "Gene", "Cytoband"), names(dt))
  sym <- as.character(dt[[sym_col]])
  sym <- sub("\\|.*$", "", sym)
  m <- as.matrix(dt[, setdiff(names(dt), drop), drop = FALSE])
  suppressWarnings(storage.mode(m) <- "numeric")
  rownames(m) <- sym
  m <- m[!is.na(rownames(m)) & rownames(m) != "", , drop = FALSE]
  # Collapse duplicate symbols by keeping the row with the highest mean
  if (anyDuplicated(rownames(m))) {
    rm <- rowMeans(m, na.rm = TRUE)
    o <- order(rownames(m), -rm)
    m <- m[o, , drop = FALSE]
    m <- m[!duplicated(rownames(m)), , drop = FALSE]
  }
  if (!is.null(genes)) m <- m[intersect(genes, rownames(m)), , drop = FALSE]
  m
}

# TCGA barcodes: sample-level (first 15 chars) and patient-level (first 12)
tcga_sample <- function(x) substr(x, 1, 15)
tcga_patient <- function(x) substr(x, 1, 12)

# Pairwise correlation table of a target gene against others
cor_table <- function(expr, target, others, method = COR_METHOD, label = "") {
  if (!target %in% rownames(expr)) return(NULL)
  out <- lapply(intersect(others, rownames(expr)), function(g) {
    x <- expr[target, ]; y <- expr[g, ]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
    data.frame(dataset = label, gene_x = target, gene_y = g, n = sum(ok),
               rho = unname(ct$estimate), p_value = ct$p.value)
  })
  out <- do.call(rbind, out)
  if (!is.null(out) && nrow(out)) out$fdr <- p.adjust(out$p_value, "BH")
  out
}

theme_mmr <- function(base = 11) {
  ggplot2::theme_bw(base_size = base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95"))
}
