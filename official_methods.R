# Auxiliary DESeq2 and MAST workflows.
# Inputs: gene-by-cell counts, group factor (B, A), and donor factor.
suppressPackageStartupMessages({
  library(DESeq2)
  library(MAST)
  library(SingleCellExperiment)
})

deseq2_pvalues <- function(counts, group, donor) {
  sample_id <- factor(paste(donor, group, sep = "_"))
  pb_counts <- vapply(levels(sample_id), function(id) {
    rowSums(counts[, sample_id == id, drop = FALSE])
  }, numeric(nrow(counts)))
  rownames(pb_counts) <- rownames(counts)
  colnames(pb_counts) <- levels(sample_id)
  sample_index <- match(levels(sample_id), sample_id)
  sample_group <- droplevels(group[sample_index])
  sample_donor <- droplevels(donor[sample_index])
  keep <- rowSums(pb_counts) > 0
  pvalues <- rep(1, nrow(counts))
  if (sum(keep) < 2L) return(pvalues)
  coldata <- data.frame(donor = sample_donor, condition = sample_group,
                        row.names = colnames(pb_counts))
  dds <- suppressMessages(DESeqDataSetFromMatrix(
    countData = round(pb_counts[keep, , drop = FALSE]),
    colData = coldata, design = ~ donor + condition
  ))
  dds <- DESeq(dds, fitType = "parametric", sfType = "poscounts",
               minReplicatesForReplace = Inf, quiet = TRUE, parallel = FALSE)
  result <- results(dds, contrast = c("condition", "A", "B"),
                    independentFiltering = FALSE, cooksCutoff = FALSE)
  local <- result$pvalue
  local[!is.finite(local)] <- 1
  pvalues[keep] <- local
  pvalues
}

mast_pvalues <- function(counts, group, donor) {
  keep <- rowSums(counts) > 0
  pvalues <- rep(1, nrow(counts))
  if (sum(keep) < 2L) return(pvalues)
  libsize <- pmax(colSums(counts), 1)
  logexpr <- log2(t(t(counts[keep, , drop = FALSE]) / libsize * 1e4) + 1)
  cngeneson <- as.numeric(scale(colSums(counts > 0)))
  cdata <- data.frame(condition = group, donor = donor, cngeneson = cngeneson,
                      row.names = colnames(counts))
  fdata <- data.frame(primerid = rownames(logexpr), row.names = rownames(logexpr))
  sca <- suppressMessages(FromMatrix(exprsArray = logexpr, cData = cdata, fData = fdata))
  invisible(capture.output(
    fit <- suppressMessages(zlm(
      ~ donor + cngeneson + condition, sca,
      silent = TRUE, parallel = FALSE
    ))
  ))
  invisible(capture.output(
    tab <- suppressMessages(summary(fit, doLRT = "conditionA")$datatable)
  ))
  tab <- tab[tab$component == "H" & tab$contrast == "conditionA", ]
  local <- setNames(tab$`Pr(>Chisq)`, tab$primerid)
  matched <- local[rownames(logexpr)]
  matched[!is.finite(matched)] <- 1
  pvalues[keep] <- matched
  pvalues
}
