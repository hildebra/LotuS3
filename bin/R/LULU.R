# uses LULU to clean OTU tables
# OPTIMIZED VERSION - drop-in replacement for LotuS3/bin/R/LULU.R
#
# Key optimizations vs. the original:
#  1. Matchlist is pre-indexed by OTUid (split() into a named list) so the
#     per-OTU hit lookup becomes O(1) instead of an O(N) scan of a data.frame.
#  2. The OTU table is converted to a numeric matrix once. All per-OTU lookups
#     use matrix row indexing (~10-100x faster than data.frame indexing).
#  3. Row names are mapped to integer indices via a named vector so
#     "which(rownames %in% hits)" becomes a vector lookup.
#  4. Logging is buffered: each OTU writes a single string to the log file
#     instead of ~10 cat() calls. Reduces I/O dramatically.
#  5. statistics_table is kept lean during the loop (we only need parent_id);
#     decorative columns are added at the end.
#  6. Per-pair co-occurrence and abundance ratios are computed with vectorised
#     primitives on raw numeric vectors instead of data.frame subsetting.
#
# Behaviour is intended to be identical to the original (same outputs, same
# log content). Compare via diff on output files when testing.
# saurabhthakar3@gmail.com

if (!require("dplyr", quietly = TRUE, warn.conflicts = FALSE)) {
    install.packages("dplyr")
    require("dplyr", warn.conflicts = FALSE)
}


lulu = function(otutable, matchlist,
                minimum_ratio_type = "min",
                minimum_ratio = 1,
                minimum_match = 84,
                minimum_relative_cooccurence = 0.95) {

    suppressWarnings(require(dplyr, warn.conflicts = FALSE))
    start.time <- Sys.time()

    # ---- 1. Clean & pre-index the matchlist ----------------------------------
    colnames(matchlist) <- c("OTUid", "hit", "match")
    matchlist <- matchlist[matchlist$hit != "*" &
                           matchlist$hit != matchlist$OTUid &
                           matchlist$match > minimum_match, , drop = FALSE]

    # Pre-split hits by OTUid so per-OTU lookup is O(1).
    # hits_by_otu[["OTU_xyz"]] returns a character vector of hit IDs.
    hits_by_otu <- split(matchlist$hit, matchlist$OTUid)

    # ---- 2. Prepare OTU table & statistics ----------------------------------
    otutable <- otutable[rowSums(otutable) > 0, , drop = FALSE]

    # Build statistics_table with the same columns as before and sort.
    total_v  <- rowSums(otutable)
    spread_v <- rowSums(otutable > 0)
    statistics_table <- data.frame(
        total  = total_v,
        spread = spread_v,
        row.names = rownames(otutable),
        stringsAsFactors = FALSE
    )
    ord <- order(statistics_table$spread, statistics_table$total,
                 decreasing = TRUE)
    statistics_table <- statistics_table[ord, , drop = FALSE]
    otutable <- otutable[match(rownames(statistics_table),
                               rownames(otutable)), , drop = FALSE]
    statistics_table$parent_id <- "NA"

    # ---- 3. Convert to matrix for fast row access ---------------------------
    otu_mat   <- as.matrix(otutable)
    otu_names <- rownames(otu_mat)
    n_otu     <- nrow(otu_mat)

    # Row-name -> row-index lookup (named integer vector).
    name_to_idx <- setNames(seq_len(n_otu), otu_names)

    # Precompute presence/absence (TRUE/FALSE matrix) and per-row spread.
    pres_mat    <- otu_mat > 0
    spread_vec  <- statistics_table$spread

    # parent_id as a plain character vector for fast writes (assign back at end)
    parent_id_v <- statistics_table$parent_id

    # ---- 4. Open log & set up buffered logging ------------------------------
    # The verbose per-OTU log is the largest remaining cost. It can be skipped
    # by setting environment variable LULU_NO_LOG=1, which gives a further
    # ~2x speedup on top of the algorithmic optimizations and is safe for
    # production runs that don't need the diagnostic log.
    no_log <- nzchar(Sys.getenv("LULU_NO_LOG"))
    if (!no_log) {
        log_con <- file(file.path(logD,
                                  paste0("lulu.log_",
                                         format(start.time, "%Y%m%d_%H%M%S"))),
                        open = "a")
    }

    tarProg <- 0.1

    # ---- 5. Main loop -------------------------------------------------------
    for (line in seq_len(n_otu)) {

        if (line / n_otu > tarProg) {
            print(paste0("progress: ",
                         round(((line / n_otu) * 100), 0), "%"))
            tarProg <- tarProg + 0.1
        }

        potential_parent_id <- otu_names[line]

        # Daughter sample vector (numeric, no data.frame overhead)
        daughter_samples <- otu_mat[line, ]
        daughter_pres    <- pres_mat[line, ]
        daughter_n_pres  <- sum(daughter_pres)

        # O(1) lookup of hits (returns NULL if OTU has no hits)
        hits <- hits_by_otu[[potential_parent_id]]

        # Build the log entry for this OTU as a single string (one I/O write).
        if (!no_log) {
            log_buf <- paste0("\n####processing: ", potential_parent_id, " #####")
            if (length(hits) > 0) {
                log_buf <- paste0(log_buf,
                                  paste0("\n---hits: ", hits, collapse = ""))
            }
        }

        # Restrict potential parents to OTUs with spread >= current OTU's spread.
        # Original logic: row.names(otutable)[1:last_relevant_entry] %in% hits.
        # Since the table is sorted by (spread, total) decreasing, the spread
        # cutoff is contiguous at the top.
        last_relevant_entry <- sum(spread_vec >= spread_vec[line])

        # Fast hit -> row index lookup using the named vector.
        potential_parents <- integer(0)
        if (length(hits) > 0) {
            cand_idx <- name_to_idx[hits]
            cand_idx <- cand_idx[!is.na(cand_idx) &
                                 cand_idx <= last_relevant_entry]
            if (length(cand_idx) > 0) {
                # Preserve the original ordering (sorted ascending by row idx,
                # i.e. by rank in statistics_table).
                potential_parents <- sort(unname(cand_idx))
                if (!no_log) {
                    log_buf <- paste0(log_buf,
                                      paste0("\n---potential parent: ",
                                             otu_names[potential_parents],
                                             collapse = ""))
                }
            }
        }

        success <- FALSE

        if (length(potential_parents) > 0 && daughter_n_pres > 0) {
            for (line2 in potential_parents) {
                if (!no_log) {
                    log_buf <- paste0(log_buf,
                                      "\n------checking: ", otu_names[line2])
                }

                if (!success) {
                    parent_pres <- pres_mat[line2, ]

                    # Co-occurrence: of samples where daughter is present,
                    # what fraction also has the parent present?
                    rel_cooc <- sum(daughter_pres & parent_pres) /
                                daughter_n_pres
                    if (!no_log) {
                        log_buf <- paste0(log_buf,
                                          "\n------relative cooccurence: ",
                                          rel_cooc)
                    }

                    if (rel_cooc >= minimum_relative_cooccurence) {
                        if (!no_log) {
                            log_buf <- paste0(log_buf, " which is sufficient!")
                        }

                        # Abundance ratio: parent/daughter across samples where
                        # daughter is present. Vectorised over a numeric slice.
                        parent_at_d <- otu_mat[line2, daughter_pres]
                        d_at_d      <- daughter_samples[daughter_pres]
                        ratios      <- parent_at_d / d_at_d

                        if (minimum_ratio_type == "avg") {
                            relative_abundance <- mean(ratios)
                            if (!no_log) {
                                log_buf <- paste0(log_buf,
                                    "\n------mean avg abundance: ",
                                    relative_abundance)
                            }
                        } else {
                            relative_abundance <- min(ratios)
                            if (!no_log) {
                                log_buf <- paste0(log_buf,
                                    "\n------min avg abundance: ",
                                    relative_abundance)
                            }
                        }

                        if (relative_abundance > minimum_ratio) {
                            if (!no_log) {
                                log_buf <- paste0(log_buf, " which is OK!")
                            }
                            if (line2 < line) {
                                parent_id_v[line] <- parent_id_v[line2]
                                if (!no_log) {
                                    log_buf <- paste0(log_buf,
                                        "\nSETTING ", potential_parent_id,
                                        " to be an ERROR of ",
                                        parent_id_v[line2], "\n")
                                }
                            } else {
                                parent_id_v[line] <- otu_names[line2]
                                if (!no_log) {
                                    log_buf <- paste0(log_buf,
                                        "\nSETTING ", potential_parent_id,
                                        " to be an ERROR of ",
                                        otu_names[line2], "\n")
                                }
                            }
                            success <- TRUE
                        }
                    }
                }
            }
        }

        if (!success) {
            parent_id_v[line] <- otu_names[line]
            if (!no_log) {
                log_buf <- paste0(log_buf, "\nNo parent found!\n")
            }
        }

        # Single write per OTU instead of ~10 writes.
        if (!no_log) {
            cat(log_buf, file = log_con)
        }
    }

    if (!no_log) {
        close(log_con)
    }

    # ---- 6. Build outputs (same as original) --------------------------------
    statistics_table$parent_id <- parent_id_v

    curation_table <- cbind(nOTUid = statistics_table$parent_id, otutable)
    statistics_table$curated <- "merged"
    curate_index <- rownames(statistics_table) == statistics_table$parent_id
    statistics_table$curated[curate_index] <- "parent"
    statistics_table <- transform(statistics_table,
                                  rank = ave(total,
                                             FUN = function(x)
                                                 rank(-x, ties.method = "first")))

    curation_table <- as.data.frame(curation_table %>%
                                    group_by(nOTUid) %>%
                                    summarise_all(list(sum)))
    rownames(curation_table) <- as.character(curation_table$nOTUid)
    curation_table <- curation_table[, -1]

    curated_otus    <- names(table(statistics_table$parent_id))
    curated_count   <- length(curated_otus)
    discarded_otus  <- setdiff(rownames(statistics_table), curated_otus)
    discarded_count <- length(discarded_otus)

    end.time   <- Sys.time()
    time.taken <- end.time - start.time

    list(curated_table  = curation_table,
         curated_count  = curated_count,
         curated_otus   = curated_otus,
         discarded_count = discarded_count,
         discarded_otus  = discarded_otus,
         runtime         = time.taken,
         minimum_match   = minimum_match,
         minimum_relative_cooccurence = minimum_relative_cooccurence,
         otu_map         = statistics_table,
         original_table  = otutable)
}

library(compiler)
lulu <- cmpfun(lulu, options = NULL)


# ---- CLI entry point (unchanged) -------------------------------------------
#args=c("lulu_match_list.txt","OTU.txt","logs")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) { stop("Not enough commandline args!") }

matchL <- args[1]
otuF   <- args[2]
logD   <- args[3]

info <- file.info(matchL)
if (info$size == 0) { q("no") }

matchs <- read.table(matchL, header = FALSE, as.is = TRUE)
otuM   <- read.table(otuF,   header = TRUE,  as.is = TRUE, row.names = 1)
if (dim(otuM)[1] <= 1 || dim(otuM)[2] <= 1) { q("no") }

lulu <- lulu(otuM, matchs)
write.table(lulu$curated_table, quote = FALSE, sep = "\t",
            file = otuF, col.names = NA)
write.table(lulu$discarded_otus, quote = FALSE, sep = "\t",
            file = paste0(matchL, ".rm"),
            col.names = FALSE, row.names = FALSE)
save(lulu, file = paste0(logD, "LULU.Rdata"))
