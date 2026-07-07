#!/usr/bin/env Rscript

# Hardened R package installer for LotuS helper packages.
# - avoids source() from remote installer scripts
# - checks install results explicitly
# - uses HTTPS CRAN mirror
# - exits non-zero if required packages could not be installed
# - supports simple environment controls for non-interactive / CI installs

required_cran <- c("dplyr", "ape")
required_bioc <- c("phyloseq", "dada2")
required_all  <- c(required_cran, required_bioc)

cran_repo <- Sys.getenv("LOTUS_CRAN_REPO", unset = "https://cloud.r-project.org")
bioc_version <- Sys.getenv("LOTUS_BIOC_VERSION", unset = NA_character_)
stop_on_fail <- tolower(Sys.getenv("LOTUS_R_INSTALL_STRICT", unset = "true")) %in% c("1", "true", "yes", "y")

options(
  repos = c(CRAN = cran_repo),
  timeout = max(300, getOption("timeout", 60)),
  install.packages.check.source = "no"
)

message("LotuS R dependency installer")
message("R version: ", paste(R.version$major, R.version$minor, sep = "."))
message("CRAN repo: ", cran_repo)

pkg_available <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

load_pkg <- function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
}

install_cran_pkg <- function(pkg) {
  if (pkg_available(pkg)) {
    message(pkg, " is already installed.")
    return(TRUE)
  }

  message("Installing CRAN package: ", pkg)
  ok <- tryCatch({
    install.packages(pkg, repos = cran_repo, dependencies = TRUE, quiet = FALSE)
    pkg_available(pkg)
  }, warning = function(w) {
    message("WARNING while installing ", pkg, ": ", conditionMessage(w))
    pkg_available(pkg)
  }, error = function(e) {
    message("ERROR while installing ", pkg, ": ", conditionMessage(e))
    FALSE
  })

  if (ok) {
    message("Installed CRAN package: ", pkg)
  } else {
    message("Failed to install CRAN package: ", pkg)
  }
  ok
}

ensure_biocmanager <- function() {
  if (pkg_available("BiocManager")) {
    return(TRUE)
  }

  message("Installing CRAN package: BiocManager")
  ok <- tryCatch({
    install.packages("BiocManager", repos = cran_repo, dependencies = TRUE, quiet = FALSE)
    pkg_available("BiocManager")
  }, warning = function(w) {
    message("WARNING while installing BiocManager: ", conditionMessage(w))
    pkg_available("BiocManager")
  }, error = function(e) {
    message("ERROR while installing BiocManager: ", conditionMessage(e))
    FALSE
  })

  ok
}

install_bioc_pkg <- function(pkg) {
  if (pkg_available(pkg)) {
    message(pkg, " is already installed.")
    return(TRUE)
  }

  if (!ensure_biocmanager()) {
    message("Cannot install ", pkg, " because BiocManager is unavailable.")
    return(FALSE)
  }

  message("Installing Bioconductor package: ", pkg)
  ok <- tryCatch({
    if (!is.na(bioc_version) && nzchar(bioc_version)) {
      BiocManager::install(pkg, version = bioc_version, ask = FALSE, update = FALSE)
    } else {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
    pkg_available(pkg)
  }, warning = function(w) {
    message("WARNING while installing ", pkg, ": ", conditionMessage(w))
    pkg_available(pkg)
  }, error = function(e) {
    message("ERROR while installing ", pkg, ": ", conditionMessage(e))
    FALSE
  })

  if (ok) {
    message("Installed Bioconductor package: ", pkg)
  } else {
    message("Failed to install Bioconductor package: ", pkg)
  }
  ok
}

# Install CRAN packages first; Bioconductor packages may depend on them.
cran_ok <- vapply(required_cran, install_cran_pkg, logical(1))
bioc_ok <- vapply(required_bioc, install_bioc_pkg, logical(1))

# Verify all required packages by namespace availability and loadability.
missing <- required_all[!vapply(required_all, pkg_available, logical(1))]
not_loadable <- required_all[!vapply(required_all, load_pkg, logical(1))]
failed <- unique(c(missing, not_loadable))

if (length(failed) > 0) {
  msg <- paste(
    "Package", paste(failed, collapse = ", "),
    "could not be installed. Please install it manually in your R environment."
  )
  message(msg)
  if (stop_on_fail) {
    quit(status = 1, save = "no")
  }
} else {
  message("All required R packages are installed and loadable: ", paste(required_all, collapse = ", "))
}

quit(status = 0, save = "no")
