library(ape)

args <- commandArgs(trailingOnly = TRUE)

PROJECT_ROOT <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = FALSE)
}

if (basename(PROJECT_ROOT) %in% c("AA", "NT", "common")) {
  PROJECT_ROOT <- dirname(PROJECT_ROOT)
}

SEQTYPES <- if (length(args) >= 2L) {
  strsplit(args[[2]], ",", fixed = TRUE)[[1]]
} else {
  c("AA", "NT")
}

OVERWRITE <- (
  length(args) >= 3L
  && tolower(args[[3]]) %in% c("true", "1", "yes")
)


to_repository_relative <- function(path) {
  normalized_path <- normalizePath(
    path,
    mustWork = FALSE
  )

  repository_prefix <- paste0(
    PROJECT_ROOT,
    .Platform$file.sep
  )

  if (startsWith(normalized_path, repository_prefix)) {
    return(substr(
      normalized_path,
      nchar(repository_prefix) + 1L,
      nchar(normalized_path)
    ))
  }

  normalized_path
}


read_distance_matrix <- function(path) {
  x <- read.csv(
    path,
    row.names = 1,
    check.names = FALSE
  )

  D <- as.matrix(x)
  storage.mode(D) <- "double"

  if (nrow(D) != ncol(D)) {
    stop("Matrix is not square: ", path)
  }

  if (!identical(rownames(D), colnames(D))) {
    stop(
      "Row and column labels differ: ",
      path
    )
  }

  if (anyDuplicated(rownames(D)) > 0) {
    stop("Duplicate taxon labels: ", path)
  }

  if (any(!is.finite(D))) {
    stop("Non-finite distance: ", path)
  }

  if (min(D) < -1e-8) {
    stop(
      "Negative distance: ",
      path,
      "; min=",
      min(D)
    )
  }

  asymmetry <- max(abs(D - t(D)))

  if (asymmetry > 1e-8) {
    stop(
      "Asymmetric matrix: ",
      path,
      "; max asymmetry=",
      asymmetry
    )
  }

  D <- (D + t(D)) / 2
  D[D < 0] <- 0
  diag(D) <- 0

  D
}


run_nj_seqtype <- function(
  seqtype,
  selected_tags = NULL
) {
  input_root <- file.path(
    PROJECT_ROOT,
    seqtype,
    "data",
    "nj_distances"
  )

  output_root <- file.path(
    PROJECT_ROOT,
    seqtype,
    "results_nj"
  )

  input_paths <- list.files(
    input_root,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  input_paths <- input_paths[
    basename(input_paths) !=
      "export_manifest.csv"
  ]

  if (length(input_paths) == 0L) {
    stop("No distance matrices were found under: ", input_root)
  }

  records <- list()
  record_index <- 1L

  for (input_path in input_paths) {
    variant <- basename(dirname(input_path))

    tag <- tools::file_path_sans_ext(
      basename(input_path)
    )

    if (
      !is.null(selected_tags)
      && !tag %in% selected_tags
    ) {
      next
    }

    output_path <- file.path(
      output_root,
      variant,
      paste0(tag, ".nwk")
    )

    if (
      file.exists(output_path)
      && !OVERWRITE
    ) {
      records[[record_index]] <- data.frame(
        seqtype = seqtype,
        tag = tag,
        variant = variant,
        status = "skipped_existing",
        input_path = to_repository_relative(input_path),
        output_path = to_repository_relative(output_path),
        elapsed_sec = 0,
        message = "",
        stringsAsFactors = FALSE
      )

      record_index <- record_index + 1L
      next
    }

    start_time <- proc.time()[["elapsed"]]

    record <- tryCatch({
      D <- read_distance_matrix(input_path)

      tree <- ape::nj(
        as.dist(D)
      )

      if (
        !setequal(
          tree$tip.label,
          rownames(D)
        )
      ) {
        stop(
          "NJ output labels differ from input"
        )
      }

      dir.create(
        dirname(output_path),
        recursive = TRUE,
        showWarnings = FALSE
      )

      temporary_path <- paste0(
        output_path,
        ".tmp"
      )

      ape::write.tree(
        tree,
        file = temporary_path
      )

      if (
        !file.rename(
          temporary_path,
          output_path
        )
      ) {
        stop(
          "Could not write: ",
          output_path
        )
      }

      data.frame(
        seqtype = seqtype,
        tag = tag,
        variant = variant,
        status = "ok",
        input_path = to_repository_relative(input_path),
        output_path = to_repository_relative(output_path),
        elapsed_sec =
          proc.time()[["elapsed"]] -
          start_time,
        message = "",
        stringsAsFactors = FALSE
      )

    }, error = function(error) {
      data.frame(
        seqtype = seqtype,
        tag = tag,
        variant = variant,
        status = "error",
        input_path = to_repository_relative(input_path),
        output_path = to_repository_relative(output_path),
        elapsed_sec =
          proc.time()[["elapsed"]] -
          start_time,
        message = conditionMessage(error),
        stringsAsFactors = FALSE
      )
    })

    records[[record_index]] <- record
    record_index <- record_index + 1L

    message(
      "[",
      seqtype,
      "] ",
      variant,
      " ",
      tag,
      " : ",
      record$status
    )
  }

  log_data <- do.call(
    rbind,
    records
  )

  dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    log_data,
    file.path(
      output_root,
      "nj_run_log.csv"
    ),
    row.names = FALSE
  )

  print(
    table(
      log_data$variant,
      log_data$status
    )
  )

  log_data
}


invalid_seqtypes <- setdiff(SEQTYPES, c("AA", "NT"))

if (length(invalid_seqtypes) > 0L) {
  stop(
    "Unknown sequence type(s): ",
    paste(invalid_seqtypes, collapse = ", ")
  )
}

invisible(
  lapply(
    SEQTYPES,
    run_nj_seqtype
  )
)
