.cdt_require_data_frame <- function(data) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  invisible(TRUE)
}

.cdt_require_columns <- function(data, cols, arg = "data") {
  cols <- cols[!is.na(cols) & nzchar(cols)]
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop("Column(s) not found in '", arg, "': ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.cdt_make_pairs <- function(D, components, interactions) {
  C <- length(components)
  if (C < 2) return(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("pair_a", "pair_b"))))

  pairs <- matrix(integer(0), ncol = 2)
  if (is.null(interactions)) {
    for (i in seq_len(C - 1)) {
      for (j in (i + 1):C) {
        if (any(D[, i] > 0 & D[, j] > 0)) pairs <- rbind(pairs, c(i, j))
      }
    }
  } else {
    if (!is.list(interactions)) stop("'interactions' must be NULL or a list of length-2 character vectors.", call. = FALSE)
    for (pr in interactions) {
      if (!is.character(pr) || length(pr) != 2L) {
        stop("Each interaction pair must be a length-2 character vector.", call. = FALSE)
      }
      i <- match(pr[1], components)
      j <- match(pr[2], components)
      if (anyNA(c(i, j))) stop("Interaction pair names must match 'components'.", call. = FALSE)
      if (i == j) stop("A component cannot interact with itself: ", pr[1], call. = FALSE)
      pairs <- rbind(pairs, sort(c(i, j)))
    }
    if (nrow(pairs)) pairs <- unique(pairs)
  }
  colnames(pairs) <- c("pair_a", "pair_b")
  pairs
}

.cdt_reference_arms <- function(data, study_idx, studies, D, ref = NULL) {
  N <- nrow(D)
  S <- length(studies)
  if (!is.null(ref) && ref %in% names(data)) {
    raw <- data[[ref]]
    is_ref <- as.integer(as.logical(raw))
    if (anyNA(is_ref)) stop("Reference column '", ref, "' contains values that cannot be interpreted as TRUE/FALSE.", call. = FALSE)
    for (s in seq_len(S)) {
      n_ref <- sum(is_ref[study_idx == s])
      if (n_ref != 1L) {
        stop("Reference column must flag exactly one arm in study '", studies[s], "' (found ", n_ref, ").", call. = FALSE)
      }
    }
    return(is_ref)
  }

  is_ref <- integer(N)
  tot <- rowSums(D)
  for (s in seq_len(S)) {
    rows <- which(study_idx == s)
    zero <- rows[tot[rows] == 0]
    base <- if (length(zero)) zero[1] else rows[which.min(tot[rows])]
    if (!length(zero)) {
      warning("Study '", studies[s], "' has no all-zero arm; using its lowest-total-dose arm as baseline.", call. = FALSE)
    }
    is_ref[base] <- 1L
  }
  is_ref
}

#' Assemble a component dose-response network
#'
#' Turns an arm-level data frame into the design that [cdt_fit()] consumes. Each
#' row is one treatment arm. Every component named in `components` must have a
#' dose column, with zero meaning the component is absent from that arm.
#'
#' @param data A data frame with one row per study arm.
#' @param study Name of the column identifying the study.
#' @param components Character vector of component dose column names. A value of
#'   zero marks the component absent.
#' @param outcome Either `"continuous"` or `"binary"`.
#' @param y For a continuous outcome, the column of arm mean changes.
#' @param se For a continuous outcome, the column of standard errors. Supply this
#'   or both `sd` and `n`.
#' @param sd,n For a continuous outcome, columns of arm standard deviation and
#'   sample size, used as `se = sd / sqrt(n)` when `se` is absent.
#' @param events,n_binary For a binary outcome, the event-count and sample-size
#'   columns.
#' @param ref Optional column flagging the study-baseline arm. It may be logical
#'   or 0/1 and must flag exactly one arm per study. When absent, the all-zero
#'   dose arm is used. A study with no all-zero arm takes its lowest-total-dose
#'   arm as the within-study baseline.
#' @param dstar Optional named numeric vector of dose-normalisation references
#'   used by the bilinear surface. Defaults to the maximum observed dose of each
#'   component.
#' @param interactions Pairs of components to give an interaction term. Either
#'   `NULL`, the default, which uses all component pairs that co-occur in at
#'   least one arm, or a list of length-2 character vectors naming the pairs.
#'
#' @return An object of class `cdt_data`: a list holding the partial Stan data,
#'   component labels, study labels, interaction pairs, and dose normalisers.
#' @export
cdt_data <- function(data, study, components,
                     outcome = c("continuous", "binary"),
                     y = NULL, se = NULL, sd = NULL, n = NULL,
                     events = NULL, n_binary = NULL,
                     ref = NULL, dstar = NULL, interactions = NULL) {
  outcome <- match.arg(outcome)
  .cdt_require_data_frame(data)

  if (!is.character(study) || length(study) != 1L) stop("'study' must be a single column name.", call. = FALSE)
  if (!is.character(components) || !length(components)) stop("'components' must be a non-empty character vector.", call. = FALSE)
  if (anyDuplicated(components)) stop("'components' contains duplicate names.", call. = FALSE)
  .cdt_require_columns(data, c(study, components), "data")

  D <- as.matrix(data[, components, drop = FALSE])
  storage.mode(D) <- "double"
  if (anyNA(D)) stop("Component dose columns must not contain NA; use 0 for absent components.", call. = FALSE)
  if (any(D < 0)) stop("Component doses must be non-negative.", call. = FALSE)
  C <- ncol(D)
  N <- nrow(D)
  if (N < 2L) stop("At least two arms are needed.", call. = FALSE)

  sf <- factor(data[[study]])
  study_idx <- as.integer(sf)
  studies <- levels(sf)
  S <- length(studies)

  yv <- rep(0, N)
  sev <- rep(1, N)
  rv <- rep(0L, N)
  nv <- rep(0L, N)

  if (outcome == "continuous") {
    if (is.null(y)) stop("Provide 'y' for a continuous outcome.", call. = FALSE)
    .cdt_require_columns(data, y, "data")
    yv <- as.double(data[[y]])

    if (!is.null(se) && se %in% names(data)) {
      sev <- as.double(data[[se]])
    } else if (!is.null(sd) && !is.null(n) && all(c(sd, n) %in% names(data))) {
      nn <- as.double(data[[n]])
      if (anyNA(nn) || any(nn <= 0)) stop("Sample sizes in 'n' must be positive.", call. = FALSE)
      sev <- as.double(data[[sd]]) / sqrt(nn)
    } else {
      stop("Provide 'se', or both 'sd' and 'n', for a continuous outcome.", call. = FALSE)
    }
    if (anyNA(yv) || anyNA(sev)) stop("Continuous outcome columns contain NA.", call. = FALSE)
    if (any(sev <= 0)) stop("Standard errors must be positive.", call. = FALSE)
  } else {
    if (is.null(events) || is.null(n_binary)) stop("Provide 'events' and 'n_binary' for a binary outcome.", call. = FALSE)
    .cdt_require_columns(data, c(events, n_binary), "data")
    rv <- as.integer(data[[events]])
    nv <- as.integer(data[[n_binary]])
    if (anyNA(rv) || anyNA(nv)) stop("Binary outcome columns contain NA.", call. = FALSE)
    if (any(rv < 0L) || any(nv < rv)) stop("Need 0 <= events <= n_binary for every arm.", call. = FALSE)
    if (any(nv <= 0L)) stop("Binary sample sizes must be positive.", call. = FALSE)
  }

  is_ref <- .cdt_reference_arms(data, study_idx, studies, D, ref = ref)

  if (is.null(dstar)) {
    dstar <- apply(D, 2, max)
    dstar[dstar == 0] <- 1
  } else {
    if (is.null(names(dstar))) stop("'dstar' must be a named numeric vector.", call. = FALSE)
    dstar <- as.double(dstar[components])
    if (anyNA(dstar)) stop("'dstar' must name every component.", call. = FALSE)
    if (any(dstar <= 0)) stop("Every value in 'dstar' must be positive.", call. = FALSE)
  }
  names(dstar) <- components

  pairs <- .cdt_make_pairs(D, components, interactions)
  P <- nrow(pairs)
  pa <- if (P) as.integer(pairs[, 1]) else integer(0)
  pb <- if (P) as.integer(pairs[, 2]) else integer(0)

  standata <- list(
    N = as.integer(N), S = as.integer(S), C = as.integer(C), P = as.integer(P),
    study = as.integer(study_idx), D = D, is_ref = as.integer(is_ref),
    pair_a = as.array(pa), pair_b = as.array(pb),
    dstar = as.array(unname(dstar)),
    outcome = if (outcome == "continuous") 1L else 2L,
    y = yv, se = sev, r = rv, nn = nv
  )

  structure(list(
    standata = standata,
    components = components,
    studies = studies,
    outcome = outcome,
    pairs = pairs,
    dstar = dstar,
    n_arms = N,
    n_studies = S
  ), class = "cdt_data")
}

#' @export
print.cdt_data <- function(x, ...) {
  cat("<cdt_data>\n")
  cat(sprintf("  arms:       %d across %d studies\n", x$n_arms, x$n_studies))
  cat(sprintf("  components: %s\n", paste(x$components, collapse = ", ")))
  cat(sprintf("  outcome:    %s\n", x$outcome))
  np <- nrow(x$pairs)
  if (np) {
    labs <- apply(x$pairs, 1, function(p) paste(x$components[p], collapse = " x "))
    cat(sprintf("  interaction pairs: %s\n", paste(labs, collapse = "; ")))
  } else {
    cat("  interaction pairs: none identifiable from co-occurring components\n")
  }
  invisible(x)
}
