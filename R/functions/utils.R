## Utility functions for generating model matrices with smooths for TMB
# Taken from sdmTMB utils.R and smoothers.R

safe_deparse <- function (x, collapse = " ") {
  paste(deparse(x, 500L), collapse = collapse)
}

barnames <- function (bars) {
  vapply(bars, function(x) safe_deparse(x[[3]]), "")
}

check_valid_factor_levels <- function(x, .name = "") {
  assertthat::assert_that(is.factor(x),
              msg = sprintf("Random effect group column `%s` is not a factor.", .name))
  lev <- sort(levels(x))
  uni <- sort(unique(as.character(x)))
  assert_that(identical(lev, uni),
              msg = sprintf(
                "Random effect group column `%s` has extra factor levels. Please remove them.", .name))
}

rm_wsp <- function (x) {
  out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
  dim(out) <- dim(x)
  out
}

all_terms <- function (x) {
  if (!length(x)) {
    return(character(0))
  }
  if (!inherits(x, "terms")) {
    x <- terms(stats::as.formula(x))
  }
  rm_wsp(attr(x, "term.labels"))
}

get_smooth_terms <- function(terms) {
  x1 <- grep("s\\(", terms)
  x2 <- grep("t2\\(", terms)
  if (length(x2) > 0L)
    stop("t2() smoothers are not yet supported due to issues with prediction on newdata.", call. = FALSE)
  x1
}

s2rPred <- function(sm, re, data) {
  ## Function to aid prediction from smooths represented as type==2
  ## random effects. re must be the result of smooth2random(sm,...,type=2).
  X <- mgcv::PredictMat(sm, data) ## get prediction matrix for new data
  ## transform to r.e. parameterization
  if (!is.null(re$trans.U)) {
    X <- X %*% re$trans.U
  }
  X <- t(t(X) * re$trans.D)
  ## re-order columns according to random effect re-ordering...
  X[, re$rind] <- X[, re$pen.ind != 0]
  ## re-order penalization index in same way
  pen.ind <- re$pen.ind
  pen.ind[re$rind] <- pen.ind[pen.ind > 0]
  ## start return object...
  r <- list(rand = list(), Xf = X[, which(re$pen.ind == 0), drop = FALSE])
  for (i in seq_along(re$rand)) { ## loop over random effect matrices
    r$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
    attr(r$rand[[i]], "s.label") <- attr(re$rand[[i]], "s.label")
  }
  names(r$rand) <- names(re$rand)
  r
}

parse_smoothers <- function(formula, data, newdata = NULL) {
  terms <- all_terms(formula)
  if (!is.null(newdata)) {
    if (any(grepl("t2\\(", terms))) stop("Prediction on newdata with t2() still has issues.", call. = FALSE)
  }
  smooth_i <- get_smooth_terms(terms)
  basis <- list()
  Zs <- list()
  Xs <- list()
  if (length(smooth_i) > 0) {
    has_smooths <- TRUE
    smterms <- terms[smooth_i]
    ns <- 0
    ns_Xf <- 0
    for (i in seq_along(smterms)) {
      obj <- eval(str2expression(smterms[i]))
      basis[[i]] <- mgcv::smoothCon(
        object = obj, data = data,
        knots = NULL, absorb.cons = TRUE,
        diagonal.penalty = TRUE
      )
      for (j in seq_along(basis[[i]])) { # elements > 1 with `by` terms
        ns_Xf <- ns_Xf + 1
        rasm <- mgcv::smooth2random(basis[[i]][[j]], names(data), type = 2)
        if (!is.null(newdata)) {
          rasm <- s2rPred(basis[[i]][[j]], rasm, data = newdata)
        }
        for (k in seq_along(rasm$rand)) { # elements > 1 with if s(x, y) or t2()
          ns <- ns + 1
          Zs[[ns]] <- rasm$rand[[k]]
        }
        Xs[[ns_Xf]] <- rasm$Xf
      }
    }
    sm_dims <- unlist(lapply(Zs, ncol))
    Xs <- do.call(cbind, Xs) # combine 'em all into one design matrix
    b_smooth_start <- c(0, cumsum(sm_dims)[-length(sm_dims)])
  } else {
    has_smooths <- FALSE
    sm_dims <- 0L
    b_smooth_start <- 0L
    Xs <- matrix(nrow = 0L, ncol = 0L)
  }
  list(Xs = Xs, Zs = Zs, has_smooths = has_smooths,
       sm_dims = sm_dims, b_smooth_start = b_smooth_start)
}

remove_s_and_t2 <- function(formula) {
  terms <- stats::terms(formula)
  terms_labels <- attr(terms, "term.labels")
  drop <- grep("s\\(", terms_labels)
  dropt2 <- grep("t2\\(", terms_labels)
  tdrop <- terms_labels[union(drop, dropt2)]
  for (i in seq_along(tdrop)) {
    formula <- stats::update(formula, paste("~ . -", tdrop[i]))
  }
  formula
}