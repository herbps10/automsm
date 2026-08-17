.onLoad <- function(libname, pkgname) {
  ns <- asNamespace("posterior")
  for(g in c("as_draws", "as_draws_array", "as_draws_matrix", "as_draws_df",
             "as_draws_list", "as_draws_rvars")) {
    registerS3method(g, "automsm", get(paste0(g, ".automsm")), envir = ns)
  }
  invisible(NULL)
}

