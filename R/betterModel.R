betterModel <- function(statstable_pre, statstable_post_hyperbolic = NULL, statstable_post_cubic = NULL, selection_method = "SSE"){

  stopifnot(
    is.character(selection_method),
    selection_method %in% c("SSE", "RelError"),
    data.table::is.data.table(statstable_pre),
    data.table::is.data.table(statstable_post_hyperbolic) || is.null(statstable_post_hyperbolic),
    data.table::is.data.table(statstable_post_cubic) || is.null(statstable_post_cubic)
  )

  if (selection_method == "SSE"){
    outdat <- statstable_pre[,c("Name", "SSE_hyperbolic", "SSE_cubic"),with=F]
    # mark the better model: 1 = cubic, 0 = hyperbolic
    outdat[,("better_model") := ifelse(get("SSE_cubic") <= get("SSE_hyperbolic"), 1, 0)]
  } else if (selection_method == "RelError"){
    outdat <- statstable_post_hyperbolic[,c("Name", "relative_error"),with=F]
    colnames(outdat) <- c("Name", "relative_error_h")
    outdat <- merge(x = outdat,
                    y = statstable_post_cubic[,c("Name", "relative_error"), with=F],
                    by = "Name",
                    all = T,
                    suffixes = c("", "_c"))
    outdat[,("better_model") := ifelse(get("relative_error_c") <= get("relative_error_h"), 1, 0)]
  }
  return(outdat)
}
