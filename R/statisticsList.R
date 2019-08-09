
#' @title statisticsList_ helper function
#'
#' @description Internal function that converts the results_list (output of \code{regressionUtility_()})
#'   into a data.table object.
#'
#' @param resultlist A list object. The results_list output of \code{regressionUtility_()}.
#' @inheritParams BiasCorrection
#'
#' @export
#'
statisticsList_ <- function(resultlist, minmax = FALSE){
  if (isFALSE(minmax)){
    outdat <- data.table::data.table("Name" = names(resultlist),
                                      "relative_error" = NA,
                                      "SSE_hyperbolic" = NA,
                                      "R2_hyperbolic" = NA,
                                      "a" = NA,
                                      "b" = NA,
                                      "d_h" = NA,
                                      "###" = NA,
                                      "SSE_cubic" = NA,
                                      "R2_cubic" = NA,
                                      "ax3" = NA,
                                      "bx2" = NA,
                                      "cx" = NA,
                                      "d" = NA)

    outdat[, ("Name") := names(resultlist)]

    vec <- names(outdat)[-1]
    outdat[,(vec) := lapply(.SD, function(x){as.numeric(as.character(x))}), .SDcols = vec]

    for (i in names(resultlist)){
      out_names <- c("relative_error",
                     "SSE_hyperbolic",
                     "R2_hyperbolic",
                     "a",
                     "b",
                     "d_h",
                     "SSE_cubic",
                     "R2_cubic",
                     "ax3",
                     "bx2",
                     "cx",
                     "d")
      out_list <- list(resultlist[[i]][["relative_error"]],
                       resultlist[[i]][["SSE_hyper"]],
                       resultlist[[i]][["Coef_hyper"]][["R2"]],
                       resultlist[[i]][["Coef_hyper"]][["a"]],
                       resultlist[[i]][["Coef_hyper"]][["b"]],
                       resultlist[[i]][["Coef_hyper"]][["d"]],
                       resultlist[[i]][["SSE_cubic"]],
                       resultlist[[i]][["Coef_cubic"]][["R2"]],
                       resultlist[[i]][["Coef_cubic"]][["ax3"]],
                       resultlist[[i]][["Coef_cubic"]][["bx2"]],
                       resultlist[[i]][["Coef_cubic"]][["cx"]],
                       resultlist[[i]][["Coef_cubic"]][["d"]])
      outdat[get("Name") == i, (out_names) := out_list]
    }
  } else if (isTRUE(minmax)){
    outdat <- data.table::data.table("Name" = names(resultlist),
                                      "relative_error" = NA,
                                      "SSE_hyperbolic" = NA,
                                      "R2_hyperbolic" = NA,
                                      "b" = NA,
                                      "y0" = NA,
                                      "y1" = NA,
                                      "m0" = NA,
                                      "m1" = NA,
                                      "###" = NA,
                                      "SSE_cubic" = NA,
                                      "R2_cubic" = NA,
                                      "a_c" = NA,
                                      "b_c" = NA,
                                     "y0_c" = NA,
                                     "y1_c" = NA,
                                     "m0_c" = NA,
                                     "m1_c" = NA)

    outdat[, ("Name") := names(resultlist)]

    vec <- names(outdat)[-1]
    outdat[,(vec) := lapply(.SD, function(x){as.numeric(as.character(x))}), .SDcols = vec]

    for (i in names(resultlist)){
      out_names <- c("relative_error",
                     "SSE_hyperbolic",
                     "R2_hyperbolic",
                     "b",
                     "y0",
                     "y1",
                     "m0",
                     "m1",
                     "SSE_cubic",
                     "R2_cubic",
                     "a_c",
                     "b_c",
                     "y0_c",
                     "y1_c",
                     "m0_c",
                     "m1_c")
      out_list <- list(resultlist[[i]][["relative_error"]],
                       resultlist[[i]][["SSE_hyper"]],
                       resultlist[[i]][["Coef_hyper"]][["R2"]],
                       resultlist[[i]][["Coef_hyper"]][["b"]],
                       resultlist[[i]][["Coef_hyper"]][["y0"]],
                       resultlist[[i]][["Coef_hyper"]][["y1"]],
                       resultlist[[i]][["Coef_hyper"]][["m0"]],
                       resultlist[[i]][["Coef_hyper"]][["m1"]],
                       resultlist[[i]][["SSE_cubic"]],
                       resultlist[[i]][["Coef_cubic"]][["R2"]],
                       resultlist[[i]][["Coef_cubic"]][["a"]],
                       resultlist[[i]][["Coef_cubic"]][["b"]],
                       resultlist[[i]][["Coef_cubic"]][["y0"]],
                       resultlist[[i]][["Coef_cubic"]][["y1"]],
                       resultlist[[i]][["Coef_cubic"]][["m0"]],
                       resultlist[[i]][["Coef_cubic"]][["m1"]])
      outdat[get("Name") == i, (out_names) := out_list]
    }
  }
  # # mark the better model: 1 = cubic, 0 = hyperbolic
  # outdat[,("better_model") := ifelse(get("SSE_cubic") <= get("SSE_hyperbolic"), 1, 0)]
  return(outdat)
}
