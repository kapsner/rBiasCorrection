
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
    dt_list <- data.table::data.table("Name" = names(resultlist),
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

    dt_list[, ("Name") := names(resultlist)]

    vec <- names(dt_list)[-1]
    dt_list[,(vec) := lapply(.SD, function(x){as.numeric(as.character(x))}), .SDcols = vec]

    for (i in names(resultlist)){
      dt_list[get("Name") == i, ("relative_error") := resultlist[[i]][["relative_error"]]
              ][
                get("Name") == i, ("SSE_hyperbolic") := resultlist[[i]][["SSE_hyper"]]
                ][
                  get("Name") == i, ("R2_hyperbolic") := resultlist[[i]][["Coef_hyper"]][["R2"]]
                  ][
                    get("Name") == i, ("a") := resultlist[[i]][["Coef_hyper"]][["a"]]
                    ][
                      get("Name") == i, ("b") := resultlist[[i]][["Coef_hyper"]][["b"]]
                      ][
                        get("Name") == i, ("d_h") := resultlist[[i]][["Coef_hyper"]][["d"]]
                        ][
                          get("Name") == i, ("SSE_cubic") := resultlist[[i]][["SSE_cubic"]]
                          ][
                            get("Name") == i, ("R2_cubic") := resultlist[[i]][["Coef_cubic"]][["R2"]]
                            ][
                              get("Name") == i, ("ax3") := resultlist[[i]][["Coef_cubic"]][["ax3"]]
                              ][
                                get("Name") == i, ("bx2") := resultlist[[i]][["Coef_cubic"]][["bx2"]]
                                ][
                                  get("Name") == i, ("cx") := resultlist[[i]][["Coef_cubic"]][["cx"]]
                                  ][
                                    get("Name") == i, ("d") := resultlist[[i]][["Coef_cubic"]][["d"]]
                                    ]
    }
  } else if (isTRUE(minmax)){
    dt_list <- data.table::data.table("Name" = names(resultlist),
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
                                      "ax3" = NA,
                                      "bx2" = NA,
                                      "cx" = NA,
                                      "d" = NA)

    dt_list[, ("Name") := names(resultlist)]

    vec <- names(dt_list)[-1]
    dt_list[,(vec) := lapply(.SD, function(x){as.numeric(as.character(x))}), .SDcols = vec]

    for (i in names(resultlist)){
      dt_list[get("Name") == i, ("relative_error") := resultlist[[i]][["relative_error"]]
              ][
                get("Name") == i, ("SSE_hyperbolic") := resultlist[[i]][["SSE_hyper"]]
                ][
                  get("Name") == i, ("R2_hyperbolic") := resultlist[[i]][["Coef_hyper"]][["R2"]]
                  ][
                    get("Name") == i, ("b") := resultlist[[i]][["Coef_hyper"]][["b"]]
                    ][
                      get("Name") == i, ("y0") := resultlist[[i]][["Coef_hyper"]][["y0"]]
                      ][
                        get("Name") == i, ("y1") := resultlist[[i]][["Coef_hyper"]][["y1"]]
                        ][
                          get("Name") == i, ("m0") := resultlist[[i]][["Coef_hyper"]][["m0"]]
                          ][
                            get("Name") == i, ("m1") := resultlist[[i]][["Coef_hyper"]][["m1"]]
                            ][
                              get("Name") == i, ("SSE_cubic") := resultlist[[i]][["SSE_cubic"]]
                              ][
                                get("Name") == i, ("R2_cubic") := resultlist[[i]][["Coef_cubic"]][["R2"]]
                                ][
                                  get("Name") == i, ("ax3") := resultlist[[i]][["Coef_cubic"]][["ax3"]]
                                  ][
                                    get("Name") == i, ("bx2") := resultlist[[i]][["Coef_cubic"]][["bx2"]]
                                    ][
                                      get("Name") == i, ("cx") := resultlist[[i]][["Coef_cubic"]][["cx"]]
                                      ][
                                        get("Name") == i, ("d") := resultlist[[i]][["Coef_cubic"]][["d"]]
                                        ]
    }
  }

  # mark the better model: 1 = cubic, 0 = hyperbolic
  dt_list[,("better_model") := ifelse(get("SSE_cubic") <= get("SSE_hyperbolic"), 1, 0)]

  return(dt_list)
}
