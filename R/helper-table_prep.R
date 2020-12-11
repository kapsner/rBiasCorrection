table_prep <- function(table) {
  vec <- colnames(table)[unlist(table[, lapply(.SD, is.numeric)])]
  ret <- table[
    ,
    lapply(
      X = .SD,
      FUN = function(x) {
        round(x, 4)
      }
    ),
    .SDcols = vec
  ]
  return(ret)
}
