# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2022 Lorenz Kapsner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
