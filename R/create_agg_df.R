# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2025 Lorenz Kapsner
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

# create aggregated datatable for calibration data
create_agg_df <- function(datatable,
                          index) {
  df <- datatable[, c("true_methylation", index), with = FALSE]
  colnames(df)[2] <- "CpG"

  df_out <- df[, ("true_methylation") := as.numeric(
    as.character(
      get("true_methylation")
    )
  )]

  df_out <- df[
    ,
    list(
      "CpG" = mean(get("CpG"), na.rm = TRUE),
      "sd" = stats::sd(get("CpG"), na.rm = TRUE)
    ),
    by = "true_methylation"
  ]

  if (df_out[is.na(get("sd")), .N] == nrow(df_out)) {
    df_out[, ("sd") := NULL]
  }

  # calculate relative error
  df_out[, ("CpG_true_diff") := abs(get("CpG") - get("true_methylation"))]
  df_out[, ("relative_error") := ifelse(
    get("true_methylation") != 0,
    (get("CpG_true_diff") / get("true_methylation")) * 100,
    NA
  )]

  return(df_out[!is.na(get("CpG")), ])
}

# create aggregated datatable for experimental data
create_agg_df_exp <- function(datatable,
                              index,
                              type) {
  if (type == 1) {
    df <- datatable[, c("sample_id", index), with = FALSE]
    colnames(df)[2] <- "CpG"

    df_out <- df[
      ,
      list(
        "CpG" = mean(get("CpG"), na.rm = TRUE),
        "sd" = stats::sd(get("CpG"), na.rm = TRUE)
      ),
      by = "sample_id"
    ]
  } else if (type == 2) {
    df <- datatable[, c("locus_id", index), with = FALSE]
    colnames(df)[2] <- "CpG"
    df_out <- df[
      ,
      list(
        "CpG" = mean(get("CpG"), na.rm = TRUE),
        "sd" = stats::sd(get("CpG"), na.rm = TRUE)
      ),
      by = "locus_id"
    ]
  }

  if (df_out[is.na(get("sd")), .N] == nrow(df_out)) {
    df_out[, ("sd") := NULL]
  }

  return(df_out[!is.na(get("CpG")), ])
}
