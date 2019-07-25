library(data.table)
dat <- fread("./tests/testthat/testdata/cal_type_1.csv")

dat <- dat[,c(1,2)]
colnames(dat) <- c("true_methylation", "CpG1")

dat[,true_methylation:=factor(true_methylation)][,CpG1:=as.numeric(gsub(",", ".", CpG1))]

df_agg <- dat[,.(mean(CpG1)),by=true_methylation]
colnames(df_agg) <- c("true_methylation", "CpG")
df_agg[,true_methylation:=as.numeric(as.character(gsub(",", ".", true_methylation)))]


hyperbolic_equation <- function(x, b, a, d){
  # old equation (16.01.2019) with min = 0, max = 100
  #return((((y1 * b) - y0) * x + 100 * y0) / ((b * x) - x + 100))
  # new equation (17.01.2019) with data-dependent min and max
  # return((((b * y1) - y0) * (x - m0) + (m1 - m0) * y0) / ((b - 1) * (x - m0) + (m1 - m0)))
  # new equation (24.06.2019) without min and max
  return(((a * x) + b) / (x + d))
}

true_levels <- df_agg[,get("true_methylation")]

fn <- function(bias, a, d){
  fitted_vals <- hyperbolic_equation(true_levels, b = bias, a = a, d = d)
  # optimize biasfactor with minimizing sum of squares error
  return(sum(I(df_agg[,get("CpG")] - fitted_vals)^2))
}

# optimization function of built in R -> based on Nelder-Mead
# by default, optim performs minimization
# bias_factor <- optim(1, fn, method = "Nelder-Mead")$par
b <- stats::optim(c(1, 1, 1), fn, lower = 0, upper = 50)$par # due to error with Nelder-Mead

model <- nls(CpG ~ hyperbolic_equation(x = true_methylation, b, a, d),
             data = df_agg,
             start = c(b = 0, a = 1, d = 1))

true_levels <- df_agg[,get("true_methylation")]
