#####################################################################################################
# Function borrowed from the R package 'ivx'                                                        #
#                                                                                                   #
# Reference: Yang, B., Long, W., Peng, L., & Cai, Z. (2020).                                        #
# Testing the predictability of US housing price index returns based on an IVX-AR model (IVX-AR).   #
# Journal of the American Statistical Association, 115(532), 1598-1619.                             #
#####################################################################################################

ivx_ar_fit <- function(y, x, horizon = 1, offset = NULL, ar = "auto", ar_max = 5, ar_ic = "bic",
                       ar_grid =  function(x) seq(x - 0.3, x + 0.3, by = 0.02), ...) 
{# begin-of-function
  
  mdl_ivx <- ivx_fit(y, x, horizon = horizon)
  
  if (ar == "auto") 
  {
    mdl_ar <- auto_ar(mdl_ivx$ols$residuals, d = 0, max.p = ar_max, ar_ic = ar_ic, ...)
  } 
  else if( ar == "forecast") 
  {
    requireNamespace("forecast", quietly = TRUE)
    mdl_ar <- forecast::auto.arima(mdl_ivx$ols$residuals, d = 0, max.p = ar_max, max.q = 0, ic = ar_ic, ...)
  } 
  else if (ar == 0)
  {
    message("Using `ivx` instead.")
    return(mdl_ivx)
  }
  else 
  {
    mdl_ar <- arima2(mdl_ivx$ols$residuals, order = c(ar, 0, 0), include.mean = FALSE, ...)
  }
  
  ar_coefs <- coefficients(mdl_ar)
  # in case the arima does not converge
  if (is_numeric0(ar_coefs))
  {
    return( list(
            coefficients = numeric(), residuals = y,
            fitted = 0 * y, df.residuals = length(y),
            ar_method = ar, ar_aic = ar_ic
                )
          )
   }
  
  res_ar   <- residuals(mdl_ar)
  q        <- length(ar_coefs)
  grid_seq <- sapply(ar_coefs, ar_grid)
  ngrid    <- nrow(grid_seq)
  
  res_ivx <- vector("list", length = ngrid)
  rse     <- vector("numeric", length = ngrid)
 
  for (i in 1:ngrid) 
  {
    y_adj <- tilt(y, grid_seq[i,], q)
    x_adj <- tilt(x, grid_seq[i,], q)
    res_ivx[[i]] <- ivx::ivx_fit(y_adj, x_adj, horizon = horizon)
    eps     <- y_adj - sum(x_adj * res_ivx[[i]]$coefficients)
    rse[i]  <- var(eps[!is.infinite(eps)])
  }
  
  rse_min <- which.min(rse)
  z <- res_ivx[[rse_min]]
  z$rse <- rse[rse_min]
  z$coefficients_ar <- grid_seq[rse_min, ]
  z$ar_method <- ar
  z$ar_ic <- ar_ic
  
  z$Wald_AR <- ac_test_wald(mdl_ivx$ols$residuals, q)
  z$q <- q
  z
  
}# end-of-function

####################################################
# Outputs
####################################################

print.ivx_ar <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{# begin-of-function
  
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")
  cat("\n\nLag Selection:\n",
      if(x$ar_method == "auto") paste0("Auto (", x$ar_ic, ")") else "Fixed",
      " with AR terms q = ", x$q, "\n\n", sep ="")
  res <- x$coefficients
  
  if (length(res)) 
  {
    cat("Coefficients:\n")
    print.default(format(res, digits = digits), print.gap = 2L, quote = FALSE)
  } 
  else 
  {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(x)
  
}# end-of-function


####################################################
# Summarizing IVX-AR Model Fits
####################################################

summary.ivx_ar <- function(object,  ...) 
{# begin-of-function

  z <- object
  
  if (is.null(z$terms))
    stop("invalid 'ivx' object: no 'terms' components")
  if (!inherits(object, "ivx_ar"))
    stop("calling summary.ivx(<fake-ivx-object>) ...")
  
  ans <- z[c("call", "terms")]
  
  ans$aliased <- is.na(z$coefficients)
  
  p_value_ivx                <- 1 - pchisq(z$Wald_Ind, 1)
  ans$coefficients           <- cbind(z$coefficients, z$Wald_Ind, p_value_ivx)
  dimnames(ans$coefficients) <- list(z$cnames, c("Estimate", "Wald Ind", "Pr(> chi)"))
  
  ans$vcov           <- z$vcov
  dimnames(ans$vcov) <- dimnames(ans$coefficients)[c(1, 1)]
  
  ans$delta           <- z$delta
  colnames(ans$delta) <- z$cnames
  
  ans$residuals <- z$residuals
  ans$fitted    <- z$fitted
  
  ans$horizon      <- z$horizon
  ans$Wald_Joint   <- z$Wald_Joint
  ans$pv_waldjoint <- 1 - pchisq(z$Wald_Joint, z$df)
  ans$df <- z$df
  
  ans$df.residuals <- z$df.residuals
  rss <- sum(ans$residuals^2)
  mss <- sum(ans$fitted^2)
  n   <- NROW(ans$residuals)
  ans$r.squared     <- mss / (rss + mss)
  ans$adj.r.squared <- 1 - (1 - ans$r.squared) * n / ans$df.residuals
  
  ans$q         <- z$q
  ans$chi_crit  <- qchisq(0.95, ans$q)
  ans$Wald_AR   <- z$Wald_AR
  ans$pv_waldar <- 1 - pchisq(ans$Wald_AR, ans$q)
  
  ans$rse       <- z$rse
  ans$ar_coefs  <- z$ar_coefs
  ans$ar_method <- z$ar_method
  ans$ar_aic    <- z$ar_aic
  
  if (is.null(z$na.action)) ans$na.action <- z$na.action
  class(ans) <- "summary.ivx_ar"
  
  ans
  
}# end-of-function

####################################################
# Print Summary
####################################################

print.summary.ivx_ar <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"),
                                 ...)
{# begin-of-function
 
 cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")
  cat("\n\n",
      if(x$ar_method == "auto") paste0("Auto (", x$ar_aic, ")") else "Fixed",
      " with AR terms q = ", x$q, "\n\n", sep ="")
  
  if (length(x$aliased) == 0L)
  {
    cat("No coefficients\n")
  }
  else
  {    
    coefs_ols <- x$coefficients_ols
    coefs <- x$coefficients
    aliased <- x$aliased
    
    if (!is.null(aliased) && any(aliased)) 
    {
      cn <- names(aliased)
      civx <- x$coefficients_ols
      coefs <- matrix(NA, NROW(civx), 5, dimnames = list(cn , colnames(civx)))
      coefs[!aliased, ] <- civx
    }
    
    cat("Coefficients:\n")
    
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 signif.legend = TRUE, has.Pvalue = TRUE, P.values = TRUE,
                 na.print = "NA", ...)
    
    cat("\nJoint Wald statistic: ", formatC(x$Wald_Joint, digits = digits),
        "on", x$df, "DF, p-value",
        format.pval(x$pv_waldjoint, digits = digits))
    
    cat("\nMultiple R-squared: ", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits))
    
    cat("\nWald AR statistic:", formatC(x$Wald_AR, digits = digits),
        "on", x$q, "DF, p-value",
        format.pval(x$pv_waldar, digits = digits))
    
    cat("\n")
  }
  
  cat("\n")
  invisible(x)
}# end-of-function


####################################################
# importFrom stats embed
####################################################

tilt <- function(x, grid_vec, ar_length) 
{# begin-of-function

  x <- as.matrix(x)
  out <- matrix(NA, nrow = NROW(x) - ar_length, ncol = NCOL(x))
  
  for (i in 1:NCOL(x))
  {
    ds <- embed(x[,i], ar_length + 1)
    out[,i] <- ds[,1] - ds[,-1, drop = FALSE] %*% grid_vec
  }
  colnames(out) <- colnames(x)
  out
  
}# end-of-function

#####################################################################################################
#####################################################################################################
