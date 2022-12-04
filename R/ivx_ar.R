#####################################################################################################
# Function borrowed from the R package 'ivx'                                                      ###
#                                                                                                 ###
# Reference: Yang, B., Long, W., Peng, L., & Cai, Z. (2020).                                      ###
# Testing the predictability of US housing price index returns based on an IVX-AR model (IVX-AR). ###
# Journal of the American Statistical Association, 115(532), 1598-1619.                           ###
#####################################################################################################

ivx_ar <- function(formula, data, horizon, ar = "auto", ar_ic = c("bic", "aic", "aicc"),
                   ar_max = 5, ar_grid =  function(x) seq(x - 0.3, x + 0.3, by = 0.02),
                   na.action, contrasts = NULL, offset, model = TRUE, x = FALSE, y = FALSE,
                   ...)
{# begin-of-function
  
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  
  if (missing(horizon))
    horizon <- cl$horizon <- 1
  
  if (! ar %in% c("auto", "forecast", c(0:24))) 
  {
    stop("`ar` should be either 'auto' or a non-negative integer.",
         call. = FALSE)
  }
  
  ar_ic <- match.arg(ar_ic)
  
  if( ar_max != trunc(ar_max) || ar_max <= 0 )
  {
    stop("`ar_max` should be a positive integer.", call. = FALSE)
  }
  
  if(!is.function(ar_grid))
  {
    stop("`ar_grid` should be function with sequence generation see `?seq`",
         call. = FALSE)
  }
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "horizon", "na.action", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf["horizon"] <- NULL
  mf <- eval.parent(mf)
  mt <- attr(mf, "terms")
  
  if (attr(mt, "intercept") == 0) 
  {
    warning("ivx estimation does not include an intercept by construction",
            call. = FALSE)
  }
  
  attr(mt, "intercept") <- 0
  
  y <- model.response(mf, "numeric")
  if(is.matrix(y)) {
    stop("multivariate model are not available",call. = FALSE)
  }
  
  ny <- length(y)
  offset <- model.offset(mf)
  x <- model.matrix(mt, mf, contrasts)
  z <- ivx_ar_fit(y, x, horizon = horizon, ar = ar, ar_max = ar_max, ar_ic = ar_ic, ar_grid = ar_grid, offset = offset, ...)
  
  class(z) <- if (ar == 0) "ivx" else c("ivx_ar", "ivx")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$assign <- attr(x, "assign")
  
  if (model)
  {
    z$model <- mf
  }
  
  if (ret.x) 
  {
    z$x <- x
  }
  
  if (ret.y)
  {
    z$y <- y
  }
  z

}# end-of-function
