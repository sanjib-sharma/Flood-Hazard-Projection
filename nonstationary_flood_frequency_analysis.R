library(splines)
library(dplyr)
library(extRemes)
library(lubridate)
library(gdata)
library(geosphere)
library(ggplot2)
# ---------------------------------------------------------------
#load("svs6308/flow projections")
#stations <- c('01553005','01551500') # Select stations to fit
#pens <- 10                       # Penalty parameter
#cprob <- c(0.90,0.95, 0.98, 0.99)     # Threshold quantiles to use
######################################
#' Fit GEV or PP-GPD where parameters vary smoothly in time
#'
#' @param y flow vector
#' @param dates date vector that matches the order of y
#' @param threshold censoring threshold if fitting PP-GPD
#' @param fit.type character string: "GEV" or "PP" (default)
#' @param deg degree of spline
#' @param decl declustering parameter
#' @param lam penalty parameter for spline smoothness (Inf -> linear fit)
#'
#' @return Fitted scale, shape, location, parameters at each date in pred.dates
#'         pred.dates are the first date of each month
fitevd <-
  function(y,
           dates,
           threshold,
           fit.type = c("PP", "GEV"),
           deg = 11,
           decl = 0,
           lam = 100) {
    ord <- order(dates)
    y <- y[ord]
    dates <- dates[ord]
    if (decl > 0) {
      y <- decluster(excs.data$flow, threshold, r = decl)
    }
    y[is.na(y)] <- 0
    bsb <- cbind(1, bs(dates, df = deg))
    if (fit.type == "PP") {
      pred.dates <- dates[(day(dates) == 1) & (month(dates) == 1)]
      bsm <- cbind(1, bs(pred.dates, df = deg))
    }
    else{
      pred.dates <- dates
      bsm <- bsb
    }
    out <-
      tryCatch(
        pspp(
          y,
          bsb,
          lam = lam,
          fit.type = fit.type,
          threshold = threshold
        )
        ,
        finally =  function(x)
          NA
      )
    pgams <- out$par
    deg1 <- ncol(bsm)
    fscale <-  bsm %*% pgams[1:deg1]
    fshape <-  bsm %*% pgams[(deg1 + 1):(2 * deg1)]
    floc <- bsm %*% pgams[(2 * deg1 + 1):(3 * deg1)]
    return(list(
      fscale = fscale,
      fshape = fshape,
      floc = floc,
      pred.dates = pred.dates
    ))
  }

#' Penalized likelihood function
#'
#' @param y flow vector
#' @param bsb matrix of basis coefficients
#' @param lam penalty parameter
#' @param fit.type fit type: one of "PP" or "GEV"
#' @param threshold If using the PP (i.e. GPD), which threshold to use
#' @param max.it Maximum number of iterations to be used by optimization routine
#'
#' @return penalized spline likelihood for passed parameters
pspp <- function(y,
                 bsb,
                 lam = 100,
                 fit.type = c("PP", "GEV"),
                 threshold = threshold,
                 max.it = 1000) {
  lam1 <- lam2 <- lam3 <- lam
  plik <- function(p, y, bsb, lam1, lam2, lam3, tol = 1e-5) {
    nb <- ncol(bsb)
    scl <-  bsb %*% p[1:nb]
    shp <-  bsb %*% p[(nb + 1):(2 * nb)]
    loc <- bsb %*% p[(2 * nb + 1):(3 * nb)]
    if (any(scl < 0)) {
      return(Inf)
    }
    else{
      return(
        levd(
          y,
          threshold = threshold,
          location = loc,
          scale = scl,
          shape = shp,
          type = fit.type,
          log = TRUE,
          negative = TRUE
        ) +
          lam1 * sum(diff(p[2:nb], 2) ^ 2) +
          lam2 * sum(diff(p[(nb + 2):(2 * nb)], 2) ^ 2) +
          lam3 * sum(diff(p[(2 * nb + 2):(3 * nb)], 2) ^ 2)
      )     # Penalty
      
    }
  }
  scl.gams <- shp.gams <- loc.gams <- rep(0, ncol(bsb))
  init <- fevd(y, threshold = threshold, type = fit.type)
  scl.gams[1] <- init$results$par["scale"]
  shp.gams[1] <- init$results$par["shape"]
  loc.gams[1] <- init$result$par["location"]
  init <- c(scl.gams, shp.gams, loc.gams)
  out <-
    optim(
      init,
      fn = plik,
      control = list(maxit = max.it),
      y = y,
      bsb = bsb,
      lam1 = lam1,
      lam2 = lam2,
      lam3 = lam3
    )
  out
}


#' Fit Wrapper that finds both the GEV and PP fits
#'
#' @param y flow Vector
#' @param dates 
#' @param threshold threshold to be used for GPD
#' @param pen penalty parameter to be used on splines
#' @param rp return period to be estimated
#'
#' @return Fitted quantiles "return levels" from GEV and GPD-PP fits that vary smoothly in time. 
#'         The idea of a return level in the traditional sense doesn't apply here. 
#'         
#'         GEV and GPD-PP location, scale, and shape estimates for each year.
fitwrap <- function(y, dates, threshold, pen, rp = 100) {
  out <-
    fitevd(
      y,
      dates,
      threshold,
      fit.type = "PP",
      deg = 5,
      lam = pen
    )
  years <- out$pred.dates
  pploc <- loc <- out$floc
  ppscale <- scale <- out$fscale
  ppshape <- shape <- out$fshape
  pprl <- mapply(function(l, sc, sh) {
    rlevd(
      rp,
      loc = l,
      scale = sc,
      shape = sh,
      npy = 365.25,
      type = "PP"
    )
  },
  l = loc,
  sc = scale,
  sh = shape)
  df <- data.frame(y = y, dates = dates)
  ams <- df  %>%
    group_by(year = year(dates)) %>%
    summarise(y = max(y))
  subyr <- data.frame(year = unique(year(dates[(day(dates) == 1) & (month(dates) == 1)])))
  ams <- as.data.frame(inner_join(ams, subyr, by= "year"))
  gev <-
    fitevd(
      ams$y,
      ams$year,
      threshold,
      fit.type = "GEV",
      deg = 5,
      lam = pen
    )
  gevloc <- loc <- gev$floc
  gevscale <- scale <- gev$fscale
  gevshape <- shape <- gev$fshape
  gevrl <- mapply(function(l, sc, sh) {
    rlevd(
      rp,
      loc = l,
      scale = sc,
      shape = sh,
      npy = 365.25,
      type = "GEV"
    )
  },
  l = loc,
  sc = scale,
  sh = shape)
  return(list(
    years = years,
    ams = ams$y,
    pprl = pprl,
    gevrl = gevrl,
    gevloc = gevloc,
    gevscale = gevscale,
    gevshape = gevshape,
    pploc = pploc,
    ppscale = ppscale,
    ppshape = ppshape
  ))
}

##########################################
# Fit and make plots ------------------------------------------------------
for(k in 1:length(stations)) {
  inid <- stations[k]
  sub <- dailyflow %>% filter(ID == inid)
  y <- sub$flow
  dates <- sub$date
  pdf(paste0("./daily_", inid, ".pdf"))
  par(mfrow = c(2, 2))
  plot(
    1,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  legend(
    "bottomleft",
    legend = c(
      "GPD RL50",
      "GPD RL100",
      "GPD RL200",
      "GEV RL50",
      "GEV RL100",
      "GEV RL200"
 ),
    lty = c(1, 1, 1, 2, 2, 2),
    col = c(
      "chocolate",
      "dodgerblue",
      "seagreen",
      "chocolate",
      "dodgerblue",
      "seagreen"
    )
  )
  for (i in 1:length(pens)) {
    for (j in 1:length(cprob)) {
      thlds <- unique(quantile(y[y > 0], cprob))
      fits <-
        fitwrap(y,
                dates,
                thlds[j],
                pen = pens[i],
                rp = c(50, 100, 200))
      years <- fits$years
      ams <- fits$ams
      pprl <- t(fits$pprl)
      gevrl <- t(fits$gevrl)
      #################

gevloc <- t(fits$gevloc)
 pploc <-t(fits$pploc)
gevscale<-t(fits$gevscale)
 ppscale <-t(fits$ppscale)
gevshape <-t(fits$gevshape)
ppshape <-t(fits$ppshape)
################   
      minlen <- min(length(ams), length(years))
      plot(
        years[1:minlen],
        ams[1:minlen],
        ylim = range(c(ams, pprl, gevrl)),
        xlab = "year",
        ylab = "AMS (in)",
        main = paste("Thrsh. Quant.: ", cprob[j])
      )
      lines(years, pprl[, 1], col = "chocolate")
      lines(years, pprl[, 2], col = "dodgerblue")
      lines(years, pprl[, 3], col = "seagreen")
      lines(years[1:minlen], gevrl[1:minlen, 1], lty = 2, col = "chocolate")
      lines(years[1:minlen], gevrl[1:minlen, 2], lty = 2, col = "dodgerblue")
      lines(years[1:minlen], gevrl[1:minlen, 3], lty = 2, col = "seagreen")
    }
  }
  dev.off()
}

