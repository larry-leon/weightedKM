#' Weighted number of events up to time x
#'
#' @param x Time point.
#' @param error Vector of event/censoring times.
#' @param W Weights for each observation (default: all 1).
#' @return Weighted count of events up to time x.
#' @export
Nkm_Weighted <- function(x, error, W = rep(1, length(error))) {
  sum(W * (error <= x))
}

#' Weighted number at risk at time x
#'
#' @param x Time point.
#' @param error Vector of event/censoring times.
#' @param W Weights for each observation (default: all 1).
#' @return Weighted count of subjects at risk at time x.
#' @export
Rkm_Weighted <- function(x, error, W = rep(1, length(error))) {
  sum(W * (error >= x))
}

#' Weighted Nelson-Aalen and Kaplan-Meier Estimator
#'
#' @param time Vector of event/censoring times.
#' @param Delta Event indicator (1=event, 0=censored).
#' @param W.n Weights for numerator (default: all 1).
#' @param W.d Weights for denominator (default: all 1).
#' @param at.points Time points for estimation.
#' @param se.type Standard error type ("greenwood" or "tsiatis").
#' @param tpoints.add Additional time points to include (default: NULL).
#' @return List with survival estimates, cumulative hazard, standard errors, risk set sizes, and increments.
#' @export
NA_CHRkm_Weighted <- function(time, Delta, W.n = rep(1, length(time)), W.d = rep(1, length(time)),
                              at.points = sort(time), se.type = "greenwood", tpoints.add = NULL) {
  if (!is.null(tpoints.add)) at.points <- sort(unique(c(at.points, tpoints.add)))
  if (!se.type %in% c("greenwood", "tsiatis")) stop("Invalid se type -- greenwood or tsiatis allowed")
  is.sorted <- !is.unsorted(time)
  if (!is.sorted) {
    id <- order(time); time <- time[id]; Delta <- Delta[id]; W.n <- W.n[id]; W.d <- W.d[id]
  }
  risk <- vapply(at.points, Rkm_Weighted, numeric(1), error = time, W = W.d)
  Hmart.chf <- ifelse(risk > 0, 1 / risk, 0)
  counting <- vapply(at.points, Nkm_Weighted, numeric(1), error = time, W = W.n * (Delta == 1))
  counting <- c(0, counting)
  dN <- diff(counting)
  dN.risk <- ifelse(risk > 0, dN / risk, 0.0)
  chf <- cumsum(dN.risk)
  var.chf <- cumsum(ifelse(risk > 0, dN / (risk^2), 0.0))
  S.KM <- cumprod(1 - dN.risk)
  S.KM[S.KM < 0] <- 0.0
  S.NA <- exp(-chf)
  var.NA <- (S.NA^2) * var.chf
  if (se.type == "greenwood") {
    aa <- dN
    bb <- risk * (risk - dN)
    var.KM <- (S.KM^2) * cumsum(ifelse(risk > 0, aa / bb, 0.0))
    se.KM <- sqrt(var.KM)
  } else {
    var.KM <- (S.KM^2) * var.chf
    se.KM <- sqrt(var.KM)
  }
  list(time = time, at.points = at.points, S.NA = S.NA, S.KM = S.KM, chf = chf, se.chf = sqrt(var.chf),
       se.NA = sqrt(var.NA), dN.risk = dN.risk, n.risk = risk, dN = dN, Hmart.chf = Hmart.chf, se.KM = se.KM)
}

#' Calculate weighted median survival time
#'
#' @param times Time points.
#' @param surv Survival probabilities at those times.
#' @param quant Quantile to compute (default: 0.5 for median).
#' @return Time at which survival drops below the specified quantile.
#' @export
weighted_median_time <- function(times, surv, quant = 0.5) {
  idx <- which(surv <= quant)
  if (length(idx) == 0) return(Inf)
  min(times[idx])
}

#' Calculate risk set at specified time points
#'
#' @param times Event/censoring times.
#' @param weights Weights for each subject.
#' @param risk_points Time points for risk set calculation.
#' @return Vector of risk set sizes at each time point.
#' @export
calculate_risk_set <- function(times, weights, risk_points) {
  vapply(risk_points, Rkm_Weighted, numeric(1), error = times, W = weights)
}

#' Format p-value for display
#'
#' @param pval P-value to format.
#' @param eps Threshold for displaying as "<eps".
#' @param digits Number of decimal places.
#' @return Formatted p-value string.
#' @export
format_pval <- function(pval, eps = 0.001, digits = 3) {
  if (is.na(pval)) return(NA)
  if (pval < eps) return(paste0("<", eps))
  format(round(pval, digits), nsmall = digits)
}

#' Validate input data frame and columns
#'
#' @param df Data frame.
#' @param required_cols Vector of required column names.
#' @return Stops with error if columns are missing; otherwise, returns invisibly.
#' @export
validate_input <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  invisible(NULL)
}

#' Prepare group data (extract vectors)
#'
#' @param df Data frame.
#' @param tte.name Time-to-event column name.
#' @param event.name Event indicator column name.
#' @param treat.name Treatment group column name.
#' @param weight.name Weight column name (optional).
#' @param treat_value Value indicating the group (default: 1).
#' @return List with vectors for time, event, and weights.
#' @export
prepare_group_data <- function(df, tte.name, event.name, treat.name, weight.name = NULL, treat_value = 1) {
  idx <- which(df[[treat.name]] == treat_value)
  time <- df[[tte.name]][idx]
  event <- df[[event.name]][idx]
  weights <- if (!is.null(weight.name)) df[[weight.name]][idx] else rep(1, length(idx))
  list(time = time, event = event, weights = weights)
}

#' Compute KM estimates for a group
#'
#' @param time Event/censoring times.
#' @param event Event indicators.
#' @param weights Weights for each subject.
#' @param at.points Time points for estimation.
#' @return List with KM survival, standard error, risk set, and event increments.
#' @export
compute_km_group <- function(time, event, weights, at.points) {
  fit <- NA_CHRkm_Weighted(time = time, Delta = (event == 1), W.n = weights, W.d = weights, at.points = at.points)
  list(S.KM = fit$S.KM, se.KM = fit$se.KM, n.risk = fit$n.risk, dN = fit$dN)
}

#' Get censoring times for a group
#'
#' @param time Event/censoring times.
#' @param event Event indicators.
#' @param at.points Time points for estimation.
#' @return List with censoring times and their indices in at.points.
#' @export
get_censor_times <- function(time, event, at.points) {
  cens <- time[event == 0]
  events <- sort(unique(time[event == 1]))
  cens <- setdiff(cens, events)
  match_idx <- match(cens, at.points)
  list(times = cens, idx = match_idx)
}

#' Plot KM curves for two groups, with censor marks
#'
#' @param at.points1 Time points for group 1.
#' @param S1.KM KM estimates for group 1.
#' @param at.points2 Time points for group 2.
#' @param S2.KM KM estimates for group 2.
#' @param col.1 Color for group 1.
#' @param col.2 Color for group 2.
#' @param ltys Line types.
#' @param lwds Line widths.
#' @param Xlab X-axis label.
#' @param Ylab Y-axis label.
#' @param ylim Y-axis limits.
#' @param xlim X-axis limits.
#' @param show.ticks Show censor marks (default: FALSE).
#' @param censor1 Censoring info for group 1.
#' @param censor2 Censoring info for group 2.
#' @param S1.KM_full Full KM for group 1.
#' @param S2.KM_full Full KM for group 2.
#' @param censor.cex Censor mark size.
#' @param ... Additional graphical parameters.
#' @return Plots the curves; no return value.
#' @export
plot_km_curves <- function(at.points1, S1.KM, at.points2, S2.KM, col.1, col.2, ltys, lwds, Xlab, Ylab, ylim, xlim,
                           show.ticks = FALSE, censor1 = NULL, censor2 = NULL, S1.KM_full = NULL, S2.KM_full = NULL, censor.cex = 1.0, ...) {
  plot(at.points1, S1.KM, type = "s", ylim = ylim, xlim = xlim, lty = ltys[1], col = col.1, lwd = lwds[1], xlab = Xlab, ylab = Ylab, ...)
  lines(at.points2, S2.KM, type = "s", lty = ltys[2], col = col.2, lwd = lwds[2])
  if (show.ticks && !is.null(censor1) && !is.null(S1.KM_full)) {
    points(censor1$times, S1.KM_full[censor1$idx], pch = 3, col = col.1, cex = censor.cex)
  }
  if (show.ticks && !is.null(censor2) && !is.null(S2.KM_full)) {
    points(censor2$times, S2.KM_full[censor2$idx], pch = 3, col = col.2, cex = censor.cex)
  }
}

#' Plot weighted KM curves for two groups (main function)
#'
#' @param df Data frame.
#' @param tte.name Time-to-event column name.
#' @param event.name Event indicator column name.
#' @param treat.name Treatment group column name.
#' @param weight.name Weight column name (optional).
#' @param strata.name Optional stratification.
#' @param show.cox Show Cox model results (default: FALSE).
#' @param cox.cex Cox legend size.
#' @param show.logrank Show logrank test results (default: FALSE).
#' @param logrank.cex Logrank legend size.
#' @param cox.eps Cox p-value threshold.
#' @param lr.eps Logrank p-value threshold.
#' @param show_arm_legend Show arm legend (default: TRUE).
#' @param arms Arm names for legend.
#' @param put.legend.arms Legend position for arms.
#' @param stop.onerror Stop on error (default: TRUE).
#' @param check.KM Check KM estimates (default: TRUE).
#' @param put.legend.cox Cox legend position.
#' @param put.legend.lr Logrank legend position.
#' @param lr.digits Digits for logrank p-value.
#' @param cox.digits Digits for Cox p-value.
#' @param tpoints.add Additional time points.
#' @param by.risk Interval for risk set display.
#' @param Xlab X-axis label.
#' @param Ylab Y-axis label.
#' @param col.1 Color for group 1.
#' @param col.2 Color for group 2.
#' @param show.med Show median survival (default: FALSE).
#' @param choose_ylim Choose y-axis limits (default: FALSE).
#' @param arm.cex Arm legend size.
#' @param quant Quantile for median (default: 0.5).
#' @param qlabel Label for quantile.
#' @param med.cex Median text size.
#' @param ymed.offset Y offset for median text.
#' @param del.med Y offset between medians.
#' @param xmed.offset X offset for median text.
#' @param risk.cex Risk set text size.
#' @param plotit Whether to plot (default: TRUE).
#' @param ltys Line types.
#' @param lwds Line widths.
#' @param censor.mark.all Show all censor marks (default: TRUE).
#' @param censor.cex Censor mark size.
#' @param show.ticks Show censor marks (default: TRUE).
#' @param risk.set Show risk set (default: TRUE).
#' @param ymin Minimum y-axis value.
#' @param ymax Maximum y-axis value.
#' @param ymin.del Y offset for risk set.
#' @param ymin2 Alternative minimum y-axis value.
#' @param risk_offset Y offset for risk set.
#' @param risk_delta Y offset between risk sets.
#' @param y.risk2 Y position for risk set 2.
#' @param show.Y.axis Show y-axis (default: TRUE).
#' @param cex_Yaxis Y-axis text size.
#' @param y.risk1 Y position for risk set 1.
#' @param add.segment Add segment (default: FALSE).
#' @param risk.add Additional risk set points.
#' @param xmin Minimum x-axis value.
#' @param xmax Maximum x-axis value.
#' @param x.truncate Truncate x-axis.
#' @param time.zero Time zero.
#' @param time.zero.label Label for time zero.
#' @param prob.points Probability points.
#' @return List with KM estimates, medians, risk sets, time points, and (if requested) Cox model results.
#' @export
KM_plot_2sample_weighted <- function(df, tte.name, event.name, treat.name, weight.name = NULL,
                                              strata.name = NULL, show.cox = FALSE, cox.cex = 0.725, show.logrank = FALSE,
                                              logrank.cex = 0.725, cox.eps = 0.001, lr.eps = 0.001, show_arm_legend = TRUE,
                                              arms = c("treat", "control"), put.legend.arms = "left", stop.onerror = TRUE,
                                              check.KM = TRUE, put.legend.cox = "topright", put.legend.lr = "top",
                                              lr.digits = 2, cox.digits = 2, tpoints.add = c(0), by.risk = NULL, Xlab = "time",
                                              Ylab = "proportion surviving", col.1 = "black", col.2 = "blue", show.med = FALSE,
                                              choose_ylim = FALSE, arm.cex = 0.7, quant = 0.5, qlabel = "median =", med.cex = 0.725,
                                              ymed.offset = 0.25, del.med = 0.075, xmed.offset = 4, risk.cex = 0.725, plotit = TRUE,
                                              ltys = c(1, 1), lwds = c(2, 2), censor.mark.all = TRUE, censor.cex = 1.0,
                                              show.ticks = TRUE, risk.set = TRUE, ymin = 0, ymax = 1, ymin.del = 0.035,
                                              ymin2 = NULL, risk_offset = 0.075, risk_delta = 0.025, y.risk2 = NULL,
                                              show.Y.axis = TRUE, cex_Yaxis = 1, y.risk1 = NULL, add.segment = FALSE,
                                              risk.add = NULL, xmin = 0, xmax = NULL, x.truncate = NULL, time.zero = 0.0,
                                              time.zero.label = 0.0, prob.points = NULL) {

  validate_input(df, c(tte.name, event.name, treat.name))
  group1 <- prepare_group_data(df, tte.name, event.name, treat.name, weight.name, treat_value = 1)
  group2 <- prepare_group_data(df, tte.name, event.name, treat.name, weight.name, treat_value = 0)
  evs1 <- sort(unique(group1$time[group1$event == 1]))
  evs2 <- sort(unique(group2$time[group2$event == 1]))
  tpoints <- sort(unique(c(time.zero, tpoints.add, evs1, evs2, max(df[[tte.name]]))))
  at.points1 <- sort(unique(c(group1$time, tpoints)))
  at.points2 <- sort(unique(c(group2$time, tpoints)))
  fit1 <- compute_km_group(group1$time, group1$event, group1$weights, at.points1)
  fit2 <- compute_km_group(group2$time, group2$event, group2$weights, at.points2)
  med.1 <- weighted_median_time(at.points1, fit1$S.KM, quant)
  med.2 <- weighted_median_time(at.points2, fit2$S.KM, quant)
  risk.points <- sort(unique(c(seq(time.zero.label, max(c(at.points1, at.points2)), by = ifelse(is.null(by.risk), 1, by.risk)), risk.add)))
  risk.1 <- calculate_risk_set(group1$time, group1$weights, risk.points)
  risk.2 <- calculate_risk_set(group2$time, group2$weights, risk.points)
  censor1 <- get_censor_times(group1$time, group1$event, at.points1)
  censor2 <- get_censor_times(group2$time, group2$event, at.points2)

  # KM curve checks
  if (check.KM) {
    check_km_curve <- function(S.KM, group_name = "Group") {
      if (any(S.KM < 0 | S.KM > 1)) {
        msg <- paste0(group_name, " KM curve has values outside [0,1].")
        if (stop.onerror) stop(msg) else warning(msg)
      }
      if (any(diff(S.KM) > 0)) {
        msg <- paste0(group_name, " KM curve is not non-increasing.")
        if (stop.onerror) stop(msg) else warning(msg)
      }
      if (abs(S.KM[1] - 1) > 1e-6) {
        msg <- paste0(group_name, " KM curve does not start at 1.")
        if (stop.onerror) stop(msg) else warning(msg)
      }
    }
    check_km_curve(fit1$S.KM, "Group 1")
    check_km_curve(fit2$S.KM, "Group 2")
  }


  cox_results <- NULL
  if (show.cox) {
    if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required for Cox model.")
    treat.name.strata <- treat.name
    if (!is.null(strata.name)) {
      treat.name.strata <- paste(treat.name,"+",paste("strata(",eval(strata.name),")"))
    }
   cox_formula <- as.formula(paste0("survival::Surv(", tte.name, ",", event.name, ") ~ ", treat.name.strata))
    if (!is.null(weight.name)) {
      cox_fit <- survival::coxph(cox_formula, data = df, weights = df[[weight.name]])
    } else {
      cox_fit <- survival::coxph(cox_formula, data = df)
    }
    cox_summary <- summary(cox_fit)
    hr <- exp(cox_fit$coef)
    hr_ci <- exp(confint(cox_fit))
    pval <- cox_summary$coefficients[1, "Pr(>|z|)"]
    cox_text <- paste0("HR = ", round(hr, cox.digits),
                       " (", round(hr_ci[1], cox.digits), ", ", round(hr_ci[2], cox.digits), ")",
                       ", p = ", format_pval(pval, eps = cox.eps, digits = cox.digits))
    cox_results <- list(
      cox_fit = cox_fit,
      hr = hr,
      hr_ci = hr_ci,
      pval = pval,
      cox_text = cox_text
    )
  }

  logrank_results <- NULL
  if (show.logrank) {
    if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required for logrank test.")
    surv_obj <- survival::Surv(df[[tte.name]], df[[event.name]])
    group <- df[[treat.name]]
    if (!is.null(strata.name)) {
      strata_var <- df[[strata.name]]
      logrank_formula <- as.formula(paste0("surv_obj ~ group + strata(strata_var)"))
      logrank_fit <- eval(bquote(survival::survdiff(.(logrank_formula))))
    } else {
      logrank_fit <- survival::survdiff(surv_obj ~ group)
    }

    chisq <- logrank_fit$chisq
    pval <- 1 - pchisq(chisq, df = 1)
    logrank_text <- paste0("Logrank p = ", format_pval(pval, eps = lr.eps, digits = lr.digits))
    logrank_results <- list(
      chisq = chisq,
      pval = pval,
      logrank_text = logrank_text
    )
  }


  if (plotit) {
    plot_km_curves(at.points1, fit1$S.KM, at.points2, fit2$S.KM, col.1, col.2, ltys, lwds, Xlab, Ylab,
                   ylim = c(ifelse(is.null(ymin2), ymin - risk_offset, ymin2), ymax),
                   xlim = c(xmin, ifelse(is.null(xmax), max(c(at.points1, at.points2)), xmax)),
                   show.ticks = show.ticks, censor1 = censor1, censor2 = censor2,
                   S1.KM_full = fit1$S.KM, S2.KM_full = fit2$S.KM, censor.cex = censor.cex)
    if (risk.set) {
      text(risk.points, rep(ifelse(is.null(y.risk1), ymin - risk_offset, y.risk1), length(risk.points)),
           round(risk.1), col = col.1, cex = risk.cex)
      text(risk.points, rep(ifelse(is.null(y.risk2), ymin - risk_offset + risk_delta, y.risk2), length(risk.points)),
           round(risk.2), col = col.2, cex = risk.cex)
    }

    # horizontal line just below ymin
    abline(h=ymin-ymin.del,lty=1,col=1)

    if (show.cox) {
      legend(put.legend.cox, legend = cox_text, cex = cox.cex, bty = "n")
    }

    if (show.logrank) {
      legend(put.legend.lr, legend = logrank_text, cex = logrank.cex, bty = "n")
    }


    if (show.med && med.1 != Inf && med.2 != Inf) {
      text(max(c(at.points1, at.points2)) - xmed.offset, ymax - ymed.offset, paste(qlabel, round(med.1, 1)), col = col.1, cex = med.cex)
      text(max(c(at.points1, at.points2)) - xmed.offset, ymax - (ymed.offset + del.med), paste(qlabel, round(med.2, 1)), col = col.2, cex = med.cex)
    }
    if (show_arm_legend) {
      legend(put.legend.arms, legend = arms, col = c(col.1, col.2), lty = ltys, lwd = lwds, cex = arm.cex, bty = "n")
    }
  }
  list(
    S1.KM = fit1$S.KM, S2.KM = fit2$S.KM,
    med.1 = med.1, med.2 = med.2,
    risk.1 = risk.1, risk.2 = risk.2, risk.points = risk.points,
    at.points1 = at.points1, at.points2 = at.points2,
    cox = cox_results,logrank = logrank_results
  )
}
