# =============================================================================
# 03_analysis.R  --  Data Analytics semestral project, topic: Specialty Cafes
#
# Main analysis: research questions A1-A10 of the Specialty Cafes topic.
# Reads analysis/dataset.csv (written by 01_cleaning.R) and exports one table
# and one figure per analysis to report/tables and report/figures.
#
# Run from the REPOSITORY ROOT:  source("analysis/03_analysis.R")
# =============================================================================


# ---- Settings fixed in advance ----------------------------------------------

ALPHA <- 0.05          # significance level for every test
SEED  <- 2026          # seed for every bootstrap
B     <- 5000          # bootstrap resamples

# A1 benchmark: 4 US dollars converted to CZK.
# Exchange rate: Czech National Bank fixing of 5 November 2026,
# 1 USD = 21.646 CZK  ->  4 USD = 86.58 CZK.
USD_CZK   <- 21.646
BENCH_CZK <- 4 * USD_CZK


# ---- Hypotheses (written BEFORE looking at the data) -------------------------
#
# Direction rule: one-tailed only where the question names a direction
# ("charge more", "how much more"); two-tailed otherwise.
#
# A1  One sample. Is the average price of a filter coffee different from 4 USD?
#     H0: mu_filter  = 86.58 CZK
#     H1: mu_filter != 86.58 CZK                                  (two-tailed)
#
# A2  Two groups. Do cafes that offer filter coffee charge more for an espresso
#     than cafes that don't?  Difference = offers filter MINUS does not offer.
#     H0: mu_espresso(filter) - mu_espresso(no filter)  = 0
#     H1: mu_espresso(filter) - mu_espresso(no filter)  > 0       (one-tailed)
#
# A3  Paired. How much more does a filter coffee cost than an espresso in the
#     same cafe?  Difference = filter price MINUS espresso price, per cafe.
#     H0: mu_d  = 0
#     H1: mu_d  > 0                                               (one-tailed)
#
# A4  More groups. Does the espresso price differ between cafes with none, one
#     or both specialty features (roastery, filter coffee)?
#     H0: mu_none = mu_one = mu_both
#     H1: at least one group mean differs                         (omnibus)
#     Post hoc: all pairwise comparisons, adjusted for multiple testing.
#
# A5  Two categorical. Is a roastery more likely to offer filter coffee than
#     other cafes, and by how many percentage points?
#     H0: offering filter coffee is independent of being a roastery
#         (P(filter | roastery) - P(filter | other) = 0)
#     H1: the two are associated                                  (two-sided test)
#     The direction and size are reported as the difference in percentage
#     points, roastery MINUS other cafes.
#
# A6  Correlation. Is the filter coffee price correlated with the espresso
#     price, and in which direction and how strongly?
#     H0: rho  = 0
#     H1: rho != 0                                                (two-tailed)
#
# A7  Simple regression. By how many CZK does the filter coffee price change
#     per 10 CZK increase in the espresso price?
#     Model: filter_price = b0 + b1 * (espresso_price / 10) + e
#     H0: b1  = 0
#     H1: b1 != 0                                                 (two-tailed)
#
# A8  Multiple regression. Does the espresso price differ between cafes that
#     offer filter coffee and those that don't once you account for the
#     distance from the Orloj, and if so, how?
#     Model 1: espresso_price = b0 + b1 * offers_filter + e
#     Model 2: espresso_price = b0 + b1 * offers_filter + b2 * distance_km + e
#     H0: b1  = 0 in Model 2
#     H1: b1 != 0 in Model 2                                      (two-tailed)
#     Also: nested F-test of Model 2 vs Model 1, and the change in b1.
#
# A9  Interaction. Does the relationship between the espresso price and the
#     filter coffee price differ between independent cafes and chains?
#     Model: filter_price = b0 + b1*espresso_price + b2*chain
#                           + b3*espresso_price:chain + e
#     H0: b3  = 0
#     H1: b3 != 0                                                 (two-tailed)
#
# A10 Logistic regression. Does the probability that a cafe offers filter
#     coffee change with each additional km east of the Orloj, and by how
#     many percentage points?
#     Model: logit P(offers_filter = 1) = b0 + b1 * ew_km
#     H0: average marginal effect of ew_km  = 0
#     H1: average marginal effect of ew_km != 0                   (two-tailed)
#
# =============================================================================


# Method choice (from 02_descriptive.R, D_table.tex): every variable and group
# that is compared is either close to normal or large enough for the central
# limit theorem, so the PRIMARY method of every question is parametric. The
# rank-based or bootstrap method is reported as the SECOND method (Appendix D).
# Regressions use HC3 robust standard errors throughout, because the residual
# spread need not be equal across cafes.
#
# Outputs (all written by this script):
#   report/tables/A1_table.tex ... A10_table.tex   one table per analysis
#   report/figures/A1_figure.pdf ... A10_figure.pdf one figure per analysis
#   report/tables/summary.tex        main results of all ten analyses
#   report/tables/D_robustness.tex   second method for every question
#   report/tables/D_posthoc.tex      full A4 post hoc tables
#   report/tables/D_holm.tex         Holm adjustment across A1-A10
#   report/tables/C_influence.tex    influential observations, VIF
#   report/figures/C_diagnostics.pdf residual plots of the linear models
#   report/tables/results_macros.tex every number quoted in the report text
#
# Packages: sandwich, lmtest, car, marginaleffects (install.packages(...)).
# =============================================================================


# ---- 0. Packages, data, helpers ---------------------------------------------

suppressPackageStartupMessages({
  library(sandwich)         # HC3 robust covariance matrices
  library(lmtest)           # coeftest(), coefci(), waldtest() with robust SEs
  library(car)              # vif()
  library(marginaleffects)  # avg_slopes(): average marginal effect (A10)
})

fig_dir <- file.path("report", "figures")
tab_dir <- file.path("report", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path("analysis", "dataset.csv"), stringsAsFactors = FALSE,
              fileEncoding = "UTF-8")
d$specialty_group <- factor(d$specialty_group, levels = c("None", "One", "Both"))
d$chain_group     <- factor(d$chain_group, levels = c("Independent", "Chain"))
d$filter_group    <- factor(d$offers_filter, levels = 0:1,
                            labels = c("No filter", "Offers filter"))
d$roastery_group  <- factor(d$is_roastery, levels = 0:1,
                            labels = c("Other cafe", "Roastery"))
d$esp10 <- d$espresso_price / 10          # espresso price in units of 10 CZK

# Colours for groups (colour-blind safe pair + third hue), grey for single series
COL  <- c("#2a78d6", "#eb6834", "#1baf7a")
GREY <- "grey25"
alpha_col <- function(col, a = 0.45) adjustcolor(col, alpha.f = a)

# --- formatting: numbers go into LaTeX, so a minus sign is typeset properly
fmt <- function(x, digits = 1) {
  s <- formatC(x, format = "f", digits = digits)
  s[round(x, digits) == 0] <- formatC(0, format = "f", digits = digits)  # no "-0.0"
  sub("^-", "\\\\ensuremath{-}", s)
}
fmt_p <- function(p) {
  ifelse(p < 0.001, "< .001", sub("^0\\.", ".", formatC(p, format = "f", digits = 3)))
}
p_tex <- function(p) if (p < 0.001) "\\plt{.001}" else paste0("\\pval{", fmt_p(p), "}")
ci_txt <- function(est, lo, hi, digits = 1) {
  paste0(fmt(est, digits), " [", fmt(lo, digits), ", ", fmt(hi, digits), "]")
}

# --- every number quoted in the report text is stored as a LaTeX macro
macros <- character(0)
put <- function(name, value) {
  stopifnot(grepl("^[A-Za-z]+$", name))
  macros[name] <<- as.character(value)
}

# --- booktabs table writer (no escaping: cells already contain LaTeX)
write_tab <- function(df, file, align, rule_after = integer(0), extra = NULL) {
  body <- apply(as.matrix(df), 1, paste, collapse = " & ")
  body <- gsub("= < ", "< ", paste(body, "\\\\"), fixed = TRUE)  # "p = < .001"
  for (r in rev(rule_after)) body <- append(body, "\\midrule", after = r)
  if (!is.null(extra)) {
    body <- c(body, "\\midrule",
              paste0("\\multicolumn{", ncol(df), "}{l}{", extra, "} \\\\"))
  }
  writeLines(c(paste0("\\begin{tabular}{", align, "}"), "\\toprule",
               paste(paste(names(df), collapse = " & "), "\\\\"), "\\midrule",
               body, "\\bottomrule", "\\end{tabular}"),
             file.path(tab_dir, file))
  cat("Written:", file.path(tab_dir, file), "\n")
}

# --- CI of a standardised effect from the noncentral t distribution
ncp_ci <- function(tval, df, level = 0.95) {
  a <- (1 - level) / 2
  # suppressWarnings: pt() notes reduced precision for far-off trial values of ncp
  f <- function(ncp, prob) suppressWarnings(pt(tval, df, ncp)) - prob
  c(uniroot(f, tval + c(-10, 10), prob = 1 - a, extendInt = "yes")$root,
    uniroot(f, tval + c(-10, 10), prob = a,     extendInt = "yes")$root)
}

# --- CI of eta squared from the noncentral F distribution
eta2_ci <- function(Fval, df1, df2, N, level = 0.95) {
  a <- (1 - level) / 2
  f <- function(ncp, prob) pf(Fval, df1, df2, ncp) - prob
  lo <- if (f(0, 1 - a) <= 0) 0 else uniroot(f, c(0, 1000), prob = 1 - a)$root
  hi <- if (f(0, a)     <= 0) 0 else uniroot(f, c(0, 1000), prob = a)$root
  c(lo, hi) / (c(lo, hi) + N)
}

# --- percentile bootstrap: stat_fun(indices) is evaluated B times
boot_pct <- function(n, stat_fun, seed) {
  set.seed(seed)
  bs <- replicate(B, stat_fun(sample.int(n, replace = TRUE)))
  quantile(bs, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
}

# --- mean with t-based 95% CI
mean_ci <- function(x) {
  x <- x[!is.na(x)]; n <- length(x); se <- sd(x) / sqrt(n)
  c(mean = mean(x), lo = mean(x) - qt(.975, n - 1) * se,
    hi = mean(x) + qt(.975, n - 1) * se, n = n)
}

# --- robust (HC3) coefficient table of a linear model
robust <- function(m) {
  V  <- vcovHC(m, type = "HC3")
  ct <- coeftest(m, vcov. = V)
  ci <- coefci(m, vcov. = V)
  data.frame(term = rownames(ct), b = ct[, 1], se = ct[, 2], t = ct[, 3],
             p = ct[, 4], lo = ci[, 1], hi = ci[, 2], row.names = NULL)
}

summary_rows <- list()   # one row per analysis for summary.tex
robust_rows  <- list()   # second method per analysis for D_robustness.tex
primary_p    <- c()      # primary p-values, for the Holm check
add_rob <- function(analysis, method, result, agrees) {
  robust_rows[[length(robust_rows) + 1]] <<- data.frame(
    Analysis = analysis, `Second method` = method, `Result` = result,
    `Agrees?` = ifelse(agrees, "Yes", "No"), check.names = FALSE)
}
sig <- function(p) p < ALPHA

# Strip plot used by A1 and A3: every cafe as a point, mean with 95% CI above,
# and a vertical reference line.
strip_plot <- function(x, ref, ref_label, xlab, main) {
  set.seed(SEED)
  mc <- mean_ci(x)
  plot(x, jitter(rep(0, length(x)), amount = 0.25), pch = 16,
       col = alpha_col(GREY, 0.5), ylim = c(-0.45, 1.05), yaxt = "n", ylab = "",
       xlab = xlab, main = main, xlim = range(c(x, ref)), bty = "l")
  abline(v = ref, lty = 2, lwd = 1.5, col = COL[2])
  text(ref, 0.98, ref_label, pos = 4, col = "grey20", cex = 0.85)
  segments(mc["lo"], 0.6, mc["hi"], 0.6, lwd = 2.5, col = COL[1])
  points(mc["mean"], 0.6, pch = 21, bg = COL[1], col = "white", cex = 1.8)
  text(mc["mean"], 0.75, sprintf("Mean %.1f, 95%% CI [%.1f, %.1f]",
                                 mc["mean"], mc["lo"], mc["hi"]), cex = 0.85)
}

# Group plot used by A2 and A4: cafes as jittered points, mean and 95% CI.
group_plot <- function(y, g, ylab, main) {
  set.seed(SEED)
  lv <- levels(droplevels(g)); k <- length(lv)
  plot(NA, xlim = c(0.5, k + 0.5), ylim = range(y), xaxt = "n", xlab = "",
       ylab = ylab, main = main, bty = "l")
  for (i in seq_len(k)) {
    yi <- y[g == lv[i]]
    points(jitter(rep(i, length(yi)), amount = 0.18) - 0.1, yi, pch = 16,
           col = alpha_col(COL[i], 0.4))
    mc <- mean_ci(yi)
    segments(i + 0.22, mc["lo"], i + 0.22, mc["hi"], lwd = 2.5, col = COL[i])
    points(i + 0.22, mc["mean"], pch = 21, bg = COL[i], col = "white", cex = 1.8)
  }
  axis(1, at = seq_len(k), tick = FALSE,
       labels = paste0(lv, "\n(n = ", table(droplevels(g))[lv], ")"), padj = 0.5)
}

f <- d[!is.na(d$filter_price), ]          # cafes that offer filter coffee (54)
cat("Cafes with a filter price:", nrow(f), "\n")
put("BenchCZK", fmt(BENCH_CZK, 2)); put("USDCZK", fmt(USD_CZK, 3))
put("Nall", nrow(d)); put("Nfilter", nrow(f)); put("Boot", format(B, big.mark = ","))
put("Seed", SEED)


# ---- A1: one-sample t-test, filter coffee vs 4 USD ---------------------------

x  <- f$filter_price; n <- length(x)
t1 <- t.test(x, mu = BENCH_CZK)
mc <- mean_ci(x)
d1 <- (mean(x) - BENCH_CZK) / sd(x)
d1_ci <- ncp_ci(t1$statistic, n - 1) / sqrt(n)
b1_ci <- boot_pct(n, function(i) mean(x[i]), SEED)
w1 <- wilcox.test(x, mu = BENCH_CZK, conf.int = TRUE, exact = FALSE)

put("AoneN", n)
put("AoneEst", sprintf("\\meanci{%s}{%s}{%s}{CZK}", fmt(mc["mean"]), fmt(mc["lo"]), fmt(mc["hi"])))
put("AoneDiff", sprintf("\\mdiffci{%s}{%s}{%s}{CZK}", fmt(mean(x) - BENCH_CZK),
                        fmt(t1$conf.int[1] - BENCH_CZK), fmt(t1$conf.int[2] - BENCH_CZK)))
put("AoneD", sprintf("\\dci{%s}{%s}{%s}", fmt(d1, 2), fmt(d1_ci[1], 2), fmt(d1_ci[2], 2)))
put("AoneT", sprintf("\\tstat{%d}{%s}, %s", n - 1, fmt(t1$statistic, 2), p_tex(t1$p.value)))
put("AoneMedian", fmt(median(x)))
put("AoneWil", sprintf("$V = %s$, %s", fmt(w1$statistic, 0), p_tex(w1$p.value)))
put("AoneBoot", ci_txt(mean(x), b1_ci[1], b1_ci[2]))
primary_p["A1"] <- t1$p.value

write_tab(data.frame(
  `$n$` = n, `Mean [95\\% CI] (CZK)` = ci_txt(mc["mean"], mc["lo"], mc["hi"]),
  `Median [IQR]` = sprintf("%s [%s, %s]", fmt(median(x)), fmt(quantile(x, .25)), fmt(quantile(x, .75))),
  `Mean $-$ benchmark [95\\% CI]` = ci_txt(mean(x) - BENCH_CZK, t1$conf.int[1] - BENCH_CZK, t1$conf.int[2] - BENCH_CZK),
  `$d$ [95\\% CI]` = ci_txt(d1, d1_ci[1], d1_ci[2], 2),
  `$t(df)$` = sprintf("$t(%d) = %s$", n - 1, fmt(t1$statistic, 2)),
  `$p$` = fmt_p(t1$p.value), check.names = FALSE), "A1_table.tex", "rrrrrrr")

pdf(file.path(fig_dir, "A1_figure.pdf"), width = 7, height = 3.6)
par(mar = c(4.5, 1, 2.5, 1))
strip_plot(x, BENCH_CZK, sprintf("4 USD = %.1f CZK", BENCH_CZK),
           "Price of a filter coffee (CZK)", "")
dev.off()

summary_rows$A1 <- c("A1", "Filter coffee price vs 4 USD", "One-sample $t$", n,
                     paste0("$M = ", fmt(mc["mean"]), "$ [", fmt(mc["lo"]), ", ", fmt(mc["hi"]), "] CZK"),
                     paste0("$d = ", fmt(d1, 2), "$ [", fmt(d1_ci[1], 2), ", ", fmt(d1_ci[2], 2), "]"),
                     fmt_p(t1$p.value))
add_rob("A1", "Wilcoxon signed-rank; bootstrap mean",
        sprintf("$p = %s$; bootstrap CI %s CZK", fmt_p(w1$p.value),
                ci_txt(mean(x), b1_ci[1], b1_ci[2])),
        sig(w1$p.value) == sig(t1$p.value))


# ---- A2: Welch t-test, espresso price by filter coffee -----------------------

a2 <- d[!is.na(d$offers_filter), ]
y1 <- a2$espresso_price[a2$offers_filter == 1]; n1 <- length(y1)
y0 <- a2$espresso_price[a2$offers_filter == 0]; n0 <- length(y0)
t2  <- t.test(y1, y0, alternative = "greater")       # Welch, one-tailed (H1)
t2b <- t.test(y1, y0)                                # two-sided 95% CI
diff2 <- mean(y1) - mean(y0)
sp <- sqrt(((n1 - 1) * var(y1) + (n0 - 1) * var(y0)) / (n1 + n0 - 2))
d2 <- diff2 / sp
d2_ci <- ncp_ci(d2 / sqrt(1 / n1 + 1 / n0), n1 + n0 - 2) * sqrt(1 / n1 + 1 / n0)
set.seed(SEED)
bd2 <- replicate(B, mean(sample(y1, replace = TRUE)) - mean(sample(y0, replace = TRUE)))
b2_ci <- quantile(bd2, c(.025, .975), names = FALSE)
w2  <- wilcox.test(y1, y0, alternative = "greater", exact = FALSE)
w2b <- wilcox.test(y1, y0, conf.int = TRUE, exact = FALSE)

put("AtwoNone", n1); put("AtwoNzero", n0)
put("AtwoMone", fmt(mean(y1))); put("AtwoMzero", fmt(mean(y0)))
put("AtwoEst", sprintf("\\diffci{%s}{%s}{%s}{CZK}{%s}{%s}{%s}", fmt(diff2),
                       fmt(t2b$conf.int[1]), fmt(t2b$conf.int[2]),
                       fmt(d2, 2), fmt(d2_ci[1], 2), fmt(d2_ci[2], 2)))
put("AtwoT", sprintf("\\tstat{%s}{%s}, %s (one-tailed)", fmt(t2$parameter),
                     fmt(t2$statistic, 2), p_tex(t2$p.value)))
put("AtwoLower", fmt(t2$conf.int[1]))
put("AtwoWil", sprintf("$W = %s$, %s (one-tailed)", fmt(w2$statistic, 0), p_tex(w2$p.value)))
put("AtwoHL", ci_txt(w2b$estimate, w2b$conf.int[1], w2b$conf.int[2]))
put("AtwoBoot", ci_txt(diff2, b2_ci[1], b2_ci[2]))
primary_p["A2"] <- t2$p.value

lev2 <- leveneTest(espresso_price ~ filter_group, data = a2)
put("AtwoLevene", sprintf("\\Fstat{%d}{%d}{%s}, %s", lev2[1, "Df"], lev2[2, "Df"],
                          fmt(lev2[1, "F value"], 2), p_tex(lev2[1, "Pr(>F)"])))
write_tab(data.frame(
  `$n$ (filter / no filter)` = sprintf("%d / %d", n1, n0),
  `Means (CZK)` = sprintf("%s / %s", fmt(mean(y1)), fmt(mean(y0))),
  `Difference [95\\% CI] (CZK)` = ci_txt(diff2, t2b$conf.int[1], t2b$conf.int[2]),
  `$d$ [95\\% CI]` = ci_txt(d2, d2_ci[1], d2_ci[2], 2),
  `Welch $t(df)$` = sprintf("$t(%s) = %s$", fmt(t2$parameter), fmt(t2$statistic, 2)),
  `$p$ (one-tailed)` = fmt_p(t2$p.value), check.names = FALSE), "A2_table.tex", "rrrrrr")

pdf(file.path(fig_dir, "A2_figure.pdf"), width = 5.5, height = 4.5)
par(mar = c(4.5, 4.5, 1, 1))
group_plot(a2$espresso_price, a2$filter_group, "Price of an espresso (CZK)", "")
dev.off()

summary_rows$A2 <- c("A2", "Espresso: filter minus no filter", "Welch $t$ (one-tailed)",
                     n1 + n0,
                     paste0("$M_\\Delta = ", fmt(diff2), "$ [", fmt(t2b$conf.int[1]), ", ",
                            fmt(t2b$conf.int[2]), "] CZK"),
                     paste0("$d = ", fmt(d2, 2), "$ [", fmt(d2_ci[1], 2), ", ", fmt(d2_ci[2], 2), "]"),
                     fmt_p(t2$p.value))
add_rob("A2", "Mann--Whitney (one-tailed); bootstrap difference",
        sprintf("$p = %s$; bootstrap CI %s CZK", fmt_p(w2$p.value),
                ci_txt(diff2, b2_ci[1], b2_ci[2])),
        sig(w2$p.value) == sig(t2$p.value))


# ---- A3: paired t-test, filter minus espresso in the same cafe ---------------

dd <- f$filter_price - f$espresso_price; n <- length(dd)
t3  <- t.test(f$filter_price, f$espresso_price, paired = TRUE, alternative = "greater")
t3b <- t.test(f$filter_price, f$espresso_price, paired = TRUE)
dz  <- mean(dd) / sd(dd)
dz_ci <- ncp_ci(t3$statistic, n - 1) / sqrt(n)
b3_ci <- boot_pct(n, function(i) mean(dd[i]), SEED)
w3  <- wilcox.test(dd, alternative = "greater", exact = FALSE)
w3b <- wilcox.test(dd, conf.int = TRUE, exact = FALSE)
# Sensitivity: the one extreme difference (75 CZK) removed
dd_trim <- dd[dd != max(dd)]
t3s <- t.test(dd_trim)

put("AthreeN", n)
put("AthreeEst", sprintf("\\meanci{%s}{%s}{%s}{CZK}", fmt(mean(dd)), fmt(t3b$conf.int[1]),
                         fmt(t3b$conf.int[2])))
put("AthreeDz", sprintf("$d_z = %s$, 95\\%% CI [%s, %s]", fmt(dz, 2), fmt(dz_ci[1], 2), fmt(dz_ci[2], 2)))
put("AthreeT", sprintf("\\tstat{%d}{%s}, %s (one-tailed)", n - 1, fmt(t3$statistic, 2), p_tex(t3$p.value)))
put("AthreeMedian", fmt(median(dd)))
put("AthreeZero", sum(dd == 0)); put("AthreeMax", fmt(max(dd), 0))
put("AthreeTrim", ci_txt(mean(dd_trim), t3s$conf.int[1], t3s$conf.int[2]))
put("AthreeWil", sprintf("$V = %s$, %s (one-tailed)", fmt(w3$statistic, 0), p_tex(w3$p.value)))
put("AthreeBoot", ci_txt(mean(dd), b3_ci[1], b3_ci[2]))
primary_p["A3"] <- t3$p.value

write_tab(data.frame(
  `$n$` = n, `Mean difference [95\\% CI] (CZK)` = ci_txt(mean(dd), t3b$conf.int[1], t3b$conf.int[2]),
  `Median [IQR]` = sprintf("%s [%s, %s]", fmt(median(dd)), fmt(quantile(dd, .25)), fmt(quantile(dd, .75))),
  `$d_z$ [95\\% CI]` = ci_txt(dz, dz_ci[1], dz_ci[2], 2),
  `$t(df)$` = sprintf("$t(%d) = %s$", n - 1, fmt(t3$statistic, 2)),
  `$p$ (one-tailed)` = fmt_p(t3$p.value), check.names = FALSE), "A3_table.tex", "rrrrrr")

pdf(file.path(fig_dir, "A3_figure.pdf"), width = 7, height = 3.6)
par(mar = c(4.5, 1, 2.5, 1))
strip_plot(dd, 0, "No difference", "Filter coffee price minus espresso price, same cafe (CZK)", "")
dev.off()

summary_rows$A3 <- c("A3", "Filter minus espresso, same caf\\'e", "Paired $t$ (one-tailed)", n,
                     paste0("$M_d = ", fmt(mean(dd)), "$ [", fmt(t3b$conf.int[1]), ", ",
                            fmt(t3b$conf.int[2]), "] CZK"),
                     paste0("$d_z = ", fmt(dz, 2), "$ [", fmt(dz_ci[1], 2), ", ", fmt(dz_ci[2], 2), "]"),
                     fmt_p(t3$p.value))
add_rob("A3", "Wilcoxon signed-rank (one-tailed); bootstrap mean",
        sprintf("$p = %s$; bootstrap CI %s CZK", fmt_p(w3$p.value), ci_txt(mean(dd), b3_ci[1], b3_ci[2])),
        sig(w3$p.value) == sig(t3$p.value))


# ---- A4: Welch ANOVA, espresso price by number of specialty features ---------
# Primary: Welch's one-way ANOVA with Games-Howell post hoc tests, chosen in
# advance because the three groups are very unequal in size (125 / 49 / 13).

a4 <- d[!is.na(d$specialty_group), ]
N4 <- nrow(a4); k <- nlevels(a4$specialty_group)
w4 <- oneway.test(espresso_price ~ specialty_group, data = a4, var.equal = FALSE)
c4 <- anova(lm(espresso_price ~ specialty_group, data = a4))
eta2 <- c4[1, "Sum Sq"] / sum(c4[, "Sum Sq"])
eta2_int <- eta2_ci(c4[1, "F value"], c4[1, "Df"], c4[2, "Df"], N4)
lev <- leveneTest(espresso_price ~ specialty_group, data = a4)
kw  <- kruskal.test(espresso_price ~ specialty_group, data = a4)

games_howell <- function(y, g) {
  lv <- levels(g); k <- length(lv)
  out <- list()
  for (j in 1:(k - 1)) for (i in (j + 1):k) {
    yi <- y[g == lv[i]]; yj <- y[g == lv[j]]
    vi <- var(yi) / length(yi); vj <- var(yj) / length(yj)
    se <- sqrt(vi + vj)
    df <- (vi + vj)^2 / (vi^2 / (length(yi) - 1) + vj^2 / (length(yj) - 1))
    dif <- mean(yi) - mean(yj)
    q  <- abs(dif) / se * sqrt(2)
    hw <- qtukey(0.95, k, df) / sqrt(2) * se
    out[[length(out) + 1]] <- data.frame(comparison = paste(lv[i], "$-$", lv[j]),
      diff = dif, lo = dif - hw, hi = dif + hw, df = df,
      p = ptukey(q, k, df, lower.tail = FALSE))
  }
  do.call(rbind, out)
}
gh  <- games_howell(a4$espresso_price, a4$specialty_group)
tk  <- TukeyHSD(aov(espresso_price ~ specialty_group, data = a4))$specialty_group
pw  <- pairwise.wilcox.test(a4$espresso_price, a4$specialty_group,
                            p.adjust.method = "holm", exact = FALSE)$p.value

means4 <- t(sapply(levels(a4$specialty_group),
                   function(g) mean_ci(a4$espresso_price[a4$specialty_group == g])))
put("AfourN", N4)
for (g in rownames(means4)) {
  nm <- c(None = "None", One = "One", Both = "Both")[g]
  put(paste0("AfourM", nm), sprintf("\\meanci{%s}{%s}{%s}{CZK}", fmt(means4[g, "mean"]),
                                    fmt(means4[g, "lo"]), fmt(means4[g, "hi"])))
  put(paste0("AfourN", nm), means4[g, "n"])
}
put("AfourF", sprintf("\\Fstat{%d}{%s}{%s}, %s", w4$parameter[1], fmt(w4$parameter[2]),
                      fmt(w4$statistic, 2), p_tex(w4$p.value)))
put("AfourEta", sprintf("$\\eta^2 = %s$, 95\\%% CI [%s, %s]", fmt(eta2, 3),
                        fmt(eta2_int[1], 3), fmt(eta2_int[2], 3)))
put("AfourClassic", sprintf("\\Fstat{%d}{%d}{%s}, %s", c4[1, "Df"], c4[2, "Df"],
                            fmt(c4[1, "F value"], 2), p_tex(c4[1, "Pr(>F)"])))
put("AfourLevene", sprintf("\\Fstat{%d}{%d}{%s}, %s", lev[1, "Df"], lev[2, "Df"],
                           fmt(lev[1, "F value"], 2), p_tex(lev[1, "Pr(>F)"])))
put("AfourKW", sprintf("$H(%d) = %s$, %s", kw$parameter, fmt(kw$statistic, 2), p_tex(kw$p.value)))
for (r in seq_len(nrow(gh))) {
  nm <- c("OneNone", "BothNone", "BothOne")[r]
  put(paste0("AfourGH", nm), sprintf("$M_\\Delta = %s$~CZK, 95\\%% CI [%s, %s], $p_{\\text{adj}} %s$",
    fmt(gh$diff[r]), fmt(gh$lo[r]), fmt(gh$hi[r]),
    ifelse(gh$p[r] < .001, "< .001", paste("=", fmt_p(gh$p[r])))))
}
primary_p["A4"] <- w4$p.value

write_tab(data.frame(
  Comparison = gh$comparison,
  `Difference [95\\% CI] (CZK)` = ci_txt(gh$diff, gh$lo, gh$hi),
  `$df$` = fmt(gh$df),
  `$p_{\\text{adj}}$` = fmt_p(gh$p), check.names = FALSE), "A4_table.tex", "lrrr",
  extra = sprintf("Welch's ANOVA: $F(%d, %s) = %s$, $p = %s$; $\\eta^2 = %s$, 95\\%% CI [%s, %s]; $N = %d$",
                  w4$parameter[1], fmt(w4$parameter[2]), fmt(w4$statistic, 2), fmt_p(w4$p.value),
                  fmt(eta2, 3), fmt(eta2_int[1], 3), fmt(eta2_int[2], 3), N4))

# Full post hoc tables for Appendix D
tk_names <- rownames(tk)
write_tab(data.frame(
  Comparison = gsub("-", " $-$ ", tk_names),
  `Games--Howell` = paste0(ci_txt(gh$diff, gh$lo, gh$hi), ", $p = ", fmt_p(gh$p), "$"),
  `Tukey HSD` = paste0(ci_txt(tk[, "diff"], tk[, "lwr"], tk[, "upr"]), ", $p = ", fmt_p(tk[, "p adj"]), "$"),
  `Wilcoxon, Holm $p$` = c(fmt_p(pw["One", "None"]), fmt_p(pw["Both", "None"]), fmt_p(pw["Both", "One"])),
  check.names = FALSE), "D_posthoc.tex", "llll")
stopifnot(identical(gsub(" \\$-\\$ ", "-", gh$comparison), tk_names))

pdf(file.path(fig_dir, "A4_figure.pdf"), width = 6.5, height = 4.5)
par(mar = c(4.5, 4.5, 1, 1))
group_plot(a4$espresso_price, a4$specialty_group, "Price of an espresso (CZK)", "")
dev.off()

summary_rows$A4 <- c("A4", "Espresso by specialty features", "Welch ANOVA", N4,
                     paste0("Both $-$ None: ", fmt(gh$diff[2]), " [", fmt(gh$lo[2]), ", ",
                            fmt(gh$hi[2]), "] CZK"),
                     paste0("$\\eta^2 = ", fmt(eta2, 3), "$ [", fmt(eta2_int[1], 3), ", ",
                            fmt(eta2_int[2], 3), "]"),
                     fmt_p(w4$p.value))
add_rob("A4", "Classic ANOVA + Tukey; Kruskal--Wallis",
        sprintf("$F$: $p = %s$; $H$: $p = %s$", fmt_p(c4[1, "Pr(>F)"]), fmt_p(kw$p.value)),
        sig(kw$p.value) == sig(w4$p.value) && sig(c4[1, "Pr(>F)"]) == sig(w4$p.value))


# ---- A5: chi-square test, roastery x filter coffee ---------------------------

a5 <- d[!is.na(d$offers_filter), ]
tab5 <- table(a5$roastery_group, a5$filter_group)
N5  <- sum(tab5)
cs  <- chisq.test(tab5, correct = FALSE)
V5  <- sqrt(unname(cs$statistic) / (N5 * (min(dim(tab5)) - 1)))
V5_ci <- boot_pct(N5, function(i) {
  tb <- table(a5$roastery_group[i], a5$filter_group[i])
  sqrt(suppressWarnings(chisq.test(tb, correct = FALSE)$statistic) / length(i))
}, SEED)
fe  <- fisher.test(tab5)
x_r <- tab5["Roastery", "Offers filter"]; n_r <- sum(tab5["Roastery", ])
x_o <- tab5["Other cafe", "Offers filter"]; n_o <- sum(tab5["Other cafe", ])
pt5 <- prop.test(c(x_r, x_o), c(n_r, n_o), correct = FALSE)
pp  <- 100 * (x_r / n_r - x_o / n_o); pp_ci <- 100 * pt5$conf.int
res <- cs$stdres

put("AfiveN", N5)
put("AfiveRoast", sprintf("%d of %d (%s\\%%)", x_r, n_r, fmt(100 * x_r / n_r)))
put("AfiveOther", sprintf("%d of %d (%s\\%%)", x_o, n_o, fmt(100 * x_o / n_o)))
put("AfivePP", sprintf("%s percentage points, 95\\%% CI [%s, %s]", fmt(pp), fmt(pp_ci[1]), fmt(pp_ci[2])))
put("AfiveChi", sprintf("$\\chi^2(1, N = %d) = %s$, %s", N5, fmt(cs$statistic, 2), p_tex(cs$p.value)))
put("AfiveMinExp", fmt(min(cs$expected)))
put("AfiveV", sprintf("$V = %s$, 95\\%% CI [%s, %s]", fmt(V5, 2), fmt(V5_ci[1], 2), fmt(V5_ci[2], 2)))
put("AfiveRes", fmt(res["Roastery", "Offers filter"], 2))
put("AfiveFisher", sprintf("$\\text{OR} = %s$, 95\\%% CI [%s, %s], %s", fmt(fe$estimate, 2),
                           fmt(fe$conf.int[1], 2), fmt(fe$conf.int[2], 2), p_tex(fe$p.value)))
primary_p["A5"] <- cs$p.value

pct <- function(a, b) sprintf("%d (%s\\%%)", a, fmt(100 * a / b))
write_tab(data.frame(
  ` ` = c("Other caf\\'e", "Roastery", "Total"),
  `No filter` = c(pct(tab5[1, 1], n_o), pct(tab5[2, 1], n_r), pct(sum(tab5[, 1]), N5)),
  `Offers filter` = c(pct(tab5[1, 2], n_o), pct(tab5[2, 2], n_r), pct(sum(tab5[, 2]), N5)),
  Total = c(n_o, n_r, N5),
  `Std.\\ residual (offers filter)` = c(fmt(res[1, 2], 2), fmt(res[2, 2], 2), "--"),
  check.names = FALSE), "A5_table.tex", "lrrrr", rule_after = 2,
  extra = sprintf("$\\chi^2(1, N = %d) = %s$, $p %s$; $V = %s$ [%s, %s]; difference %s pp [%s, %s]; smallest expected count %s",
                  N5, fmt(cs$statistic, 2), ifelse(cs$p.value < .001, "< .001", paste("=", fmt_p(cs$p.value))),
                  fmt(V5, 2), fmt(V5_ci[1], 2), fmt(V5_ci[2], 2), fmt(pp), fmt(pp_ci[1]), fmt(pp_ci[2]),
                  fmt(min(cs$expected))))

wil <- function(x, n) prop.test(x, n, correct = FALSE)$conf.int * 100
ci_r <- wil(x_r, n_r); ci_o <- wil(x_o, n_o)
pdf(file.path(fig_dir, "A5_figure.pdf"), width = 5.5, height = 4.5)
par(mar = c(4.5, 4.5, 1, 1))
bp <- barplot(c(100 * x_o / n_o, 100 * x_r / n_r), names.arg = c(
  sprintf("Other cafe\n(n = %d)", n_o), sprintf("Roastery\n(n = %d)", n_r)),
  ylim = c(0, 100), col = COL[1], border = NA, space = 0.6,
  ylab = "Cafes offering filter coffee (%)")
arrows(bp, c(ci_o[1], ci_r[1]), bp, c(ci_o[2], ci_r[2]), angle = 90, code = 3,
       length = 0.06, lwd = 2, col = "grey20")
text(bp, c(ci_o[2], ci_r[2]) + 4, sprintf("%.0f%%", c(100 * x_o / n_o, 100 * x_r / n_r)))
dev.off()

summary_rows$A5 <- c("A5", "Roastery $\\times$ offers filter", "$\\chi^2$ test", N5,
                     paste0(fmt(pp), " pp [", fmt(pp_ci[1]), ", ", fmt(pp_ci[2]), "]"),
                     paste0("$V = ", fmt(V5, 2), "$ [", fmt(V5_ci[1], 2), ", ", fmt(V5_ci[2], 2), "]"),
                     fmt_p(cs$p.value))
add_rob("A5", "Fisher's exact test",
        sprintf("OR = %s [%s, %s], $p = %s$", fmt(fe$estimate, 2), fmt(fe$conf.int[1], 2),
                fmt(fe$conf.int[2], 2), fmt_p(fe$p.value)),
        sig(fe$p.value) == sig(cs$p.value))


# ---- A6: correlation, filter coffee price and espresso price -----------------

n <- nrow(f)
r6  <- cor.test(f$filter_price, f$espresso_price)
s6  <- cor.test(f$filter_price, f$espresso_price, method = "spearman", exact = FALSE)
k6  <- cor.test(f$filter_price, f$espresso_price, method = "kendall", exact = FALSE)
rb  <- boot_pct(n, function(i) cor(f$filter_price[i], f$espresso_price[i]), SEED)
sb  <- boot_pct(n, function(i) cor(f$filter_price[i], f$espresso_price[i], method = "spearman"), SEED)
kb  <- boot_pct(n, function(i) cor(f$filter_price[i], f$espresso_price[i], method = "kendall"), SEED)

put("AsixN", n)
put("AsixR", sprintf("$r = %s$, 95\\%% CI [%s, %s]", fmt(r6$estimate, 2), fmt(r6$conf.int[1], 2), fmt(r6$conf.int[2], 2)))
put("AsixT", sprintf("\\tstat{%d}{%s}, %s", r6$parameter, fmt(r6$statistic, 2), p_tex(r6$p.value)))
put("AsixRho", sprintf("$\\rho = %s$, 95\\%% CI [%s, %s], %s", fmt(s6$estimate, 2), fmt(sb[1], 2), fmt(sb[2], 2), p_tex(s6$p.value)))
put("AsixTau", sprintf("$\\tau = %s$, 95\\%% CI [%s, %s], %s", fmt(k6$estimate, 2), fmt(kb[1], 2), fmt(kb[2], 2), p_tex(k6$p.value)))
put("AsixBoot", sprintf("[%s, %s]", fmt(rb[1], 2), fmt(rb[2], 2)))
primary_p["A6"] <- r6$p.value

put("AsixRsq", fmt(r6$estimate^2, 2))
write_tab(data.frame(
  `$n$` = n, `Pearson $r$ [95\\% CI]` = ci_txt(r6$estimate, r6$conf.int[1], r6$conf.int[2], 2),
  `$r^2$` = fmt(r6$estimate^2, 2),
  `$t(df)$` = sprintf("$t(%d) = %s$", r6$parameter, fmt(r6$statistic, 2)),
  `$p$` = fmt_p(r6$p.value), check.names = FALSE), "A6_table.tex", "rrrrr")

pdf(file.path(fig_dir, "A6_figure.pdf"), width = 5.5, height = 4.8)
par(mar = c(4.5, 4.5, 1, 1))
set.seed(SEED)
plot(jitter(f$espresso_price, amount = 0.8), jitter(f$filter_price, amount = 0.8),
     pch = 21, bg = alpha_col(COL[1], 0.6), col = "white", cex = 1.4, bty = "l",
     xlab = "Price of an espresso (CZK)", ylab = "Price of a filter coffee (CZK)")
legend("topleft", bty = "n", legend = sprintf("r = %.2f, 95%% CI [%.2f, %.2f], n = %d",
                                              r6$estimate, r6$conf.int[1], r6$conf.int[2], n))
dev.off()

summary_rows$A6 <- c("A6", "Filter and espresso price", "Pearson $r$", n,
                     paste0("$r = ", fmt(r6$estimate, 2), "$ [", fmt(r6$conf.int[1], 2), ", ",
                            fmt(r6$conf.int[2], 2), "]"),
                     paste0("$\\rho = ", fmt(s6$estimate, 2), "$"), fmt_p(r6$p.value))
add_rob("A6", "Spearman $\\rho$; Kendall $\\tau$",
        sprintf("$\\rho = %s$, $p = %s$; $\\tau = %s$, $p = %s$", fmt(s6$estimate, 2),
                fmt_p(s6$p.value), fmt(k6$estimate, 2), fmt_p(k6$p.value)),
        sig(s6$p.value) == sig(r6$p.value) && sig(k6$p.value) == sig(r6$p.value))


# ---- A7: simple regression, filter price on espresso price (per 10 CZK) ------

m7 <- lm(filter_price ~ esp10, data = f)
r7 <- robust(m7)
cl7 <- confint(m7)
s7 <- summary(m7)
b7 <- r7[r7$term == "esp10", ]
set.seed(SEED)
b7_boot <- quantile(replicate(B, {
  i <- sample.int(nrow(f), replace = TRUE); coef(lm(filter_price ~ esp10, data = f[i, ]))[2]
}), c(.025, .975), names = FALSE)
cook7 <- cooks.distance(m7); infl7 <- which(cook7 > 4 / nrow(f))
m7s <- lm(filter_price ~ esp10, data = f[-infl7, ]); r7s <- robust(m7s)

put("AsevenN", nrow(f))
put("AsevenB", sprintf("$b = %s$~CZK per 10~CZK, robust 95\\%% CI [%s, %s]", fmt(b7$b, 2), fmt(b7$lo, 2), fmt(b7$hi, 2)))
put("AsevenT", sprintf("\\tstat{%d}{%s}, %s", m7$df.residual, fmt(b7$t, 2), p_tex(b7$p)))
put("AsevenRsq", fmt(s7$r.squared, 3))
put("AsevenInt", fmt(r7$b[1]))
put("AsevenClassic", sprintf("[%s, %s]", fmt(cl7[2, 1], 2), fmt(cl7[2, 2], 2)))
put("AsevenBoot", sprintf("[%s, %s]", fmt(b7_boot[1], 2), fmt(b7_boot[2], 2)))
put("AsevenInfl", length(infl7))
put("AsevenNoInfl", sprintf("$b = %s$, 95\\%% CI [%s, %s]", fmt(r7s$b[2], 2), fmt(r7s$lo[2], 2), fmt(r7s$hi[2], 2)))
primary_p["A7"] <- b7$p

write_tab(data.frame(
  Term = c("Intercept", "Espresso price (per 10 CZK)"),
  `$b$` = fmt(r7$b, 2), `Robust SE` = fmt(r7$se, 2),
  `Robust 95\\% CI` = sprintf("[%s, %s]", fmt(r7$lo, 2), fmt(r7$hi, 2)),
  `$t(%d)$` = fmt(r7$t, 2), `$p$` = fmt_p(r7$p), check.names = FALSE) |>
  setNames(c("Term", "$b$", "Robust SE", "Robust 95\\% CI", sprintf("$t(%d)$", m7$df.residual), "$p$")),
  "A7_table.tex", "lrrrrr",
  extra = sprintf("$n = %d$, $R^2 = %s$", nobs(m7), fmt(s7$r.squared, 3)))

pdf(file.path(fig_dir, "A7_figure.pdf"), width = 5.5, height = 4.8)
par(mar = c(4.5, 4.5, 1, 1))
gx <- seq(min(f$espresso_price), max(f$espresso_price), length.out = 100)
pr <- predict(m7, newdata = data.frame(esp10 = gx / 10), interval = "confidence")
set.seed(SEED)
plot(jitter(f$espresso_price, amount = 0.8), jitter(f$filter_price, amount = 0.8),
     pch = 21, bg = alpha_col(GREY, 0.45), col = "white", cex = 1.4, bty = "l",
     xlab = "Price of an espresso (CZK)", ylab = "Price of a filter coffee (CZK)")
polygon(c(gx, rev(gx)), c(pr[, "lwr"], rev(pr[, "upr"])), col = alpha_col(COL[1], 0.2), border = NA)
lines(gx, pr[, "fit"], lwd = 2, col = COL[1])
legend("topleft", bty = "n", legend = sprintf("b = %.2f CZK per 10 CZK, R² = %.2f",
                                              b7$b, s7$r.squared))
dev.off()

summary_rows$A7 <- c("A7", "Filter price per 10 CZK espresso", "OLS, HC3 SE", nrow(f),
                     paste0("$b = ", fmt(b7$b, 2), "$ [", fmt(b7$lo, 2), ", ", fmt(b7$hi, 2), "] CZK"),
                     paste0("$R^2 = ", fmt(s7$r.squared, 3), "$"), fmt_p(b7$p))
add_rob("A7", "Pairs bootstrap of the slope; without influential caf\\'es",
        sprintf("CI [%s, %s]; $b = %s$ [%s, %s]", fmt(b7_boot[1], 2), fmt(b7_boot[2], 2),
                fmt(r7s$b[2], 2), fmt(r7s$lo[2], 2), fmt(r7s$hi[2], 2)),
        (b7_boot[1] > 0 || b7_boot[2] < 0) == sig(b7$p) && sig(r7s$p[2]) == sig(b7$p))


# ---- A8: multiple regression, espresso price ~ filter (+ distance) -----------

a8 <- d[!is.na(d$offers_filter), ]
m81 <- lm(espresso_price ~ offers_filter, data = a8)
m82 <- lm(espresso_price ~ offers_filter + distance_km, data = a8)
r81 <- robust(m81); r82 <- robust(m82)
nest  <- anova(m81, m82)
nestR <- waldtest(m81, m82, vcov = function(m) vcovHC(m, type = "HC3"), test = "F")
vif8 <- vif(m82)
m8log <- lm(espresso_price ~ offers_filter + log(distance_km), data = a8); r8log <- robust(m8log)
cook8 <- cooks.distance(m82); infl8 <- which(cook8 > 4 / nrow(a8))
m8s <- lm(espresso_price ~ offers_filter + distance_km, data = a8[-infl8, ]); r8s <- robust(m8s)
b81 <- r81[2, ]; b82 <- r82[2, ]; bd <- r82[3, ]

put("AeightN", nrow(a8))
put("AeightBone", sprintf("$b = %s$~CZK, robust 95\\%% CI [%s, %s], %s", fmt(b81$b, 2), fmt(b81$lo, 2), fmt(b81$hi, 2), p_tex(b81$p)))
put("AeightBtwo", sprintf("$b = %s$~CZK, robust 95\\%% CI [%s, %s], \\tstat{%d}{%s}, %s",
                          fmt(b82$b, 2), fmt(b82$lo, 2), fmt(b82$hi, 2), m82$df.residual, fmt(b82$t, 2), p_tex(b82$p)))
put("AeightDist", sprintf("$b = %s$~CZK per km, robust 95\\%% CI [%s, %s], %s", fmt(bd$b, 2), fmt(bd$lo, 2), fmt(bd$hi, 2), p_tex(bd$p)))
put("AeightChange", fmt(b82$b - b81$b, 2))
put("AeightChangePct", fmt(100 * (b82$b - b81$b) / b81$b, 0))
put("AeightF", sprintf("\\Fstat{%d}{%d}{%s}, %s", nest$Df[2], nest$Res.Df[2], fmt(nest$F[2], 2), p_tex(nest$`Pr(>F)`[2])))
put("AeightFrob", sprintf("\\Fstat{%d}{%d}{%s}, %s", abs(nestR$Df[2]), nestR$Res.Df[2], fmt(nestR$F[2], 2), p_tex(nestR$`Pr(>F)`[2])))
put("AeightVIF", fmt(max(vif8), 2))
put("AeightRsqOne", fmt(summary(m81)$r.squared, 3)); put("AeightRsqTwo", fmt(summary(m82)$r.squared, 3))
put("AeightLog", sprintf("$b = %s$, 95\\%% CI [%s, %s]", fmt(r8log$b[2], 2), fmt(r8log$lo[2], 2), fmt(r8log$hi[2], 2)))
put("AeightInfl", length(infl8))
put("AeightNoInfl", sprintf("$b = %s$, 95\\%% CI [%s, %s]", fmt(r8s$b[2], 2), fmt(r8s$lo[2], 2), fmt(r8s$hi[2], 2)))
primary_p["A8"] <- b82$p

cell <- function(r) c(sprintf("%s (%s)", fmt(r$b, 2), fmt(r$se, 2)),
                      sprintf("[%s, %s]", fmt(r$lo, 2), fmt(r$hi, 2)))
blank <- c("", "")
write_tab(data.frame(
  Term = c("Intercept", "", "Offers filter coffee (1 = yes)", "",
           "Distance from the Orloj (km)", "", "$N$", "$R^2$", "Adjusted $R^2$"),
  `Model 1` = c(cell(r81[1, ]), cell(r81[2, ]), blank, nobs(m81),
                fmt(summary(m81)$r.squared, 3), fmt(summary(m81)$adj.r.squared, 3)),
  `Model 2` = c(cell(r82[1, ]), cell(r82[2, ]), cell(r82[3, ]), nobs(m82),
                fmt(summary(m82)$r.squared, 3), fmt(summary(m82)$adj.r.squared, 3)),
  check.names = FALSE), "A8_table.tex", "lcc", rule_after = 6,
  extra = sprintf("Filter coefficient, Model 2: $t(%d) = %s$, $p = %s$. Nested $F$-test: $F(%d, %d) = %s$, $p = %s$",
                  m82$df.residual, fmt(b82$t, 2), fmt_p(b82$p), nest$Df[2], nest$Res.Df[2],
                  fmt(nest$F[2], 2), fmt_p(nest$`Pr(>F)`[2])))

pdf(file.path(fig_dir, "A8_figure.pdf"), width = 6.5, height = 4.8)
par(mar = c(4.5, 4.5, 1, 1))
grp <- a8$offers_filter + 1
plot(a8$distance_km, a8$espresso_price, pch = 21, bg = alpha_col(COL[grp], 0.55),
     col = "white", cex = 1.3, bty = "l", xlab = "Distance from the Orloj (km)",
     ylab = "Price of an espresso (CZK)")
gx <- seq(0, max(a8$distance_km), length.out = 50)
for (g in 0:1) {
  pg <- predict(m82, data.frame(offers_filter = g, distance_km = gx), interval = "confidence")
  polygon(c(gx, rev(gx)), c(pg[, "lwr"], rev(pg[, "upr"])), col = alpha_col(COL[g + 1], 0.18), border = NA)
  lines(gx, pg[, "fit"], lwd = 2.5, col = COL[g + 1])
}
mtext(sprintf("Filter gap at equal distance: %.1f CZK, 95%% CI [%.1f, %.1f]", b82$b, b82$lo, b82$hi),
      side = 3, line = -1.2, adj = 0.02, cex = 0.85)
legend("bottomleft", bty = "n", pch = 21, pt.bg = COL[1:2], col = "white", pt.cex = 1.4,
       lwd = 2.5, lty = 1, legend = c(sprintf("No filter (n = %d)", sum(a8$offers_filter == 0)),
                                      sprintf("Offers filter (n = %d)", sum(a8$offers_filter == 1))))
dev.off()

summary_rows$A8 <- c("A8", "Filter gap, adjusted for distance", "OLS, HC3 SE", nrow(a8),
                     paste0("$b = ", fmt(b82$b, 2), "$ [", fmt(b82$lo, 2), ", ", fmt(b82$hi, 2), "] CZK"),
                     paste0("$R^2 = ", fmt(summary(m82)$r.squared, 3), "$"), fmt_p(b82$p))
add_rob("A8", "log distance; without influential caf\\'es",
        sprintf("$b = %s$ [%s, %s]; $b = %s$ [%s, %s]", fmt(r8log$b[2], 2), fmt(r8log$lo[2], 2),
                fmt(r8log$hi[2], 2), fmt(r8s$b[2], 2), fmt(r8s$lo[2], 2), fmt(r8s$hi[2], 2)),
        sig(r8log$p[2]) == sig(b82$p) && sig(r8s$p[2]) == sig(b82$p))


# ---- A9: interaction, espresso price x chain on filter price -----------------

m9a <- lm(filter_price ~ esp10 + chain_group, data = f)
m9  <- lm(filter_price ~ esp10 * chain_group, data = f)
r9  <- robust(m9); V9 <- vcovHC(m9, type = "HC3")
bi  <- r9[r9$term == "esp10:chain_groupChain", ]
nest9  <- anova(m9a, m9)
nest9R <- waldtest(m9a, m9, vcov = function(m) vcovHC(m, type = "HC3"), test = "F")
tq9 <- qt(0.975, m9$df.residual)
slope <- function(w) {     # slope of esp10 in a group: w = 0 independent, 1 chain
  L <- c(0, 1, 0, w); est <- sum(L * coef(m9)); se <- sqrt(drop(t(L) %*% V9 %*% L))
  c(est = est, lo = est - tq9 * se, hi = est + tq9 * se, p = 2 * pt(-abs(est / se), m9$df.residual))
}
sI <- slope(0); sC <- slope(1)
nI <- sum(f$chain_group == "Independent"); nC <- sum(f$chain_group == "Chain")
cook9 <- cooks.distance(m9); infl9 <- which(cook9 > 4 / nrow(f))
m9s <- lm(filter_price ~ esp10 * chain_group, data = f[-infl9, ]); r9s <- robust(m9s)

put("AnineN", nrow(f)); put("AnineNind", nI); put("AnineNchain", nC)
put("AnineB", sprintf("$b = %s$~CZK per 10~CZK, robust 95\\%% CI [%s, %s], \\tstat{%d}{%s}, %s",
                      fmt(bi$b, 2), fmt(bi$lo, 2), fmt(bi$hi, 2), m9$df.residual, fmt(bi$t, 2), p_tex(bi$p)))
put("AnineF", sprintf("\\Fstat{%d}{%d}{%s}, %s", nest9$Df[2], nest9$Res.Df[2], fmt(nest9$F[2], 2), p_tex(nest9$`Pr(>F)`[2])))
put("AnineFrob", sprintf("\\Fstat{%d}{%d}{%s}, %s", abs(nest9R$Df[2]), nest9R$Res.Df[2], fmt(nest9R$F[2], 2), p_tex(nest9R$`Pr(>F)`[2])))
put("AnineSind", sprintf("$b = %s$~CZK, 95\\%% CI [%s, %s]", fmt(sI["est"], 2), fmt(sI["lo"], 2), fmt(sI["hi"], 2)))
put("AnineSchain", sprintf("$b = %s$~CZK, 95\\%% CI [%s, %s]", fmt(sC["est"], 2), fmt(sC["lo"], 2), fmt(sC["hi"], 2)))
put("AnineRsq", fmt(summary(m9)$r.squared, 3))
put("AnineInfl", length(infl9))
put("AnineNoInfl", sprintf("$b = %s$, 95\\%% CI [%s, %s]", fmt(r9s$b[4], 2), fmt(r9s$lo[4], 2), fmt(r9s$hi[4], 2)))
put("AnineVIF", fmt(max(vif(m9a)), 2))
primary_p["A9"] <- bi$p

write_tab(data.frame(
  Term = c("Intercept", "Espresso price (per 10 CZK)", "Chain (vs independent)",
           "Espresso $\\times$ chain", "Slope, independent caf\\'es", "Slope, chains"),
  `$b$` = fmt(c(r9$b, sI["est"], sC["est"]), 2),
  `Robust 95\\% CI` = sprintf("[%s, %s]", fmt(c(r9$lo, sI["lo"], sC["lo"]), 2),
                              fmt(c(r9$hi, sI["hi"], sC["hi"]), 2)),
  `$p$` = fmt_p(c(r9$p, sI["p"], sC["p"])), check.names = FALSE),
  "A9_table.tex", "lrrr", rule_after = 4,
  extra = sprintf("Interaction: $t(%d) = %s$, $p = %s$; robust nested $F(%d, %d) = %s$, $p = %s$; $n = %d$, $R^2 = %s$",
                  m9$df.residual, fmt(bi$t, 2), fmt_p(bi$p), abs(nest9R$Df[2]), nest9R$Res.Df[2],
                  fmt(nest9R$F[2], 2), fmt_p(nest9R$`Pr(>F)`[2]), nobs(m9), fmt(summary(m9)$r.squared, 3)))

pdf(file.path(fig_dir, "A9_figure.pdf"), width = 6.5, height = 4.8)
par(mar = c(4.5, 4.5, 1, 1))
gi <- as.integer(f$chain_group)
set.seed(SEED)
plot(jitter(f$espresso_price, amount = 0.8), jitter(f$filter_price, amount = 0.8),
     pch = 21, bg = alpha_col(COL[gi], 0.6), col = "white", cex = 1.4, bty = "l",
     xlab = "Price of an espresso (CZK)", ylab = "Price of a filter coffee (CZK)")
for (g in levels(f$chain_group)) {
  rg <- range(f$espresso_price[f$chain_group == g]); gx <- seq(rg[1], rg[2], length.out = 50)
  pg <- predict(m9, data.frame(esp10 = gx / 10, chain_group = g), interval = "confidence")
  cg <- COL[match(g, levels(f$chain_group))]
  polygon(c(gx, rev(gx)), c(pg[, "lwr"], rev(pg[, "upr"])), col = alpha_col(cg, 0.18), border = NA)
  lines(gx, pg[, "fit"], lwd = 2.5, col = cg)
}
legend("topleft", bty = "n", pch = 21, pt.bg = COL[1:2], col = "white", pt.cex = 1.4,
       lwd = 2.5, legend = c(sprintf("Independent (n = %d)", nI), sprintf("Chain (n = %d)", nC)))
dev.off()

summary_rows$A9 <- c("A9", "Espresso slope: chain vs independent", "OLS interaction, HC3", nrow(f),
                     paste0("$b_{\\times} = ", fmt(bi$b, 2), "$ [", fmt(bi$lo, 2), ", ", fmt(bi$hi, 2), "] CZK"),
                     paste0("$R^2 = ", fmt(summary(m9)$r.squared, 3), "$"), fmt_p(bi$p))
add_rob("A9", "Classic nested $F$; without influential caf\\'es",
        sprintf("$p = %s$; $b_{\\times} = %s$ [%s, %s]", fmt_p(nest9$`Pr(>F)`[2]),
                fmt(r9s$b[4], 2), fmt(r9s$lo[4], 2), fmt(r9s$hi[4], 2)),
        sig(nest9$`Pr(>F)`[2]) == sig(bi$p) && sig(r9s$p[4]) == sig(bi$p))


# ---- A10: logistic regression, offers filter ~ km east of the Orloj ----------

a10 <- d[!is.na(d$offers_filter), ]
m10 <- glm(offers_filter ~ ew_km, family = binomial, data = a10)
ame <- avg_slopes(m10, variables = "ew_km")
or10 <- exp(coef(m10)[2]); or10_ci <- exp(suppressMessages(confint(m10))[2, ])
z10 <- summary(m10)$coefficients
set.seed(SEED)
ame_boot <- quantile(replicate(B, {
  i <- sample.int(nrow(a10), replace = TRUE)
  mb <- glm(offers_filter ~ ew_km, family = binomial, data = a10[i, ])
  p <- fitted(mb); mean(coef(mb)[2] * p * (1 - p))
}), c(.025, .975), names = FALSE)
lpm <- lm(offers_filter ~ ew_km, data = a10); rl <- robust(lpm)
cook10 <- cooks.distance(m10); infl10 <- which(cook10 > 4 / nrow(a10))
m10s <- glm(offers_filter ~ ew_km, family = binomial, data = a10[-infl10, ])
ame_s <- avg_slopes(m10s, variables = "ew_km")

put("AtenN", nrow(a10)); put("AtenEvents", sum(a10$offers_filter))
put("AtenAME", sprintf("%s percentage points per km, 95\\%% CI [%s, %s], $z = %s$, %s",
                       fmt(100 * ame$estimate), fmt(100 * ame$conf.low), fmt(100 * ame$conf.high),
                       fmt(ame$statistic, 2), p_tex(ame$p.value)))
put("AtenOR", sprintf("$\\text{OR} = %s$ per km, 95\\%% CI [%s, %s]", fmt(or10, 2), fmt(or10_ci[1], 2), fmt(or10_ci[2], 2)))
put("AtenBoot", sprintf("[%s, %s]", fmt(100 * ame_boot[1]), fmt(100 * ame_boot[2])))
put("AtenLPM", sprintf("%s percentage points, 95\\%% CI [%s, %s]", fmt(100 * rl$b[2]), fmt(100 * rl$lo[2]), fmt(100 * rl$hi[2])))
put("AtenInfl", length(infl10))
put("AtenNoInfl", sprintf("%s percentage points, 95\\%% CI [%s, %s]", fmt(100 * ame_s$estimate),
                          fmt(100 * ame_s$conf.low), fmt(100 * ame_s$conf.high)))
p_at <- function(km) predict(m10, data.frame(ew_km = km), type = "response")
put("AtenPwest", fmt(100 * p_at(-2), 0)); put("AtenPzero", fmt(100 * p_at(0), 0))
put("AtenPeast", fmt(100 * p_at(4), 0))
primary_p["A10"] <- ame$p.value

write_tab(data.frame(
  Quantity = c("Intercept (log-odds)", "East of the Orloj (log-odds per km)",
               "Odds ratio per km (check)", "Average marginal effect (pp per km)",
               "Linear probability model (pp per km, check)"),
  Estimate = c(fmt(z10[, 1], 2), fmt(or10, 2), fmt(100 * ame$estimate), fmt(100 * rl$b[2])),
  `95\\% CI` = c(sprintf("[%s, %s]", fmt(z10[, 1] - 1.96 * z10[, 2], 2), fmt(z10[, 1] + 1.96 * z10[, 2], 2)),
                 sprintf("[%s, %s]", fmt(or10_ci[1], 2), fmt(or10_ci[2], 2)),
                 sprintf("[%s, %s]", fmt(100 * ame$conf.low), fmt(100 * ame$conf.high)),
                 sprintf("[%s, %s]", fmt(100 * rl$lo[2]), fmt(100 * rl$hi[2]))),
  `$p$` = c(fmt_p(z10[, 4]), fmt_p(z10[2, 4]), fmt_p(ame$p.value), fmt_p(rl$p[2])),
  check.names = FALSE), "A10_table.tex", "lrrr", rule_after = 2,
  extra = sprintf("AME: $z = %s$, $p = %s$; $n = %d$, %d caf\\'es offer filter coffee",
                  fmt(ame$statistic, 2), fmt_p(ame$p.value), nrow(a10), sum(a10$offers_filter)))

pdf(file.path(fig_dir, "A10_figure.pdf"), width = 6.5, height = 4.6)
par(mar = c(4.5, 4.5, 1, 1))
gx <- seq(min(a10$ew_km), max(a10$ew_km), length.out = 200)
pl <- predict(m10, data.frame(ew_km = gx), type = "link", se.fit = TRUE)
set.seed(SEED)
plot(a10$ew_km, jitter(a10$offers_filter * 100, amount = 2.5), pch = 16,
     col = alpha_col(GREY, 0.35), ylim = c(-3, 103), bty = "l",
     xlab = "Distance east of the Orloj (km; negative = west)",
     ylab = "Probability of offering filter coffee (%)")
polygon(c(gx, rev(gx)), 100 * plogis(c(pl$fit - 1.96 * pl$se.fit, rev(pl$fit + 1.96 * pl$se.fit))),
        col = alpha_col(COL[1], 0.2), border = NA)
lines(gx, 100 * plogis(pl$fit), lwd = 2.5, col = COL[1])
abline(v = 0, lty = 3, col = "grey50")
dev.off()

summary_rows$A10 <- c("A10", "Offers filter, per km east", "Logit, AME", nrow(a10),
                      paste0(fmt(100 * ame$estimate), " pp [", fmt(100 * ame$conf.low), ", ",
                             fmt(100 * ame$conf.high), "]"),
                      paste0("OR $= ", fmt(or10, 2), "$ [", fmt(or10_ci[1], 2), ", ", fmt(or10_ci[2], 2), "]"),
                      fmt_p(ame$p.value))
add_rob("A10", "Bootstrap AME; linear probability model",
        sprintf("CI [%s, %s] pp; %s pp [%s, %s]", fmt(100 * ame_boot[1]), fmt(100 * ame_boot[2]),
                fmt(100 * rl$b[2]), fmt(100 * rl$lo[2]), fmt(100 * rl$hi[2])),
        sig(rl$p[2]) == sig(ame$p.value) && ((ame_boot[1] > 0 || ame_boot[2] < 0) == sig(ame$p.value)))


# ---- Summary, robustness, Holm, diagnostics ----------------------------------

summ <- as.data.frame(do.call(rbind, summary_rows))
summ <- cbind(summ[, 1:4], Lost = nrow(d) - as.integer(summ[, 4]), summ[, 5:7])
names(summ) <- c("", "Question", "Method", "$n$", "Lost", "Estimate [95\\% CI]",
                 "Effect size [95\\% CI]", "$p$")
write_tab(summ, "summary.tex", "lllrrlll")

write_tab(do.call(rbind, robust_rows), "D_robustness.tex", "lp{4.8cm}p{7.4cm}c")

holm <- p.adjust(primary_p, method = "holm")
write_tab(data.frame(Analysis = names(primary_p), `Primary $p$` = fmt_p(primary_p),
                     `Holm-adjusted $p$` = fmt_p(holm),
                     `Significant after adjustment?` = ifelse(holm < ALPHA, "Yes", "No"),
                     check.names = FALSE), "D_holm.tex", "lrrc")
put("HolmSig", sum(holm < ALPHA)); put("RawSig", sum(primary_p < ALPHA))

# Diagnostics of the three linear models (Appendix C)
diag_models <- list(A7 = m7, A8 = m82, A9 = m9)
pdf(file.path(fig_dir, "C_diagnostics.pdf"), width = 11, height = 10)
par(mfrow = c(3, 3), mar = c(4.5, 4.5, 3, 1))
for (nm in names(diag_models)) {
  m <- diag_models[[nm]]
  plot(fitted(m), resid(m), pch = 16, col = alpha_col(GREY, 0.5), bty = "l",
       xlab = "Fitted value (CZK)", ylab = "Residual (CZK)",
       main = paste(nm, "- residuals vs fitted"))
  abline(h = 0, lty = 2); lines(lowess(fitted(m), resid(m)), col = COL[2], lwd = 2)
  qqnorm(rstandard(m), pch = 16, col = alpha_col(GREY, 0.6), main = paste(nm, "- normal Q-Q of residuals"))
  qqline(rstandard(m), lwd = 2)
  cd <- cooks.distance(m)
  plot(cd, type = "h", lwd = 2, col = ifelse(cd > 4 / nobs(m), COL[2], "grey50"), bty = "l",
       xlab = "Observation", ylab = "Cook's distance", main = paste(nm, "- Cook's distance"))
  abline(h = 4 / nobs(m), lty = 2)
}
dev.off()

bp7 <- bptest(m7); bp8 <- bptest(m82); bp9 <- bptest(m9)
write_tab(data.frame(
  Model = c("A7: filter $\\sim$ espresso", "A8: espresso $\\sim$ filter + distance",
            "A9: filter $\\sim$ espresso $\\times$ chain", "A10: logit filter $\\sim$ east km"),
  `$n$` = c(nobs(m7), nobs(m82), nobs(m9), nobs(m10)),
  `Cook's $D > 4/n$` = c(length(infl7), length(infl8), length(infl9), length(infl10)),
  `Max VIF` = c("--", fmt(max(vif8), 2), fmt(max(vif(m9a)), 2), "--"),
  `Breusch--Pagan $p$` = c(fmt_p(bp7$p.value), fmt_p(bp8$p.value), fmt_p(bp9$p.value), "--"),
  `Key estimate, all` = c(fmt(b7$b, 2), fmt(b82$b, 2), fmt(bi$b, 2), fmt(100 * ame$estimate)),
  `Key estimate, without influential` = c(fmt(r7s$b[2], 2), fmt(r8s$b[2], 2), fmt(r9s$b[4], 2),
                                          fmt(100 * ame_s$estimate)),
  check.names = FALSE), "C_influence.tex", "lrrrrrr")
put("BPseven", fmt_p(bp7$p.value)); put("BPeight", fmt_p(bp8$p.value)); put("BPnine", fmt_p(bp9$p.value))

# Software versions for Appendix E
put("Rversion", paste(R.version$major, R.version$minor, sep = "."))
for (pk in c("sandwich", "lmtest", "car", "marginaleffects")) {
  put(paste0("Ver", gsub("[^A-Za-z]", "", pk)), as.character(packageVersion(pk)))
}

writeLines(c("% Written by analysis/03_analysis.R -- do not edit by hand.",
             sprintf("\\newcommand{\\%s}{%s}", names(macros), macros)),
           file.path(tab_dir, "results_macros.tex"))
cat("Written:", file.path(tab_dir, "results_macros.tex"), "\n")
cat("\nPrimary p-values and Holm adjustment:\n")
print(round(cbind(p = primary_p, holm = holm), 4))
