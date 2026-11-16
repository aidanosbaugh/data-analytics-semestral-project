# =============================================================================
# 02_descriptive.R  --  Data Analytics semestral project, topic: Specialty Cafes
#
# Describes every variable used in A1-A10 BEFORE any test is run, and decides
# for each variable (and each group that is compared) whether a parametric or
# a non-parametric method is justified. Steps follow Appendix B.2:
#   1. Count   2. Summarise   3. Plot   4. Check normality
#   5. If not normal: does the central limit theorem (CLT) help?
#   6. Fill in the decision table (D_table.tex)
#
# Reads:   analysis/dataset.csv
# Writes:  report/tables/D_table.tex        decision table (main output)
#          report/tables/descriptives.tex   full numeric descriptives (App. A)
#          report/tables/D_categorical.tex  counts and percentages     (App. A)
#          report/figures/D_*.pdf           descriptive plots
#          report/figures/distributions.pdf summary figure for Appendix B
#
# Run from the REPOSITORY ROOT:  source("analysis/02_descriptive.R")
# Packages: none (base R).
# =============================================================================


# ---- 0. Settings and data ----------------------------------------------------

SEED <- 2026     # same seed as 03_analysis.R
B    <- 5000     # bootstrap resamples for the CLT check

# Rules of thumb used for the decision (stated in the report, Appendix B):
SKEW_LIMIT <- 0.5   # |skewness| below this -> distribution treated as normal
N_CLT      <- 30    # smallest group size for which the CLT can help
BOOT_LIMIT <- 0.5   # |skewness of bootstrapped means| below this -> CLT works

fig_dir <- file.path("report", "figures")
tab_dir <- file.path("report", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path("analysis", "dataset.csv"), stringsAsFactors = FALSE,
              fileEncoding = "UTF-8")

# A CSV does not store factor levels, so set them again in a fixed order.
d$specialty_group <- factor(d$specialty_group, levels = c("None", "One", "Both"))
d$chain_group     <- factor(d$chain_group, levels = c("Independent", "Chain"))
d$cafe_type       <- factor(d$cafe_type,
                            levels = c("Independent", "Local Chain",
                                       "International Chain"))
d$filter_group    <- factor(d$offers_filter, levels = 0:1,
                            labels = c("No filter", "Offers filter"))
d$roastery_group  <- factor(d$is_roastery, levels = 0:1,
                            labels = c("Other cafe", "Roastery"))

# A3 uses the within-cafe difference: filter price MINUS espresso price.
d$diff_filter_espresso <- d$filter_price - d$espresso_price

cat("Cafes in dataset.csv:", nrow(d), "\n")


# ---- Helper functions --------------------------------------------------------

# Sample skewness (g1): 0 for a symmetric distribution, > 0 = long right tail.
skewness <- function(x) {
  x <- x[!is.na(x)]
  m <- mean(x)
  mean((x - m)^3) / mean((x - m)^2)^1.5
}

# p-value in APA style: three decimals, no leading zero, "< .001" if tiny.
fmt_p <- function(p) {
  ifelse(p < 0.001, "< .001", sub("^0", "", formatC(p, format = "f", digits = 3)))
}

fmt <- function(x, digits = 1) formatC(x, format = "f", digits = digits)

# Numeric summary of one variable: n, missing, mean with t-based 95% CI, SD,
# median, IQR (Q1-Q3) and range.
describe_num <- function(x) {
  n_miss <- sum(is.na(x))
  x  <- x[!is.na(x)]
  n  <- length(x)
  m  <- mean(x)
  se <- sd(x) / sqrt(n)
  tq <- qt(0.975, df = n - 1)
  q  <- quantile(x, c(0.25, 0.75), names = FALSE)
  c(n = n, missing = n_miss, mean = m, ci_low = m - tq * se, ci_high = m + tq * se,
    sd = sd(x), median = median(x), q1 = q[1], q3 = q[2],
    min = min(x), max = max(x))
}

# Write a data frame as a booktabs LaTeX tabular (no caption: the caption is
# added in report.tex around \input{...}). Escapes &, %, _ and # in the text.
write_tex_table <- function(df, file, align = NULL) {
  if (is.null(align)) align <- paste0("l", strrep("r", ncol(df) - 1))
  esc  <- function(s) gsub("([&%_#])", "\\\\\\1", s)
  head <- paste(esc(names(df)), collapse = " & ")
  body <- apply(as.matrix(df), 1, function(r) paste(esc(r), collapse = " & "))
  writeLines(c(paste0("\\begin{tabular}{", align, "}"), "\\toprule",
               paste(head, "\\\\"), "\\midrule", paste(body, "\\\\"),
               "\\bottomrule", "\\end{tabular}"), file)
  cat("Written:", file, "\n")
}

# Histogram with the normal curve of the same mean and SD on top.
hist_normal <- function(x, main, xlab) {
  x <- x[!is.na(x)]
  hist(x, breaks = "FD", freq = FALSE, main = main, xlab = xlab,
       col = "grey85", border = "white")
  grid_x <- seq(min(x), max(x), length.out = 200)
  lines(grid_x, dnorm(grid_x, mean(x), sd(x)), lwd = 2)
}

# Normal Q-Q plot: points close to the line = close to normal.
qq_plot <- function(x, main) {
  x <- x[!is.na(x)]
  qqnorm(x, main = main, pch = 16, col = "grey30")
  qqline(x, lwd = 2)
}

# Boxplot of y by group with the individual cafes as jittered points.
box_points <- function(y, g, main, ylab) {
  boxplot(y ~ g, main = main, ylab = ylab, xlab = "", outline = FALSE,
          col = "grey90")
  stripchart(y ~ g, vertical = TRUE, method = "jitter", add = TRUE, pch = 16,
             col = adjustcolor("grey20", alpha.f = 0.5))
}


# ---- 1-2. Count and summarise ------------------------------------------------

num_vars <- c(
  "Espresso price (CZK)"                      = "espresso_price",
  "Filter coffee price (CZK)"                 = "filter_price",
  "Filter minus espresso, same cafe (CZK)"    = "diff_filter_espresso",
  "Distance from the Orloj (km)"              = "distance_km",
  "East-west distance from the Orloj (km)"    = "ew_km"
)

desc <- t(sapply(num_vars, function(v) describe_num(d[[v]])))
cat("\nNumeric variables:\n"); print(round(desc, 2))

desc_tex <- data.frame(
  Variable  = names(num_vars),
  n         = desc[, "n"],
  Missing   = desc[, "missing"],
  `Mean [95% CI]` = paste0(fmt(desc[, "mean"]), " [", fmt(desc[, "ci_low"]),
                           ", ", fmt(desc[, "ci_high"]), "]"),
  SD        = fmt(desc[, "sd"]),
  Median    = fmt(desc[, "median"]),
  IQR       = paste0(fmt(desc[, "q1"]), "--", fmt(desc[, "q3"])),
  Range     = paste0(fmt(desc[, "min"]), "--", fmt(desc[, "max"])),
  check.names = FALSE)
write_tex_table(desc_tex, file.path(tab_dir, "descriptives.tex"),
                align = "lrrrrrrr")

# Categorical variables: counts and percentages (missing shown separately).
cat_vars <- c("Offers filter coffee" = "filter_group",
              "Roastery"             = "roastery_group",
              "Specialty features"   = "specialty_group",
              "Cafe type"            = "cafe_type",
              "Independent or chain" = "chain_group")
cat_rows <- do.call(rbind, lapply(names(cat_vars), function(lab) {
  x   <- d[[cat_vars[lab]]]
  tab <- table(x, useNA = "always")
  names(tab)[is.na(names(tab))] <- "Missing"
  data.frame(Variable = c(lab, rep("", length(tab) - 1)),
             Category = names(tab),
             n        = as.integer(tab),
             Percent  = fmt(100 * as.integer(tab) / length(x)))
}))
cat("\nCategorical variables:\n"); print(cat_rows, row.names = FALSE)
write_tex_table(cat_rows, file.path(tab_dir, "D_categorical.tex"), align = "llrr")

# A5 contingency table (roastery x filter coffee), shown here for description.
cat("\nRoastery x offers filter coffee:\n")
print(table(d$roastery_group, d$filter_group, useNA = "ifany"))


# ---- 3. Plot every numeric variable -----------------------------------------

for (lab in names(num_vars)) {
  v <- num_vars[lab]
  x <- d[[v]]
  pdf(file.path(fig_dir, paste0("D_", v, ".pdf")), width = 9, height = 3.2)
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  hist_normal(x, main = paste0("Histogram (n = ", sum(!is.na(x)), ")"), xlab = lab)
  boxplot(x, main = "Boxplot", ylab = lab, col = "grey90")
  qq_plot(x, main = "Normal Q-Q plot")
  dev.off()
}

# Groups side by side.
# A2: espresso price, cafes with and without filter coffee
pdf(file.path(fig_dir, "D_espresso_by_filter.pdf"), width = 9, height = 6)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))
box_points(d$espresso_price, d$filter_group, "Espresso price by group",
           "Espresso price (CZK)")
for (g in levels(d$filter_group)) {
  x <- d$espresso_price[d$filter_group %in% g]
  hist_normal(x, main = paste0(g, " (n = ", length(x), ")"),
              xlab = "Espresso price (CZK)")
}
plot.new()   # empty cell to keep the grid aligned
for (g in levels(d$filter_group)) {
  qq_plot(d$espresso_price[d$filter_group %in% g], main = paste("Q-Q:", g))
}
dev.off()

# A4: espresso price by number of specialty features
pdf(file.path(fig_dir, "D_espresso_by_specialty.pdf"), width = 12, height = 6)
par(mfrow = c(2, 4), mar = c(4.5, 4.5, 3, 1))
box_points(d$espresso_price, d$specialty_group,
           "Espresso price by specialty features", "Espresso price (CZK)")
for (g in levels(d$specialty_group)) {
  x <- d$espresso_price[d$specialty_group %in% g]
  hist_normal(x, main = paste0(g, " (n = ", length(x), ")"),
              xlab = "Espresso price (CZK)")
}
plot.new()
for (g in levels(d$specialty_group)) {
  qq_plot(d$espresso_price[d$specialty_group %in% g], main = paste("Q-Q:", g))
}
dev.off()

# A3: the paired comparison is checked on the DIFFERENCES (already plotted as
# D_diff_filter_espresso.pdf). Here: the two prices per cafe, linked by a line.
paired <- d[!is.na(d$filter_price), ]
pdf(file.path(fig_dir, "D_paired_prices.pdf"), width = 5, height = 5)
par(mar = c(4.5, 4.5, 3, 1))
plot(c(1, 2), range(c(paired$espresso_price, paired$filter_price)), type = "n",
     xaxt = "n", xlab = "", ylab = "Price (CZK)", xlim = c(0.7, 2.3),
     main = paste0("Espresso and filter coffee, same cafe (n = ",
                   nrow(paired), ")"))
axis(1, at = 1:2, labels = c("Espresso", "Filter coffee"))
segments(1, paired$espresso_price, 2, paired$filter_price,
         col = adjustcolor("grey30", alpha.f = 0.4))
points(rep(1, nrow(paired)), paired$espresso_price, pch = 16)
points(rep(2, nrow(paired)), paired$filter_price, pch = 16)
dev.off()


# ---- 4-6. Normality, CLT and the decision table ------------------------------

# One row per variable and per group that a question compares. For regressions
# (A7-A10) the assumption is about the residuals, which are checked after the
# model is fitted in 03_analysis.R; the outcome is described here anyway.
checks <- list(
  list("Filter coffee price", "Cafes offering filter", "A1, A6, A7, A9",
       d$filter_price),
  list("Espresso price", "Offers filter", "A2, A6",
       d$espresso_price[d$offers_filter %in% 1]),
  list("Espresso price", "No filter", "A2",
       d$espresso_price[d$offers_filter %in% 0]),
  list("Filter minus espresso", "Cafes offering filter", "A3",
       d$diff_filter_espresso),
  list("Espresso price", "No specialty feature", "A4",
       d$espresso_price[d$specialty_group %in% "None"]),
  list("Espresso price", "One feature", "A4",
       d$espresso_price[d$specialty_group %in% "One"]),
  list("Espresso price", "Both features", "A4",
       d$espresso_price[d$specialty_group %in% "Both"]),
  list("Espresso price", "All cafes", "A8 (outcome)", d$espresso_price)
)

boot_means <- list()   # kept for the bootstrap figure below
decision <- do.call(rbind, lapply(seq_along(checks), function(i) {
  ck <- checks[[i]]
  x  <- ck[[4]][!is.na(ck[[4]])]
  n  <- length(x)
  sk <- skewness(x)
  sw <- shapiro.test(x)$p.value          # optional; does not override plots
  set.seed(SEED + i)                     # reproducible, different per row
  bm <- replicate(B, mean(sample(x, replace = TRUE)))
  boot_means[[i]] <<- bm
  normal <- abs(sk) < SKEW_LIMIT
  clt    <- n >= N_CLT && abs(skewness(bm)) < BOOT_LIMIT
  data.frame(variable = ck[[1]], group = ck[[2]], used_in = ck[[3]], n = n,
             skew = sk, shapiro_p = sw, normal = normal,
             boot_skew = skewness(bm), clt = clt)
}))

# Judgement from the plots. The rule above is only a first screen: look at the
# Q-Q plots and histograms, and if your judgement differs, set it here with a
# comment saying why (row number = row of `decision`). Example:
#   normal_override <- c("1" = FALSE)  # Q-Q plot: long right tail
normal_override <- c()
if (length(normal_override) > 0) {
  decision$normal[as.integer(names(normal_override))] <- normal_override
}

decision$method <- ifelse(decision$normal | decision$clt,
                          "Parametric", "Non-parametric")
cat("\nDecision table:\n")
print(transform(decision, skew = round(skew, 2), boot_skew = round(boot_skew, 2),
                shapiro_p = fmt_p(shapiro_p)), row.names = FALSE)

d_tex <- data.frame(
  Variable          = decision$variable,
  Group             = decision$group,
  `Used in`         = decision$used_in,
  n                 = decision$n,
  Skewness          = fmt(decision$skew, 2),
  `Shapiro--Wilk p` = fmt_p(decision$shapiro_p),
  `Normal?`         = ifelse(decision$normal, "Yes", "No"),
  `CLT applies?`    = ifelse(decision$normal, "--",
                             ifelse(decision$clt, "Yes", "No")),
  Method            = decision$method,
  check.names = FALSE)
write_tex_table(d_tex, file.path(tab_dir, "D_table.tex"), align = "lllrrrccl")


# ---- Figures for Appendix B --------------------------------------------------

row_title <- function(i) {
  paste0(decision$variable[i], "\n", decision$group[i], " (n = ",
         decision$n[i], ")")
}

# Bootstrapped means: if these look normal, the CLT works for that row.
pdf(file.path(fig_dir, "D_bootstrap_means.pdf"), width = 12, height = 6)
par(mfrow = c(2, 4), mar = c(4.5, 4.5, 4, 1))
for (i in seq_along(boot_means)) {
  hist_normal(boot_means[[i]], main = row_title(i), xlab = "Bootstrapped mean (CZK)")
}
dev.off()

# distributions.pdf (Figure 2 in the report template): histogram and Q-Q plot
# of every row of the decision table.
pdf(file.path(fig_dir, "distributions.pdf"), width = 12, height = 12)
par(mfrow = c(4, 4), mar = c(4.5, 4.5, 4, 1))
for (i in seq_along(checks)) {
  x <- checks[[i]][[4]]
  hist_normal(x, main = row_title(i), xlab = "CZK")
  qq_plot(x, main = paste("Q-Q:", decision$group[i]))
}
dev.off()

cat("\nDone. Look at the Q-Q plots in report/figures before trusting the table.\n")
